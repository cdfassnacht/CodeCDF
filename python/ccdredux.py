"""
ccdredux.py - A library of functions to do various basic CCD image processing
              operations

High-level Functions:
   make_bias         - combines input bias or dark frames into one master file
   make_flat         - combines input flat frames into a master flat
   make_flat_files   - exactly like make_flat, but input is a file
                        containing a list of flat-field files rather than
                        an array of frame numbers
   median_combine    - does a generic median combination of a list of images
   apply_calib       - applies calibration corrections to input files
   hdr_offsets       - uses header information in the input files to derive 
                        initial guesses for the offsets between dithered 
                        exposures on a field.
   xcorr_offsets     - does a fft-based cross-correlation to get the initial
                        estimate of the offsets between dithered exposures
                        on a field

Low-level Functions:
   define_trimsec    - sets the section of an image to trim
   sigma_clip        - sigma-clips a data array
"""

try:
   from astropy.io import fits as pf
except:
   import pyfits as pf
import numpy as np
from scipy.ndimage import filters
from math import cos,sin,pi,sqrt,atan2
import imfuncs as imf
import wcs as wcsmwa
import coords

# -----------------------------------------------------------------------

def sigma_clip(data, nsig=3.,verbose=False):
   # Only compute outputs for data values that are numbers
   if verbose:
      print " sigma_clip: Full size of data       = %d" % data.size
      print " sigma_clip: Number of finite values = %d" % \
          data[np.isfinite(data)].size
   d = data[np.isfinite(data)].flatten()
   avg = d.mean()
   std = d.std()
   
   delta = 1
   while delta:
      size = d.size
      d = d[abs(d-avg)<nsig*std]
      avg = d.mean()
      std = d.std()
      delta = size-d.size
   return avg,std

# ---------------------------------------------------------------------------

def robust_sigma(data, refzero=False):
    """
    Calculate a robust estimate of the dispersion of a distribution.
    For an uncontaminated distribution, this estimate is identical to
    the standard deviation.

    This code is ported from robust_sigma.pro in IDL/astrolib

    Inputs:

    Output:
       rsig - robust sigma estimator.  Return -1 if failure.
    """

    """ Set a tolerance """
    eps = 1.e-20

    """ Set central point for cfomputing the dispersion """
    if refzero:
        dat0 = 0.0
    else:
        dat0 = np.median(data)
    datdiff = (data - dat0).flatten()

    """ Find absolute deviation about the median """
    mad = np.median(np.absolute(datdiff))/0.6745

    """ Try the mean absolute deviation if the mad is zero """
    if mad<eps:
        mad = (np.absolute(datdiff)).mean()/0.8
    if mad<eps:
        return 0.0

    """ Do the biweighted value """
    u = datdiff / (6. * mad)
    uu = u*u
    q = uu<1.
    if q.sum()<3:
        print ''
        print 'robust_sigma: input distribution is just too weird.'
        print 'returning value of -1.'
        return -1.
    ntot = data[np.isfinite(data)].sum()
    num = ((data[q] - dat0)**2 * (1.-uu[q])**4).sum()
    denom = ((1.-uu[q]) * (1. - 5.*uu[q])).sum()
    rvar = ntot * num / (denom * (denom - 1.))

    if rvar>0.:
        return np.sqrt(rvar)
    else:
        return 0.

# -----------------------------------------------------------------------

def set_param_array(hdulen,inval):
   """
   Converts an input parameter value, which may have been passed as an
   integer or float, into an array so that it can be used in general 
   image-processing functions.
   
   Inputs:
      hdulen - the number of HDUs in the associated fits file
      inval  - the value of the variable

   Output:
      valarr - the input value as an array, if necessary
   """

   if (isinstance(inval,int)) or (isinstance(inval,float)):
      if hdulen == 1:
         valarr = np.array([inval]).astype(int)
      else:
         valarr = np.zeros(hdulen-1).astype(int)
         for i in range(hdulen-1):
            valarr[i] = inval
   else:
      valarr = inval

   return valarr

# -----------------------------------------------------------------------

def define_trimsec(hdu,x1,x2,y1,y2):
   xmax = hdu.header['NAXIS1']
   ymax = hdu.header['NAXIS2']
   #xmax = pf.getval(fullfits,'NAXIS1')
   #ymax = pf.getval(fullfits,'NAXIS2')
   if x1 != 0:
      x1 = int(x1)
   if x2 != 0:
      x2 = int(x2)
   else:
      x2 = xmax
   #xtrim = slice(x1,x2,1)
   if y1 != 0:
      y1 = int(y1)
   if y2 != 0:
      y2 = int(y2)
   else:
      y2 = ymax
   #ytrim = slice(y1,y2,1)
   return x1,x2,y1,y2

# -----------------------------------------------------------------------

def divide_images(a,b,output,preserve_header=0):
   print "Dividing images: '%s' / '%s' = '%s'" % (a,b,output)
   try:
      hdua = pf.open(a)[0]
   except:
      hdua = pf.open(a,ignore_missing_end=True)[0]
   try:
      hdub = pf.open(b)[0]
   except:
      hdub = pf.open(b,ignore_missing_end=True)[0]
   if preserve_header == 1:
      hdu = pf.PrimaryHDU(hdua.data / hdub.data,header=hdua.header)
   elif preserve_header == 2:
      hdu = pf.PrimaryHDU(hdua.data / hdub.data,header=hdub.header)
   else:
      hdu = pf.PrimaryHDU(hdua.data/hdub.data)
   hdu.writeto(output,output_verify='ignore')

# -----------------------------------------------------------------------

def subtract_images(a,b,output,hexta=0,hextb=0):
   print("Subtracting images: '%s' - '%s' = '%s'" % (a,b,output))
   try:
      hdua = pf.open(a)[hexta]
   except:
      hdua = pf.open(a,ignore_missing_end=True)[hexta]
   try:
      hdub = pf.open(b)[hextb]
   except:
      hdub = pf.open(b,ignore_missing_end=True)[hextb]
   hdu = pf.PrimaryHDU(hdua.data - hdub.data)
   hdu.writeto(output,output_verify='ignore',clobber=True)

# -----------------------------------------------------------------------

def add_images(a,b,output,preserve_header=0):
   print("Adding images: '%s' + '%s' = '%s'" % (a,b,output))
   try:
      hdua = pf.open(a)[0]
   except:
      hdua = pf.open(a,ignore_missing_end=True)[0]
   try:
      hdub = pf.open(b)[0]
   except:
      hdub = pf.open(b,ignore_missing_end=True)[0]

   if preserve_header == 1:
      hdu = pf.PrimaryHDU(hdua.data + hdub.data,header=hdua.header)
   elif preserve_header == 2:
      hdu = pf.PrimaryHDU(hdua.data + hdub.data,header=hdub.header)
   else:
      hdu = pf.PrimaryHDU(hdua.data + hdub.data)

   hdu.writeto(output,output_verify='ignore',clobber=True)

# -----------------------------------------------------------------------

def read_calfile(filename, file_description):
   print 'Reading in %s file: %s' % (file_description, filename)
   try:
      calhdulist = pf.open(filename)
   except:
      try:
         calhdulist = pf.open(filename,ignore_missing_end=True)
      except:
         print " ERROR: Requested %s file %s does not exist" % \
             (file_description, filename)
         print ""
         return -1
   return calhdulist

# -----------------------------------------------------------------------

def median_combine(input_files, output_file, method='median', x1=0, x2=0,
                   y1=0, y2=0, biasfile=None, gain=-1.0, normalize=False,
                   zeromedian=False, NaNmask=False, hdu0only=False):
   """ 
   Given a list of input file names, this function will:
      1. Subtract a bias frame (if the optional biasfile parameter is set)
      2. Multiply by the gain, required to be in e-/ADU (if the optional
         gain parameter is set)
      3. Normalize the frame (if the optional normalize parameter is set)
      4. Subtract the median (if the optional zeromedian parameter is set)
      5. Median combine the resulting data
      6. Write the output to a file.
   """

   print "median_combine: Inputs:"
   print "-----------------------"
   if (biasfile != None):
      # NB: Should not need to trim the input bias frame, since it should
      #  have been created out of trimmed files.
      bias = pf.open(biasfile)
      bias.info()
   else:
      print "  bias frame: [No bias file]"

   print ""
   print "median_combine: Loading files"
   print "-----------------------------"
   files = []
   for filename in input_files:
      print " %s" % filename
      try:
         f = pf.open(filename)
      except:
         f = pf.open(filename, ignore_missing_end=True)
      files.append(f)

   # Use first file to check for multiple HDUs, and set some variables
   #  depending on the result
   print ""
   print "median_combine: Getting info on first file"
   print "------------------------------------------"
   files[0].info()
   hdulen = len(files[0])
   if (hdulen == 1) or (hdu0only == True):
      imhdu = np.arange(1)
   else:
      imhdu = np.arange(1, hdulen)
      phdu = pf.PrimaryHDU()
      phdu.header.add_comment(
         'This file contains %d image extensions' % (hdulen-1), after='extend')
      hdulist = pf.HDUList([phdu])

   # Set up the trim variables as numpy arrays, if they are not
   #  already in that format

   x1 = set_param_array(hdulen, x1)
   x2 = set_param_array(hdulen, x2)
   y1 = set_param_array(hdulen, y1)
   y2 = set_param_array(hdulen, y2)
   gain = set_param_array(hdulen, gain)

   # Start outer loop on the extension number.  For old-school FITS files
   #  there is only one extension, which is the main image.

   for j in imhdu:

      # Set trim section for this HDU in the input files, and use that
      #  to define the container for the stack of images

      k = j-1
      xt1, xt2, yt1, yt2 = define_trimsec((files[0])[j], x1[k], x2[k], y1[k],
                                          y2[k])
      print ""
      print "median_combine: setting up stack for images (HDU %d)" % j
      print "----------------------------------------------------"
      print "Stack will have dimensions (%d, %d, %d)" \
          %(len(input_files), yt2 - yt1, xt2 - xt1)
      stack = np.zeros((len(input_files), yt2 - yt1, xt2 - xt1))

      # Inner loop is on the input files

      count = 0
      for i in range(len(input_files)):
         
         print " %s" % files[i].filename()

         # Process the data (bias and gain only), if desired
         if (biasfile == None):
            process_data(files[i], j, gain=gain[k], x1=x1[k], x2=x2[k], \
                            y1=y1[k], y2=y2[k])
         else:
            process_data(files[i], j, bias, gain=gain[k], x1=x1[k], \
                            x2=x2[k], y1=y1[k], y2=y2[k])

         tmpf = (files[i])[j].data

         # Normalize or set to zero median, if desired
         if(normalize == True):
            frame_med = np.median(tmpf,axis=None)
            print "    Normalizing %s by %f" % (files[i].filename(),frame_med)
            tmpf *= 1.
            tmpf /= (1. * frame_med)
         if(zeromedian == True):
            print "    Subtracting the median from %s" % files[i].filename()
            tmpf -= np.median(tmpf,axis=None)
         stack[count] = tmpf.copy()
         count += 1
      
      print ""

      # Actually form the median (or sum, if that was requested)
      if method == 'sum':
         if(NaNmask == True):
            print "median_combine: Computing summed frame using NaN masking"
            print "   Can take a while..."
            outdat = np.nansum(stack, axis=0)
         else:
            print "median_combine: Computing summed frame (can take a while)..."
            outdat = np.sum(stack, axis=0)
      else:
         if(NaNmask == True):
            print "median_combine: Computing median frame using NaN masking"
            print "   Can take a while..."
            outdat = np.nanmedian(stack, axis=0)
         else:
            print "median_combine: Computing median frame (can take a while)..."
            outdat = np.median(stack, axis=0)
      del stack

      # Save the median HDU
      # For multi-extension files, use pf.ImageHDU(outdat)
      if(hdulen == 1) or (hdu0only):
         phdu = pf.PrimaryHDU(outdat)
         hdulist = pf.HDUList([phdu])
      else:
         hdu = pf.ImageHDU(outdat)
         hdulist.append(hdu)

   # Write the output median file

   hdulist.writeto(output_file, output_verify='ignore', clobber=True)
   print "   ... Writing output to %s." % output_file

   # Clean up

   for i in range(len(input_files)):
      files[i].close()

# -----------------------------------------------------------------------

def apply_rough_wcs(hdu, pixscale, rakey='ra', deckey='dec', phdu=None):
   """
   Takes the RA and Dec pointing info from the fits header, along with
   a pixel scale, and converts that information into the standard WCS
   headers (CD1_1, etc.).

   For multi-extension FITS files may contain the RA and Dec information
   only in the PHDU and not in the image extension HDUs.  In this case,
   the function can be passed the PHDU via the optional phdu parameter.
   If phdu=None then the function will look for the RA and Dec header
   cards in the hdu that has been passed via the hdu parameter.
   """

   """ Set up phdu in case there is not a separate one """
   if phdu is None:
      phdu = hdu

   if pixscale > 0.0:

      """ 
      Get the size of the data array.  Note that if this fails, the
      probable reason is that the selected hdu is the PHDU of a
      multiextension fits file.  In that case, don't set the CRPIXes
      """

      containsdata = True
      try:
         shape = hdu.data.shape
      except:
         containsdata = False

      """
      Read the header cards
      """
      hdr = hdu.header
      phdr = phdu.header

      wcsread = True
      if wcsmwa.is_degree(hdr[rakey]):
         ra = hdr[rakey]
      else:
         try:
            ra  = wcsmwa.ra2deg(hdr[rakey].strip())
         except:
            try:
               ra  = wcsmwa.ra2deg(phdr[rakey].strip())
            except:
               print 'ERROR. Attempts to read RA header card (%s) failed.' % \
                   rakey.upper()
               print 'No wcs information'
               wcsread = False

      if wcsmwa.is_degree(hdr[deckey]):
         dec = hdr[deckey]
      else:
         try:
            dec = wcsmwa.dec2deg(hdr[deckey].strip())
         except:
            try:
               dec = wcsmwa.dec2deg(phdr[deckey].strip())
            except:
               print 'ERROR. Attempts to read Dec header card failed.'
               print 'No wcs information'
               wcsread = False

      """ Clean out old WCS and mosaic info """
      #foo = hdr.ascardlist()[4:]
      foo = hdr.keys()[4:]
      for k in range(0,len(foo)):
         #keyname = foo[k].key
         keyname = foo[k]
         if keyname[0:4] == 'CD1_' or keyname[0:4] == 'CD2_' or \
                keyname[0:4] == 'PANE' or keyname[0:5] == 'CRVAL' \
                or keyname[0:5] == 'CRPIX' or keyname[0:5] == 'CDELT' \
                or keyname[0:5] == 'CTYPE' or keyname[0:5] == 'CUNIT' \
                or keyname[0:5] == 'CRDER' or keyname[0:5] == 'CSYER' \
                or keyname[0:7] == 'WCSNAME' or keyname[0:3] == 'ADC':
            del hdr[keyname]

      """ Apply the rough WCS info """
      if wcsread:
         hdr.update('CRVAL1',ra)
         hdr.update('CRVAL2',dec)
         if containsdata:
            hdr.update('CRPIX1',shape[1]/2.)
            hdr.update('CRPIX2',shape[0]/2.)
         hdr.update('CD1_1',-pixscale/3600.)
         hdr.update('CD1_2',0.)
         hdr.update('CD2_1',0.)
         hdr.update('CD2_2',pixscale/3600.)
         hdr.update('CTYPE1','RA---TAN')
         hdr.update('CTYPE2','DEC--TAN')
         hdr.update('EQUINOX',2000.0)
         hdr.update('RADECSYS','FK5')

# -----------------------------------------------------------------------

def make_bias(infiles, outfile="Bias.fits", x1=0, x2=0, y1=0, y2=0, 
              hdu0only=False):
   """ 

       This function takes an input list of dark frames (either true darks, 
       or bias frames) and median-combines them to create a master dark/bias

       Required inputs:
        infiles     - a list of input files

       Optional inputs:
        outfile     - output filename (default="masterdark.fits")
        x1          - to set a trim section that is smaller than the full frame
        x2          - to set a trim section that is smaller than the full frame
        y1          - to set a trim section that is smaller than the full frame
        y2          - to set a trim section that is smaller than the full frame

   """

   median_combine(infiles,outfile,x1,x2,y1,y2,hdu0only=hdu0only)

# -----------------------------------------------------------------------

def make_bias_frames(bias_frames, raw_prefix, rawdir="../Raw", rawext='.fits',
                     outfile="Bias.fits", x1=0,x2=0,y1=0,y2=0,hdu0only=False):
   """ 

       This function takes as input the frame numbers of the
       dark frames (either true darks, or bias frames) and
       median-combines them to create a master dark

       Required inputs:
        dark_frames - an array of frame numbers (e.g., [101,102,103,105])
        raw_prefix  - prefix for frame numbers (e.g., "lred0")

       Optional inputs:
        rawdir      - directory with raw files (default="../Raw")
        rawext      - extension for raw filenames (default = 'fits')
        outfile     - output filename (default="masterdark.fits")
        x1          - to set a trim section that is smaller than the full frame
        x2          - to set a trim section that is smaller than the full frame
        y1          - to set a trim section that is smaller than the full frame
        y2          - to set a trim section that is smaller than the full frame

   """

   # Make file list
   filenames = []
   for i in bias_frames:
      filenames.append('%s/%s%d.%s'%(rawdir,raw_prefix,i,rawext))

   # Call median_combine
   median_combine(filenames,outfile,x1,x2,y1,y2,hdu0only=hdu0only)

# -----------------------------------------------------------------------

def make_flat(flat_frames, raw_prefix, rawdir="../Raw", rawext='.fits',
      outfile="Flat.fits", biasfile=None, gain=-1.0, normalize=True,
      x1=0, x2=0, y1=0, y2=0, framesig=0):
   """ 

       This function takes as input the frame numbers of the flat-field
       frames and median-combines them to create a master flat

       Required inputs:
        flat_frames - array of input frame numbers (e.g., [101,102,103,105])
        raw_prefix  - prefix for frame numbers (e.g., "lred0")

       Optional inputs:
        rawdir      - directory with raw files (default="../Raw")
        rawext      - extension for raw filenames (default = '.fits')
        outfile     - output filename (default="Flat.fits")
        biasfile    - input bias file to subtract before combining. 
                      (default=None)
        gain        - gain to divide by before combining (default=1.0)
        normalize   - normalize flat-field frames before combining?
                      (default = True)
        x1          - to set a trim section that is smaller than the full frame
        x2          - to set a trim section that is smaller than the full frame
        y1          - to set a trim section that is smaller than the full frame
        y2          - to set a trim section that is smaller than the full frame
        framesig    - number of significant digits in the frame numbers.  Keep
                      at the default (framesig=0) if all of the frame numbers
                      have the same number of significant digits (e.g., if all of
                      them are between 100 and 999).  However, if there is a
                      range of significant figures (e.g., frame numbers between
                      75 and 125), then set this to the higher number of figures.

   """

   # Make file list
   filenames = []
   for i in flat_frames:
      if framesig == 2:
         filenames.append('%s/%s%02d%s'%(rawdir,raw_prefix,i,rawext))
      elif framesig == 3:
         filenames.append('%s/%s%03d%s'%(rawdir,raw_prefix,i,rawext))
      elif framesig == 4:
         filenames.append('%s/%s%04d%s'%(rawdir,raw_prefix,i,rawext))
      else:
         filenames.append('%s/%s%d%s'%(rawdir,raw_prefix,i,rawext))

   # Call median_combine
   median_combine(filenames,outfile,biasfile=biasfile,gain=gain,
                  normalize=normalize,x1=x1,x2=x2,y1=y1,y2=y2)

# -----------------------------------------------------------------------

def make_flat_files(infiles, outfile="Flat.fits", biasfile=None, gain=1.0, 
                    normalize=True, x1=0, x2=0, y1=0, y2=0):
   """ 

       This function takes as input the frame numbers of the flat-field
       frames and median-combines them to create a master flat

       Required inputs:
        infiles - list of the input files that contain the individual
                  flat-field exposures.

       Optional inputs:
        outfile     - output filename (default="Flat.fits")
        biasfile    - input bias file to subtract before combining. 
                      (default=None)
        gain        - gain to divide by before combining (default=1.0)
        normalize   - normalize flat-field frames before combining?
                      (default = True)
        x1          - to set a trim section that is smaller than the full frame
        x2          - to set a trim section that is smaller than the full frame
        y1          - to set a trim section that is smaller than the full frame
        y2          - to set a trim section that is smaller than the full frame

   """

   """  Call median_combine """
   median_combine(infiles,outfile,biasfile=biasfile,gain=gain,
                  normalize=normalize,x1=x1,x2=x2,y1=y1,y2=y2)

# -----------------------------------------------------------------------

def make_fringe(fringe_frames, in_prefix, indir=None, inext='fits',
      outfile="Fringe.fits", biasfile=None, gain=1.0, normalize=False,
      zeromedian=True, x1=0, x2=0, y1=0, y2=0):
   """ 

       This function takes an array containing the flat-fielded frame 
       numbers and subtracts from each flat-fielded science frame
       its individual median.  The resulting data are median-combined
        to create a master fringe file

       Required inputs:
        fringe_frames - array of input frame numbers (e.g., [101,102,103,105])
        in_prefix     - prefix for frame numbers (e.g., "ff")

       Optional inputs:
        indir       - directory containing input files (default=None)
        inext       - extension for input filenames (default = 'fits')
        outfile     - output filename (default="Fringe.fits")
        biasfile    - input bias file to subtract before combining. 
                      (default=None)
        gain        - gain to divide by before combining (default=1.0)
        zeromedian  - subtract median from each science frame before combining?
                      (default = True)
        normalize   - normalize science frames before combining?
                      (default = False)
        x1          - to set a trim section that is smaller than the full frame
        x2          - to set a trim section that is smaller than the full frame
        y1          - to set a trim section that is smaller than the full frame
        y2          - to set a trim section that is smaller than the full frame

   """

   # Make file list
   if indir is None:
      indir = '.'
   filenames = []
   for i in fringe_frames:
      filenames.append('%s/%s%d.%s'%(indir,in_prefix,i,inext))

   # Call median_combine
   median_combine(filenames,outfile,biasfile=biasfile,gain=gain,
    normalize=normalize,zeromedian=zeromedian,x1=x1,x2=x2,y1=y1,y2=y2)

# -----------------------------------------------------------------------

def make_fringe_files(fringe_frames, in_prefix, indir=None, 
      outfile="Fringe.fits", biasfile=None, gain=1.0, normalize=False,
      zeromedian=True, x1=0, x2=0, y1=0, y2=0):
   """ 

       This function takes a list of the flat-fielded input files
       and subtracts their individual medians before median-combining
       them to create a master fringe file

       Required inputs:
        file_with_filelist - file containing list of input files 
          (e.g., "filelist.txt").  The input files should be
          flat-fielded science files.

       Optional inputs:
        indir       - directory containing input files (default=None)
        outfile     - output filename (default="Flat.fits")
        biasfile    - input bias file to subtract before combining. 
                      (default=None)
        gain        - gain to divide by before combining (default=1.0)
        zeromedian  - subtract median from each science frame before combining?
                      (default = True)
        normalize   - normalize science frames before combining?
                      (default = False)
        x1          - to set a trim section that is smaller than the full frame
        x2          - to set a trim section that is smaller than the full frame
        y1          - to set a trim section that is smaller than the full frame
        y2          - to set a trim section that is smaller than the full frame

   """

   # Extract base file names from input file
   basenames = open(file_with_filelist).read().split()

   # Make file list
   if indir is None:
      filenames = basenames
   else:
      filenames = []
      for i in basenames:
         filenames.append('%s/%s'%(indir,i))

   # Call median_combine
   median_combine(filenames,outfile,biasfile=biasfile,gain=gain,
    normalize=normalize,zeromedian=zeromedian,x1=x1,x2=x2,y1=y1,y2=y2)

# -----------------------------------------------------------------------

def process_data(hdulist, hdunum, bias=None, flat=None, fringe=None, 
                 darksky=None, gain=-1.0, texp_key=None, skysub=False,
                 flip=0, pixscale=0.0, rakey='ra', deckey='dec',
                 x1=0, x2=0, y1=0, y2=0):

   """ This function applies calibration corrections to the passed HDU.  All
        of the calbration steps are by default turned off (keywords set to None).
        To apply a particular calibration step, set the appropriate keyword.
        The possible steps, along with their keywords are:

          Keyword     Calibration step
          ----------  ----------------------------------
          bias        Bias subtraction
          gain        Convert from ADU to electrons if set to value > 0
                      NB: Gain must be in e-/ADU
          flat        Flat-field correction
          fringe      Fringe subtraction
          darksky     Dark-sky flat correction
          skysub      Subtract mean sky level if keyword set to True
          texp_key    Divide by exposure time (set keyword to fits header name)
          flip        0 => no flip
                      1 => PFCam style (flip x then rotate -90), 
                      2 => P60 CCD13 style (not yet implemented)
                      3 => flip x-axis
          pixscale    If >0, apply a rough WCS using this pixel scale (RA and
                       Dec come from telescope pointing info in fits header)
          rakey       FITS header keyword for RA of telescope pointing.
                      Default = 'ra'
          deckey      FITS header keyword for Dec of telescope pointing.
                      Default = 'dec'
      
       Required inputs:
        in_frames   - list of input file names
        in_prefix   - input prefix that will be replaced by the output prefix
        out_prefix  - prefix for the flattened output frames

       Optional inputs (in addition to those listed above in the keyword list):
        rawdir      - directory with raw files (default="../Raw")
        x1          - to set a trim section that is smaller than the full frame
        x2          - to set a trim section that is smaller than the full frame
        y1          - to set a trim section that is smaller than the full frame
        y2          - to set a trim section that is smaller than the full frame
   """

   # Make the name of the passed HDU more tractable
   tmp = hdulist[hdunum]
   phdu = hdulist[0]

   # Trim the data if requested
   xt1,xt2,yt1,yt2 = define_trimsec(tmp,x1,x2,y1,y2)
   tmp.data = tmp.data[yt1:yt2,xt1:xt2].astype(float)
   if xt2-xt1 != tmp.header['naxis1'] or yt2-yt1 != tmp.header['naxis2']:
      tmp.header['trim'] = \
          'Trim data section is [%d:%d,%d:%d] ([xrange,yrange])' % \
          (xt1,xt2,yt1,yt2)
      print "   Trimmed data to section [xrange,yrange] [%d:%d,%d:%d]" \
          % (xt1,xt2,yt1,yt2)
   
   # Set up a string for use in header keywords
   if hdunum==0:
      hdustr = 'these data'
   else:
      hdustr = 'HDU %d' % hdunum
   
   # Bias-subtract if requested
   if bias is not None:
      biasdata = bias[hdunum].data
      tmp.data -=  biasdata
      biasmean = biasdata.mean()
      if (hdunum == 0):
         keystr = 'biassub'
      else:
         keystr = 'biassub'+str(hdunum)
      tmp.header[keystr] = 'Bias frame for %s is %s with mean %f' % \
          (hdustr,bias.filename(),biasmean)
      print "   Subtracted bias frame %s" % bias.filename()
   
   # Convert to electrons if requested
   if gain>0:
      tmp.data *= gain
      tmp.header['gain'] =  (1.0, 'Units are now electrons')
      if (hdunum == 0):
         keystr = 'gainorig'
      else:
         keystr = 'gainori'+str(hdunum)
      tmp.header.set(keystr, gain, 'Original gain for %s in e-/ADU' % hdustr,
                            after='gain')
      tmp.header['bunit'] =  ('Electrons','Converted from ADU in raw image')
      if (hdunum == 0):
         keystrb1 = 'binfo_1'
      else:
         keystrb1 = 'binfo'+str(hdunum)+'_1'
      keystrb1 = keystrb1.upper()
      tmp.header[keystrb1] = \
          'Units for %s changed from ADU to e- using gain=%6.3f e-/ADU' % \
          (hdustr,gain)
      print "   Converted units to e- using gain = %f" % gain
   
   # Divide by the exposure time if requested
   if texp_key is not None:
      texp_good = True
      try:
         texp = tmp.header[texp_key]
      except:
         try:
            texp = phdu.header[texp_key]
         except:
            print("")
            print("ERROR: No exposure time keyword called %s found in header"
                  % texp_key)
            print(" Setting texp = 1.0")
            texp_good = False
            texp = 1.0
      tmp.data /= texp
      if texp_good:
         if (hdunum == 0):
            keystr = 'binfo_2'
         else:
            keystr = 'binfo'+str(hdunum)+'_2'
         keystr = keystr.upper()
         tmp.header['gain'] = (texp, 'If units are e-/s then gain=t_exp')
         tmp.header['bunit'] = ('Electrons/sec','See %s header' % keystr,
                                keystr)
         tmp.header.set(keystr,
                        'Units for %s changed from e- to e-/s using texp=%7.2f'
                        % (hdustr,texp), after=keystrb1)
         print "   Converted units from e- to e-/sec using exposure time %7.2f" \
             % texp
      else:
         tmp.header['bunit'] = \
             ('Electrons', 'There was an ERROR in converting to e-/s')
         tmp.header[keystr] = 'ERROR: Exposure time query failed.'
   
   # Apply the flat-field correction if requested
   if flat is not None:
      flatdata = flat[hdunum].data
      tmp.data /= flatdata
      flatmean = flatdata.mean()
      # Set up a bad pixel mask based on places where the flat frame = 0,
      #  since dividing by zero gives lots of problems
      zeromask = flatdata==0
      # Correct for any zero pixels in the flat-field frame
      tmp.data[zeromask] = 0
      if (hdunum == 0):
         keystr = 'flatcor'
      else:
         keystr = 'flatcor'+str(hdunum)
      tmp.header[keystr] = \
          'Flat field image for %s is %s with mean=%f' % \
          (hdustr,flat.filename(),flatmean)
      print "   Divided by flat-field image: %s" % flat.filename()
   
   # Apply the fringe correction if requested
   if fringe is not None:
      tmp.data -= fringe
      if (hdunum == 0):
         keystr = 'fringcor'
      else:
         keystr = 'frngcor'+str(hdunum)
      tmp.header[keystr] = \
          'Fringe image for %s is %s with mean=%f' % \
          (hdustr,fringe.filename(),fringemean)
      print "   Subtracted fringe image: %s" % fringe.filename()
   
   # Apply the dark sky flat-field correction if requested
   if darksky is not None:
      tmp.data /= darksky
      # Correct for any zero pixels in the flat-field frame
      tmp.data[dszeromask] = 0
      if (hdunum == 0):
         keystr = 'darksky'
      else:
         keystr = 'darksky'+str(hdunum)
      tmp.header[keystr] = \
          'Dark-sky flat image for %s is %s with mean=%f' % \
          (hdustr,darksky.filename(),darkskymean)
      print "   Divided by dark-sky flat: %s" % darksky.filename()

   # Subtract the sky level if requested
   if skysub:
      m,s = sigma_clip(tmp.data)
      tmp.data -= m
      if (hdunum == 0):
         keystr = 'skysub'
      else:
         keystr = 'skysub'+str(hdunum)
      tmp.header[keystr] = ('For %s, subtracted mean sky level of %f' % \
                               (hdustr,m))
      print '   Subtracted mean sky level of %f' % m
   
   # Flip if requested
   if flip == 1:
      tmp.data = tmp.data.T[::-1,::-1]
   if flip == 3:
      tmp.data = tmp.data[:,::-1]
   
   # Add a very rough WCS if requested
   if pixscale > 0.0:
      if hdunum>0:
         apply_rough_wcs(tmp,pixscale,rakey,deckey,hdulist[0])
      else:
         apply_rough_wcs(tmp,pixscale,rakey,deckey)

#   return hdu

# -----------------------------------------------------------------------

def apply_calib(in_frames, in_prefix, out_prefix, split=False,
                biasfile=None, flatfile=None, fringefile=None, darkskyfile=None,
                skysub=False, gain=-1.0, texp_key=None, 
                flip=0, pixscale=0.0, rakey='ra', deckey='dec',
                rawdir="../Raw", rawext='.fits', x1=0, x2=0, y1=0, y2=0):
   """ This function applies calibration corrections to the input files,
        which are designated by an array of frame numbers.  All of the
        calbration steps are by default turned off (keywords set to None).
        To apply a particular calibration step, set the appropriate keyword.
        The possible steps, along with their keywords are:

          Keyword     Calibration step
          ----------  ----------------------------------
          biasfile    Bias subtraction
          gain        Convert from ADU to electrons if set to value > 0
                      NB: gain must be in e-/ADU
          flatfile    Flat-field correction
          fringefile  Fringe subtraction
          darkskyfile Dark-sky flat correction
          skysub      Subtract mean sky level if keyword set to True
          texp_key    Divide by exposure time (set keyword to fits header name)
          flip        0 => no flip
                      1 => PFCam-style (flip x then rotate -90), 
                      2 => P60 CCD13 style (not yet implemented)
                      3 => flip x-axis
          pixscale    If >0, apply a rough WCS using this pixel scale (RA and
                       Dec come from telescope pointing info in fits header)
          rakey       FITS header keyword for RA of telescope pointing.
                      Default = 'ra'
          deckey      FITS header keyword for Dec of telescope pointing.
                      Default = 'dec'
      
       Required inputs:
        in_frames   - list of input file names
        in_prefix   - input prefix that will be replaced by the output prefix
        out_prefix  - prefix for the flattened output frames

       Optional inputs (in addition to those listed above in the keyword list):
        rawdir      - directory with raw files (default="../Raw")
        rawext      - extension for raw filenames (default = '.fits')
        x1          - to set a trim section that is smaller than the full frame
        x2          - to set a trim section that is smaller than the full frame
        y1          - to set a trim section that is smaller than the full frame
        y2          - to set a trim section that is smaller than the full frame
   """

   """ Read in calibration frames if they have been selected """

   print ''
   if biasfile is not None:
      bias = read_calfile(biasfile,'bias')
      if (bias == -1):
         return
   else:
      bias = None

   if flatfile is not None:
      print 'Reading in flat-field file: %s' % flatfile
      flat = read_calfile(flatfile,'flat-field')
      if (flat == -1):
         return
   else:
      flat = None

   if fringefile is not None:
      print 'Reading in fringe file: %s' % fringefile
      fringe = pf.getdata(fringefile)
      fringemean = fringe.mean()
   else:
      fringe = None

   if darkskyfile is not None:
      print 'Reading in dark-sky flat file: %s' % darkskyfile
      darksky = pf.getdata(darkskyfile)
      darkskymean = darksky.mean()
      """
      Set up a bad pixel mask based on places where the flat frame = 0,
       since dividing by zero gives lots of problems
      """
      dszeromask = darksky==0
   else:
      darksky = None

   print ""
   print "Processing files..."
   print "-------------------"

   write_one_output_file = True
   for i in in_frames:
      filename = '%s/%s%s%s'%(rawdir,in_prefix,i,rawext)
      print "%s:" % filename

      """
      Open the input file.  
      """
      try:
         hdu = pf.open(filename)
      except:
         hdu = pf.open(filename,ignore_missing_end=True)

      """ Set things up depending on whether there are extensions or not """
      hdulen = len(hdu)
      if hdulen == 1:
         imhdu = np.arange(1)
      else:
         imhdu = np.arange(1,hdulen)
         phdu = hdu[0]
         phdu.header.add_comment(
            'This file contains %d image extensions' % (hdulen-1),after='extend')
         outhdu = pf.HDUList([phdu])

      """
      Set up the trim and flip variables as numpy arrays, if they are not
       already in that format
      """

      x1 = set_param_array(hdulen,x1)
      x2 = set_param_array(hdulen,x2)
      y1 = set_param_array(hdulen,y1)
      y2 = set_param_array(hdulen,y2)
      gain = set_param_array(hdulen,gain)
      flip = set_param_array(hdulen,flip)

      """
       Loop over fits extensions.  This loop assumes that if the input file
        is a multi-extension fits file, then HDU 0 does not contain any
        data, and only HDU's 1-N do.  Of course, if this is a straightforward
        fits file with no extensions, then HDU 0 will contain the data
      """

      for j in imhdu:

         # Check that extension is valid (to be implemented)

         # Read in the data
         if hdulen>0:
            print " Processing image extension %d" % j

         # Process the data
         k = j-1
         process_data(hdu,j,bias,flat,fringe,darksky,gain[k],texp_key,
                      skysub,flip[k],pixscale,rakey,deckey,
                      x1[k],x2[k],y1[k],y2[k])

         # If multiple extensions, append the current extension to the HDU list
         if hdulen>1:
            if(split):
               write_one_output_file = False
               print " Splitting image extension %d to output file."
            else:
               outhdu.append(hdu[j])
         else:
            outhdu = pf.PrimaryHDU(hdu[j].data,hdu[j].header)

      # Write out final file
      # if hdulen>1:
      #    newhdu = pf.PrimaryHDU(tmp.data)
      #    foo = phdu.header.ascardlist()[4:]
      #    for k in range(0,len(foo)):
      #       if foo[k].key[0:5] != 'NAXIS' and foo[k].key[0:4] != 'PANE' and \
      #              foo[k].key != 'COMMENT':
      #          newhdu.header.update(foo[k].key,foo[k].value,foo[k].comment)
      # else:
      #    newhdu = tmp

      if(write_one_output_file):
         outname = '%s%s.fits'%(out_prefix,i)
         outname = outname.strip()
         outhdu.writeto(outname,output_verify='ignore',clobber=True)
         print " Writing output file: %s" % outname

      hdu.close()

# For LRIS B
#  x1 = [400, 51, 51, 400]
#  x2 = 1068
#  y1 = 775
#  y2 = 3200

# -----------------------------------------------------------------------

def make_wht_from_pixval(infile, maxgood, inwhtfile=None, outsuff='_wht',
                         goodval=1, hdunum=0):
   """
   Creates a weight file that is essentially a bad pixel mask, where the
   bad pixels are defined by their pixel values.  This function can be used,
   for example, to flag cosmic rays and saturated pixels.

   The user can either pass an existing weight file, in which case that file
   will be modified to incorporate the new information, or can create a
   new weight file from scratch.

   Inputs:
      infile    - input science file
      maxgood   - maximum good pixel value. Pixels with values greater than
                  maxgood are flagged as bad.
      inwhtfile - optional input weight file.  This file will get modified to
                  include information about the bad pixels by this function.
                  If inwhtfile is None (the default), then a new weight file
                  will be created from scratch.
      outsuff   - suffix for brand new weight file if inwhtfile is None.
                  Default = '_wht', so if infile is foo.fits the new weight
                  file will be foo_wht.fits
      goodval   - Pixel value in the weight file indicating a good pixel.
                  This can either be goodval=1 (the default), in which case
                  bad pixels will be marked with 0, or goodval=0, in which
                  case bad pixels will be marked with 1 (as in a bad pixel mask)
   """

   """ Get science data """
   try:
      indat = pf.getdata(infile,hdunum)
   except:
      print ''
      print 'ERROR: Could not open input file %s' % infile
      print ''
      return

   """ 
   Open input weight file, if one has been requested.  If not, then create
   an array to be used as the basis for the new weight file.
   """
   if inwhtfile is not None:
      try:
         whthdu = pf.open(inwhtfile,mode='update')
      except:
         del indat
         print ''
         print 'ERROR. Could not open input weight file %s' % inwhtfile
         print ''
         return
      whtdat = whthdu[hdunum].data
   else:
      if goodval == 1:
         whtdat = np.ones(indat.shape)
      else:
         whtdat = np.zeros(indat.shape)

   """ Flag the bad pixels """
   if goodval == 1:
      whtdat[indat>maxgood] = 0
   else:
      whtdat[indat>maxgood] = 1

   """ Save the result """
   if inwhtfile is not None:
      whthdu.flush()
   else:
      try:
         objname = 'foo' # FIX THIS
      except:
         objname = 'Weight file for %s' % infile
      outname = infile.replace('.fits','%s.fits' % outsuff)
      pf.PrimaryHDU(whtdat).writeto(outname)
      del whtdat

   """ Clean up """
   del indat

# ------------------------------------------------------------------------------

def make_wht_for_final(infiles, medfile, nsig, inwht_suff='.weight.fits',
                       outwht_suff='_wht.fits', flag_posonly=False, 
                       medwhtfile='default'):
   """
   Creates a weight file for each input image to a final swarp call.
   This weight file will assign zero weight to pixels that differ by
   more than nsig*rms from the median-stacked image that was created
   in an earlier step.
   Therefore, this function is a lot like the "blot" step in multidrizzle.

   Inputs:
      infiles      - list of individual input fits files that will be compared
                      to the median stack.  These will probably be called 
                      *resamp.fits
      medfile      - median-stacked fits file
      nsig         - minimum sigma difference between an individual input image
                     and the median stacked image that will cause a pixel to be
                     flagged.
      flag_posonly - Sets flagging behavior.  If flag_posonly=True then only
                     flag pixels where the value in the individual image is
                     more than nsig*rms _larger_ than the value in the
                     median-stacked image.
                     If False, then flag pixels for which the absolute
                     value of the difference is more than nsig*rms
                     Explanation: Set to True to flag, e.g., cosmic rays in
                     the individual image and to avoid flagging
      medwhtfile   - Name of the weight file associated with medfile.  The
                     default value ('default') will take the name of the
                     median-stacked image (medfile) and replace '.fits' with
                     '_wht.fits'
   """

   import imfuncs as imf
   from math import sqrt

   """ Make sure that the input is either a list or a single file """

   if type(infiles) is str:
      print ""
      print "Single input file"
      tmplist = [infiles,]
   elif type(infiles) is list:
      print ""
      print "Input file list with %d members" % (len(infiles))
      tmplist = infiles
   else:
      print ""
      print "Warning.  Input frames need to be either a list of files "
      print " (python type==list) or a single input file name."
      print ""
      return

   """ First loop through the input files to get information """
   gain = np.zeros(len(tmplist))
   bkgd = np.zeros(len(tmplist))
   fscal = np.zeros(len(tmplist))
   x1 = np.zeros(len(tmplist),dtype=int)
   y1 = np.zeros(len(tmplist),dtype=int)
   rmssky = np.zeros(len(tmplist))

   print ''
   print ' Input file                  gain  <bkgd>   fscal  rms_sky rms_scal'
   print '---------------------------- ---- -------- ------- ------- --------'
   for i in range(len(tmplist)):

      """ Load the header information """
      f = tmplist[i]
      try:
         hdr = pf.getheader(f)
      except:
         print ""
         print "ERROR: Could not open %s" % f
         print ""
         return

      origfile = f.replace('resamp.fits','fits')
      try:
         orighdr = pf.getheader(origfile)
      except:
         print ""
         print "ERROR: Could not open %s" % origfile
         print ""
         return

      """ Set the relevant values """
      x1[i]     = int(hdr['comin1'] - 1)
      y1[i]     = int(hdr['comin2'] - 1)
      bkgd[i]   = hdr['backmean']
      fscal[i]  = hdr['flxscale']
      gain[i]   = orighdr['gain']
      rmssky[i] = sqrt(bkgd[i] / gain[i])

      """ Print out the relevant information and clean up """
      print '%-28s %4.2f %8.2f %7.5f %7.3f %8.3f' \
          % (f[:-5],gain[i],bkgd[i],fscal[i],rmssky[i],rmssky[i]*fscal[i])
      del hdr,orighdr

   gainmean = gain.mean()
   print '-------------------------------------------------------------------'
   print ' %3d input files             gain  <bkgd>   fscal  rms_sky rms_scal'\
       % (len(tmplist))
   print '---------------------------- ---- -------- ------- ------- --------'
   print '    Mean values              %4.2f %8.2f %7.5f %7.3f %8.3f' % \
       (gainmean,bkgd.mean(),fscal.mean(),rmssky.mean(), \
           rmssky.mean()*fscal.mean())

   """ Set up the weight file associated with the median-stacked image """
   if medwhtfile == 'default':
      medwhtfile = medfile.replace('.fits','_wht.fits')

   """ Open up the median file and its associated weight image"""
   medall = pf.getdata(medfile)
   medwhtall = pf.getdata(medwhtfile)

   """ Loop through the input files """
   epsil = 1.e-20
   count = 1
   print ''
   print 'Updating weight files, using diff > %d sigma' % nsig
   print '--------------------------------------------------------------------'
   for i in range(len(tmplist)):
      """ Set up file names """
      f = tmplist[i]
      inwhtfile = f.replace('.fits',inwht_suff)
      outwhtfile = f.replace('.fits',outwht_suff)

      """ Load input resamp.fits file """
      print '%s -- File %d of %d' % (f[:-5],count,len(tmplist))
      print '-----------------------------------------'
      print 'Loading individual exposure: %s' % f
      try:
         indat,hdr = pf.getdata(f,header=True)
      except:
         print ""
         print "ERROR: Could not open %s" % f
         print ""
         return

      """ Load the associated weight file data """
      try:
         inwht,whthdr = pf.getdata(inwhtfile,header=True)
      except:
         print ""
         print "ERROR: Could not open %s" % inwhtfile
         print ""
         return


      """ Set up the relevant region to be examined """
      x2 = int(x1[i] + indat.shape[1])
      y2 = int(y1[i] + indat.shape[0])

      """ 
      Make the cutout of the median-stack image and then take the 
      difference between this file and the scaled input individual file
      """
        #meddat = (pf.getdata(medfile))[y1[i]:y2,x1[i]:x2]
        #medwht = (pf.getdata(medwhtfile))[y1[i]:y2,x1[i]:x2]
        #meddat = medhdu[0].section[y1[i]:y2,x1[i]:x2]
        #medwht = medwhthdu[0].section[y1[i]:y2,x1[i]:x2]
      print 'Selecting data from full stacked median file'
      meddat = medall[y1[i]:y2,x1[i]:x2]
      medwht = medwhtall[y1[i]:y2,x1[i]:x2]
      print 'Calculating difference image'
      diff = np.zeros(indat.shape)
      whtmask = inwht>0
      diff[whtmask] = indat[whtmask] * fscal[i] - meddat[whtmask]
      del whtmask

      """ 
      Get an estimate of the RMS noise in the data.  Since we are creating
      a difference, indat - meddat, the final rms will be
          sigma_diff^2 = sigma_ind^2 + sigma_med^2
      The input variances will come from:

          sigma_ind^2 = (inddat+bkgd) / gain
            - inddat is a background-subtracted image with units of ADU.
              Therefore, add back in the background and then divide by
              the gain to get the variance in each pixel.
              This is because sigma_ADU^2 = ADU/gain, assuming that
              N_e = gain * ADU is a Poisson process.

          sigma_med^2 = (1. / medwht) + [for SNR>1: meddat/gain] 
            - Swarp produces an inverse-variance weight map
              NOTE: However, this does NOT include Poisson noise from the
              objects in the image.  This is because SExtractor computes this
              Poisson noise when it does its object detections, and therefore 
              expects any input weight file to not include the Poisson noise.

      """
      print 'Masking and then calculating variance from stacked median'
      medvar = np.zeros(indat.shape)
      medmask = (inwht>0) & (np.absolute(medwht>epsil))
      medvar[medmask] = 1. / medwht[medmask]
      del medmask
      print 'Masking and then adding Poisson noise from stacked median'
      medmask = meddat > np.sqrt(medvar)
      medvar[medmask] += meddat[medmask] / gainmean
      del medmask
      #del medwht,meddat
      print 'Masking and then calculating variance in individual image'
      indvar = np.zeros(indat.shape)
      mask = (inwht>0) & ((indat + bkgd[i]) > 0.)
      indvar[mask] = fscal[i]**2 * (indat[mask]+bkgd[i])/gain[i]
      print 'Calculating combined rms'
      rms = np.sqrt(medvar + indvar)
      del medvar,indvar,mask,indat

      """ 
      Flag pixels that deviate by more than nsig sigma from the 
      median-stacked image
      """
      print 'Flagging pixels that differ by more than %d sigma from median'\
          % nsig
      if flag_posonly:
         blotmask = diff > nsig*rms
      else:
         blotmask = np.absolute(diff) > nsig*rms
      print '%-38s %9d' \
          % (f[:-5],blotmask.sum())

      """ Debugging step(s) """
      #foodat = np.ones(diff.shape)
      #foodat[blotmask] = 0
      #fooname = f.replace('.fits','_foo.fits')
      #pf.PrimaryHDU(foodat).writeto(fooname,clobber=True)
      #del foodat

      #rmsname = f.replace('.fits','_rms.fits')
      #pf.PrimaryHDU(rms,hdr).writeto(rmsname,clobber=True)

      #snr = np.zeros(diff.shape)
      #rmsmask = rms>0
      #snr[rmsmask] = diff[rmsmask] / rms[rmsmask]
      #snrname = f.replace('.fits','_snr.fits')
      #pf.PrimaryHDU(snr,hdr).writeto(snrname,clobber=True)
      #del snr

      """ Write out a new weight file with the newly-flagged pixels """
      inwht[blotmask] = 0
      print '%s -> %s' % (inwhtfile,outwhtfile)
      pf.PrimaryHDU(inwht,whthdr).writeto(outwhtfile,clobber=True)
      count += 1
      print ''

      """ Close the files for this loop """
      del diff,blotmask,rms,inwht

   """ Clean up """
   #medhdu.close()
   #medwhthdu.close()
   del medall,medwhtall

# ---------------------------------------------------------------------------

def add_exptime(inlist, exptime, exptkey='exptime', hext=0, verbose=True):
   """
   Given a list of fits files, adds to each one an EXPTIME header keyword
   with a value based on the passed exptime parameter.
   ** NOTE ** The exptime parameter can either be a numerical value, in
     which case it is interpreted as the value for EXPTIME, or it can
     be a string in which case it is interpreted as the name of a reference
     file that contains a valid exposure time keyword (designated by
     the passed exptkey parameter) that will be copied into the input files.

   Inputs:
     inlist  - list of fits files to which the keyword will be added
     exptime - can take one of two forms:
               (1) the value of the exposure time to add to each file
               (2) the name of a reference file that contains a valid
                   exposure time in its header.  This exposure time
                   (designated by the exptkey keyword in the header) will
                   be copied into the input list.
     exptkey - if the exptime parameter is the name of a reference file, then
               this designates the keyword in the reference file header that
               contains the exposure time value.  Default is 'exptime'
     hext    - HDU to modify (default = 0)
     verbose - Set to True (the default) for some status
   """

   """ Check format of the exptime parameter """
   if type(exptime) is str:
      try:
         hdr = pf.getheader(exptime)
      except:
         print ''
         print 'ERROR: Unable to open fits file %s' % exptime
         print ''
         return
      try:
         texp = hdr[exptkey]
      except:
         print 'ERROR: Unable to read %s keyword in %s fits file' % \
             (exptkey,exptime)
         del hdr
         print ''
         return
   elif type(exptime) is float:
      texp = exptime
   elif type(exptime) is int:
      texp = float(exptime)
   else:
      print ''
      print 'ERROR: exptime needs to be a number or name of a reference file'
      print ''
      return

   """ Put the desired exposure time into the input files """
   for i in inlist:
      hdu = pf.open(i,mode='update')
      hdr = hdu[hext].header
      hdr.update('exptime',texp)
      hdu.flush()
      if verbose:
         print 'Updated %s with EXPTIME=%.2f' % (i,texp)

# ---------------------------------------------------------------------------

def make_texp_map(infiles, texp, whtext='wht', outext='texp'):
   """
   Uses the weight files in the input list to create individual exposure
   time maps for each resampled image.
   NOTE: The input weight files need to be the resampled and then 
   modified ones that are created with the make_wht_for_final function

   Inputs:
      infiles  - list of input files, perhaps created with a glob.glob call
      texp     - exposure time for the original image
      whtext   - extension on input files, which should be the weight
                 files.  The default, 'wht' means that the files in infiles
                 are expected to be called [root]_wht.fits
      outext   - extension used for output exposure-time file name. In other
                   words, for an input file of [root]_[whtext].fits the output 
                   file name will be [root]_[outext].fits
                 Default: 'texp'

   Outputs:
      Each input file called [root]_[whtext].fits will produce an output 
       exposure time file called [root]_[outext].fits
   """

   """ Make sure that the input is either a list or a single file """

   if type(infiles) is str:
      print ""
      print "Single input file"
      tmplist = [infiles,]
   elif type(infiles) is list:
      print ""
      ncount = len(infiles)
      print "Input file list with %d members" % ncount
      tmplist = infiles
   else:
      print ""
      print "Warning.  Input frames need to be either a list of files "
      print " (python type==list) or a single input file name."
      print ""
      return

   """ Run through the list """

   print ''
   print 'Making individual exposure-time maps'
   print '------------------------------------'
   count = 1
   for f in tmplist:
      """ Get input file"""
      print 'File %d of %d' % (count,ncount)
      print 'Input file:  %s' % f
      data,hdr = pf.getdata(f,header=True)

      """ Create the exposure time data as integer """
      data[data>0] = texp
      data = data.astype(int)
      hdr['object'] = 'Exposure time map'

      """ Write the output file """
      outfile = f.replace('%s.fits' % whtext,'%s.fits' % outext)
      pf.PrimaryHDU(data,hdr).writeto(outfile)
      print 'Output file: %s' % outfile
      print ''

      """ Clean up """
      del data,hdr
      count += 1

# -----------------------------------------------------------------------

def split_imext(infits, next, outroot=None):
   """
   Splits a multi-extension fits file into its separate images.
   The assumption, for now, is that the file has the following structure:
      HDU0 - General information
      HDU1 - Image HDU 1
      ...
      HDUN - Image HDU N

   Inputs:
      infits  - input multiextension fits file
      next    - number of extensions to split out
      outroot - keep at the default (None) to have the output files have the
                same root as the input.  If this is not None, the output
                files will have a different root.
                Example: split_imext('myfile.fits',2) will produce output files
                    myfile_1.fits and myfile_2.fits
                Example: split_imext('myfile.fits',2,'foo') will produce
                    foo_1.fits and foo_2fits
   Output:
      the N individual fits files
   """

   """ Set up the output naming scheme """
   if outroot is None:
      outroot = infits[:-5]

   """ Open the input file and split it """
   hdu = pf.open(infits)
   print ''
   print 'Splitting input file: %s' % infits
   print '-------------------------------------------------------------'
   for i in range(1,next+1):
      data = hdu[i].data.copy()
      outfile = '%s_%d.fits' % (outroot,i)
      pf.PrimaryHDU(data,hdu[i].header).writeto('%s' % outfile)
      del data
      print ' Split %s' % outfile

   """ Close """
   hdu.close()

# -----------------------------------------------------------------------


def read_wcsinfo(fitsfile, inhdu=0, verbose=True):
   """
   Reads the wcs information from a fits file and returns that information.

   NOT YET FUNCTIONAL
   """

   validwcs = False

   try:
      hdulist = imf.open_fits(fitsfile)
   except:
      return
   hdr = hdulist[inhdu].header

   """ Get the CTYPE, CRPIX, CRVAL pairs first """
   ctype = []
   crpix = np.zeros(hdr['naxis'])
   crval = np.zeros(hdr['naxis'])
   for i in range(hdr['naxis']):
      typecard = 'ctype%d' % (i+1)
      pixcard = 'crpix%d' % (i+1)
      valcard = 'crval%d' % (i+1)
      ctype.append(hdr[typecard])
      crpix[i] = hdr[pixcard]
      crval[i] = hdr[valcard]

   """ Search for CD matrix """
   cdmatx = np.zeros((hdr['naxis'],hdr['naxis']))
   validwcs = True
   for i in range(hdr['naxis']):
      for j in range(hdr['naxis']):
         cdcard = 'cd%d_%d' % (i+1,j+1)
         try:
            cdmatx[i,j] = hdr[cdcard]
         except:
            validwcs = False

   """ Search for CDELT and CROT headers """
   if not validwcs:
      validwcs = True
      cdelt = np.zeros(hdr['naxis'])
      crota = np.zeros(hdr['naxis'])
      for i in range(hdr['naxis']):
         deltcard = 'cdelt%d' % (i+1)
         rotcard = 'crota%d' % (i+1)
         try:
            cdelt[i] = hdr[deltcard]
         except:
            validwcs = False
         try:
            crota[i] = hdr[rotcard]
         except:
            crota[i] = 0.0
      if validwcs:
         """ Convert to CD matrix (on hold for now) """

   if not validwcs:
      print ""
      print "ERROR: read_wcsinfo.  No valid pixel scale header cards found"
      print ""
      return

   """ Print out the results, if desired """
   if(verbose and validwcs):
      print ""
      print "WCS information for %s" % fitsfile
      print "----------------------------------------------------------"
      for i in range(hdr['naxis']):
         print "%8s  %9.3f %13.9f " % (ctype[i],crpix[i],crval[i])

   """ Make the structure for returning the info """
   

# -----------------------------------------------------------------------


def make_wcs_from_tel_pointing(infile, pixscale, rotatekey=None, 
                               rakey='ra', deckey='dec', hext=0):
   """
   Many fits files produced by cameras on ground-based telescopes
   do not come with full WCS header cards.  However, most do store the
   (rough) telescope pointing information in the RA and DEC header cards.
   This function uses those and the input pixel scale to create a fake
   output header.
   This function takes an input fits file and creates a temporary header
   that contains the rough WCS information in a more useful format.
   It writes the new WCS information into the input file, and also
   returns the new header in case it is needed for other files

   Inputs:
      infile    -  input fits file
      pixscale  -  desired pixel scale in arcsec/pix for the output header 
                   (should be the expected pixel scale for the camera)
      rotatekey -  [OPTIONAL] if the image was taken with a rotation, that can
                   be incorporated into the output header.  In that case,
                   rotatekey needs to be set to the keyword for the header
                   card that contains the rotation value (these card keywords
                   are very non-standardized).
      rakey     -  [OPTIONAL] header keyword for RA of pointing, if different
                   from the default value of 'RA'
      deckey    -  [OPTIONAL] header keyword for Dec of pointing, if different
                   from the default value of 'Dec'
      hext    -  [OPTIONAL] HDU number to use.  Default=0 (i.e., the primary
                   HDU)
   """

   print ""
   """ Open the input fits file """
   try:
      hdu = pf.open(infile,mode='update')
   except:
      print 'ERROR: Could not open %s' % infile
      print ''
      exit()
   print 'Opened input fits file: %s' % infile
   try:
      inhdr = hdu[hext].header
   except:
      print ''
      print 'ERROR: Could not open header for HDU %d in %s' % (infile,hext)
      print ''
      hdu.close()
      exit()
   print ''
   hdr_error = False

   """ Read in the RA and Dec pointing values """
   try:
      ra_tel = inhdr[rakey]
   except:
      print "RA designated by keyword %s was not found in header" % rakey
      print ""
      hdr_error = True
   try:
      dec_tel = inhdr[deckey]
   except:
      print "Dec designated by keyword %s was not found in header" % deckey
      print ""
      hdr_error = True

   """ Get size of image """
   try:
      xsize = inhdr['naxis1']
   except:
      print "Problem reading x-axis size"
      print ""
      hdr_error = True
   try:
      ysize = inhdr['naxis2']
   except:
      print "Problem reading y-axis size"
      print ""
      hdr_error = True

   """ Exit if any problems before this point """
   if hdr_error:
      print ''
      print 'ERROR: Not continuing'
      print ''
      hdu.close()
      exit()

   """ Create the temporary output header """
   try:
      outhdr = wcsmwa.make_header(ra_tel,dec_tel,xsize,ysize,pixscale)
   except:
      print "ERROR. Could not create output header"
      print ""
      hdu.close()
      exit()

   """ Rotate if requested """
   if rotatekey is not None:
      try:
         rot = inhdr[rotatekey]
      except:
         print "ERROR. Could not read rotation from header with keyword %s" %\
             rotatekey
         print ""
         return
      print "Rotating header by %7.2f degrees" % rot
      rothdr = wcsmwa.rotate_header(outhdr,rot)
      outhdr = rothdr.copy()

   print "Created header from telescope pointing info."
   print outhdr
   print 'Updating fits file header'
   wcskeys = ['ctype1','ctype2','crval1','crval2','crpix1','crpix2',
              'cd1_1','cd2_2']
   for k in wcskeys:
      inhdr[k] = outhdr[k]
               
   hdu.flush()
   return outhdr

# -----------------------------------------------------------------------


def make_wcs_from_ref_tel(reffile, infile, pixscale, rotatekey=None,
                          hext=0, rakey='ra', deckey='dec'):
   """
   Creates a full WCS header from the telescope pointing information in
   the reffile, and then copies the WCS information into the infiles.

   Inputs:
      reffile   -  reference file containing the telescope pointing info in
                   its RA and DEC fits headers.
      infile    -  input file to which WCS information will be copied
      pixscale  -  desired pixel scale for the output header (should be the
                   expected pixel scale for the camera)
      rotatekey -  [OPTIONAL] if the reffile image was taken with a rotation, 
                   that can be incorporated into the output header.  In that 
                   case, rotatekey needs to be set to the keyword for the header
                   card that contains the rotation value (these card keywords
                   are very non-standardized).
      hext    -  HDU number to use.  Default=0 (i.e., the primary HDU)
      rakey     -  [OPTIONAL] header keyword for RA of pointing, if different
                   from the default value of 'RA'
      deckey    -  [OPTIONAL] header keyword for Dec of pointing, if different
                   from the default value of 'Dec'
   """

   print ""

   """ Open the files """
   try:
      refhdu = imf.open_fits(reffile)
   except:
      print "ERROR.  Could not read reference file %s" % reffile
      return
   refhdr = refhdu[hext].header

   try:
      inhdu = imf.open_fits(infile,mode='update')
   except:
      return
   inhdr = inhdu[hext].header

   """ Make the temporary header with full WCS info from reffile """
   try:
      refwcs = make_wcs_from_tel_pointing(refhdr,pixscale,rotatekey,
                                          rakey,deckey)
   except:
      print "Could not create WCS header from %s" % reffile
   refinfo = wcsmwa.parse_header(refwcs)

   """ Copy the information into the input file """
   inhdr.update('ctype1',refwcs['ctype1'])
   inhdr.update('crpix1',refwcs['crpix1'])
   inhdr.update('crval1',refwcs['crval1'])
   inhdr.update('ctype2',refwcs['ctype2'])
   inhdr.update('crpix2',refwcs['crpix2'])
   inhdr.update('crval2',refwcs['crval2'])
   inhdr.update('cd1_1',refinfo[2][0,0])
   inhdr.update('cd2_2',refinfo[2][1,1])
   inhdr.update('cd1_2',refinfo[2][0,1])
   inhdr.update('cd2_1',refinfo[2][1,0])
               
   inhdu.flush()

# -----------------------------------------------------------------------

def hdr_offsets(files, pixscale=0, rakey=None, deckey=None, rot=None,
                oformat='pix', hext=0, verbose=True):
   """
   Uses the telescope-provided RA and Dec in the fits header to calculate
   an initial guess for the offsets between the input files.

   The default behavior for each file is to:
     1. Try to get the RA and Dec values from the CRPIXn header keywords
        and the pixel scale from the CDn_m keywords
     2. If step 1 fails for the RA and Dec, then try to get the RA and
        Dec information from the telescope pointing keywords defined by
        the rakey and deckey parameters.

   Inputs:
      files:     list of files for which to calculate offsets
      pixscale:  pixel scale of the input files, in arcsec/pix.  The default
                 value of 0 means that the pixel scale should be derived
                 from the CD matrix (e.g., CD1_1, etc.) in the fits header
      rakey:     name of the fits header keyword for RA. Default=None implies
                 that the default behavior described above is followed.
                 Set to a non-None value to override default
      deckey:    name of the fits header keyword for Dec.  Default=None implies
                 that the default behavior described above is followed.
                 Set to a non-None value to override default  
      rot:       Rotation of coordinate axes, North through East, in degrees.
                 Set this parameter to override default setting, which is the
                 rotation set by the CD matrix (if it exists in the fits header)
                 or 0.0 (if the CD matrix doesn't exist)
      oformat:   Output format for offsets.  Default is 'pix' for pixel offsets.
                 The only other option is 'arcsec' for offsets in arcseconds.
      hext:      HDU extension.  Default is hext=0, i.e., the primary HDU,
                 which is also the only HDU for many files
      verbose:   set to True for verbose output.  Default=True
   """

   # Initialize containers
   ra   = np.zeros(len(files))
   dec  = np.zeros(len(files))
   dx   = np.zeros(len(files))
   dy   = np.zeros(len(files))
   cdmatx = []

   # Print out some header information
   if verbose:
      print ""
      print \
          "   File                RA       Dec     xpxscl ypxscl rotation"
      print \
          "------------------- --------- --------- ------ ------ --------"
   

   # Read in the WCS information and process it
   count = 0
   for f in files:

      # Open the file and get the appropriate header info
      hdulist = imf.open_fits(f)
      if hdulist is None:
         return None,None
      hdr = hdulist[hext].header

      # Initialize the containers
      ra[count] = np.nan
      dec[count] = np.nan

      # Start with the default behavior of trying to read in the WCS info
      # If this is successful, set the parameter values based on the WCS info
      foundwcs = True
      badwcs = False
      try:
         wcsinfo = wcsmwa.parse_header(hdr)
      except:
         foundwcs = False
      if foundwcs:
         cd = wcsinfo[2]
         if wcsinfo[4] == 'RA':
            ra[count] = wcsinfo[0][1]
         else:
            ra[count] = wcsinfo[0][0]
         if wcsinfo[4] == 'DEC':
            dec[count] = wcsinfo[0][1]
         else:
            dec[count] = wcsinfo[0][0]
      else:
         cd = np.array([[-1.0,0.0],[0.0,1.0]])

      # Calculate the RA and Dec based on whether the override parameters are set
      if rakey is not None:
         try:
            ra[count] = hdr[rakey]
         except:
            print ""
            print "ERROR: Could not read RA using keyword %s from file %s" % \
                (rakey,f)
            return None,None

      if deckey is not None:
         try:
            dec[count] = hdr[deckey]
         except:
            print ""
            print "ERROR: Could not read DEC using keyword %s from file %s" % \
                (deckey,f)
            return None,None

      # Set the CD matrix if pixscale is set to override. Default rotation = 0
      if pixscale > 0.:
         if rot is None:
            rot = 0.0
         cd = coords.rscale_to_cdmatrix(pixscale,rot,verbose=False)

      # Print out information about the input files
      if verbose:
         cd1,cd2,cr2 = coords.cdmatrix_to_rscale(cd)
         print "%-18s  %9.5f %+9.5f %+6.3f %+6.3f %+6.1f" \
             % (files[count],ra[count],dec[count],3600.*cd1,3600.*cd2,
                cr2*180./pi)

      # Finally, save the constructed CD matrix
      cdmatx.append(cd)

      count += 1
      del hdr
      del hdulist

   # Calculate the offsets
   dalpha = 3600. * ((ra - ra[0]) * cos(pi * dec[0]/180.))
   ddelta = 3600. * (dec - dec[0])
   if verbose:
      print ""
      print \
          "   File              da(asec) dd(asec)  dx(pix) dy(pix)"
      print \
          "-------------------  -------- --------  ------- -------"
   for i in range(dx.size):
      dx[i],dy[i] = coords.darcsec_to_dpix(dalpha[i],ddelta[i],cdmatx[i])
      if verbose:
         print "%-18s    %+7.3f  %+7.3f  %+7.2f %+7.2f" \
             % (files[i],dalpha[i],ddelta[i],dx[i],dy[i])

   # Return the offsets
   if oformat == 'arcsec':
      return dalpha,ddelta
   else:
      return dx,dy

# -----------------------------------------------------------------------

def plot_hdr_offsets(files, pixscale=0, rakey=None, deckey=None, rot=None, 
                     oformat='pix', hext=0, verbose=True):
   """
   Makes an initial calculation for the offsets between the input files based on
   information in the fits header cards (via hdr_offsets).  
   These offsets get plotted in the desired units, either pixels (the default)
   or arcseconds.

   Inputs:
      files:     list of files for which to calculate offsets
      pixscale:  pixel scale of the input files, in arcsec/pix.  The default
                 value of 0 means that the pixel scale should be derived
                 from the CD matrix (e.g., CD1_1, etc.) in the fits header
      rakey:     name of the fits header keyword for RA. Default=None implies
                 that the default behavior described in hdr_offsets is followed.
                 Set to a non-None value to override default
      deckey:    name of the fits header keyword for Dec.  Default=None implies
                 that the default behavior described in hdr_offsets is followed.
                 Set to a non-None value to override default           
      oformat:   Output format for offsets.  Default is 'pix' for pixel offsets.
                 The only other option is 'arcsec' for offsets in arcseconds.
      hext:      HDU extension.  Default is hext=0, i.e., the primary HDU,
                 which is also the only HDU for many files
      verbose:   set to True for verbose output.  Default=True
   """

   # Set up
   from matplotlib import pyplot as plt

   # Get the offsets
   dx,dy = hdr_offsets(files,pixscale,rakey,deckey,rot,oformat,hext,verbose)

   if dx is None or dy is None:
      print "ERROR. Failed to calculate offsets"
      return

   # Plot the offsets
   plt.scatter(dx,dy,marker='+')
   if oformat=='arcsec':
      plt.xlabel(r'd$\alpha$ (arcsec)')
      plt.ylabel(r'd$\delta$ (arcsec)')
   else:
      plt.xlabel('dx (pixels)')
      plt.ylabel('dy (pixels)')

# -----------------------------------------------------------------------------

def xcorr_offsets(infiles, hext=0):
   """
   Estimates the pixel shifts between dithered exposures on a field by doing
   a cross-correlation between the images.  The cross-correlation is done
   via fft.

   Inputs:
      infiles  - input files
      hext     - HDU of the data (default = 0)

   Output:
      xshift   - returned array of x offsets
      yshift   - returned array of y offsets
   """

   from scipy import fftpack

   """ Set the first file as the reference """
   ref = imf.Image(infiles[0])
   refdat = ref.hdu[hext].data.astype(np.float32)
   if refdat.shape[1]%2 == 0:
      xcent = refdat.shape[1] / 2
   else:
      xcent = (refdat.shape[1] - 1)/2
   if refdat.shape[0]%2 == 0:
      ycent = refdat.shape[0] / 2
   else:
      ycent = (refdat.shape[0] - 1)/2
   Fref = fftpack.fft2(refdat)

   """ Loop through the list of files """
   xshift = np.zeros(len(infiles))
   yshift = np.zeros(len(infiles))
   for i in range(1,len(infiles)):
      comp = imf.Image(infiles[i])
      compdat = comp.hdu[hext].data.astype(np.float32)
      Fprod = Fref * fftpack.fft2(compdat).conj()
      xcorr = fftpack.ifft2(Fprod)
      corrshift = fftpack.fftshift(xcorr)
      ymax,xmax = np.unravel_index(corrshift.argmax(),corrshift.shape)
      xshift[i] = xmax - xcent
      yshift[i] = ymax - ycent
      del compdat,Fprod,xcorr,corrshift
      #comp.close()
      del comp

   """ Clean up and return shifts """
   del refdat,Fref
   #ref.close()
   del ref
   return xshift,yshift
   
# ------------------------------------------------------------------------------

def coadd_intshift(infiles, xshifts, yshifts, outfile, origsize=False, hext=0):

   """ For now assume that all of the inputs files have the same size """
   ref = imf.Image(infiles[0])
   ysize,xsize = ref.hdu[hext].data.shape
   tmpxshift = xshifts - xshifts.min()
   tmpyshift = yshifts - yshifts.min()
   if origsize:
      xmed = np.median(tmpxshift)
      ymed = np.median(tmpyshift)
   del ref

   """ Modify the size based on the shifts """
   dimx = xsize + int(np.ceil(xshifts.max() - xshifts.min()))
   dimy = ysize + int(np.ceil(yshifts.max() - yshifts.min()))

   """ Create containers for the output data and nexp info """
   odat = np.zeros((len(infiles),dimy,dimx))
   nexp = np.zeros((len(infiles),dimy,dimx))

   """ Load data into the containers """
   for i in range(len(infiles)):
      tmpim = imf.Image(infiles[i])
      tmpdat = tmpim.hdu[hext].data.copy()
      xstart = int(tmpxshift[i])
      xend = int(xstart + xsize)
      ystart = int(tmpyshift[i])
      yend = int(ystart + ysize)
      odat[i,ystart:yend,xstart:xend] = tmpdat.copy()
      nexp[i,ystart:yend,xstart:xend] = 1.
      del tmpdat,tmpim

   """ Take the weighted average """
   nexp[nexp==0] = 0.1
   outdat = odat.sum(axis=0) / nexp.sum(axis=0)
   if origsize:
      outdat = outdat[int(ymed):int(ymed+ysize),int(xmed):int(xmed+xsize)]
   pf.PrimaryHDU(outdat).writeto(outfile,clobber=True)
   del odat,nexp

# ------------------------------------------------------------------------------

def fixpix_wht(datafile, whtfile, outfile=None, boxsize=11, datahdu=0,
               whthdu=0):
   """
   Replaces bad pixels in a fits image (datafile) based on the associated
   weight file.  Pixels for which the weight file is zero will be replaced
   by the median-filter value in the data file (i.e., a version of the
   data file is created where each pixel is replaced by the median pixel
   value in a boxsize x boxsize box centered on the pixel).
   If the optional outfile parameter is given, then the results are written
   out to a new file.  If not, then the data file is modified in place.
   """

   from scipy.ndimage import filters

   """ Open the input files """
   if outfile is None:
      dhdulist = imf.open_fits(datafile,'update')
      data = dhdulist[datahdu].data
   else:
      dhdulist = imf.open_fits(datafile)
      data = dhdulist[datahdu].data.copy()

   whdulist = imf.open_fits(whtfile)
   whtdat = whdulist[whthdu].data.copy()

   """ Median filter the data """
   meddat = filters.median_filter(data,boxsize)

   """
   Replace the bad data (where whtdat==0) with the median-filtered values
   """
   mask = whtdat==0
   data[mask] = meddat[mask]
   dhdulist[datahdu].header.update('fixpix','ccdredux.fixpix_wht',
                                   'Ran code to fix bad pixels')

   """ Write output """
   if outfile is None:
      dhdulist.flush()
   else:
      pf.PrimaryHDU(data,dhdulist[datahdu].header).writeto(outfile)

# -----------------------------------------------------------------------------

def fixpix_rms(datafile, rms_high=5., rms_low=5., rms_sigclip=3., outfile=None, 
               boxsize=11, datahdu=0, verbose=True):
   """
   Replaces bad pixels in a fits image (datafile), where the bad pixels are
   determined by:
    (1) running sigma_clip on the image data to get a clipped mean and sigma
    (2) finding pixels where the counts are (> mean + rms_high*sigma) or
        (< mean - rms_low*sigma)
   Bad pixels will be replaced by the median-filter value in the data file 
   (i.e., a version of the data file is created where each pixel is replaced by 
   the median pixel value in a boxsize x boxsize box centered on the pixel).
   If the optional outfile parameter is given, then the results are written
   out to a new file.  If not, then the data file is modified in place.
   """

   from scipy.ndimage import filters

   """ Open the input file """
   if outfile is None:
      dhdulist = imf.open_fits(datafile,'update')
      data = dhdulist[datahdu].data
   else:
      dhdulist = imf.open_fits(datafile)
      data = dhdulist[datahdu].data.copy()

   """ Sigma clip the input data and set the bad pixel mask """
   m, s = sigma_clip(data, verbose=verbose)
   mask = (data > m + rms_high * s) | (data < m-rms_low * s)

   """ Pass 1, replace bad pixels with the overall median value """
   data[mask] = np.median(data)

   """ Median filter the data """
   meddat = filters.median_filter(data, boxsize)

   """ Pass 2, replace the bad data with the median-filtered values """
   data[mask] = meddat[mask]

   """ Write output """
   if outfile is None:
      dhdulist.flush()
   else:
      pf.PrimaryHDU(data, dhdulist[datahdu].header).writeto(outfile)

# ---------------------------------------------------------------------------


def check_type(x, validtypes):
   """
   Checks the type of x versus a list of valid types.  If type(x) matches
   any member of the list, then return True, otherwise return False
   """
   type_ok = False
   tmp = type(x)
   for i in validtypes:
      if tmp == i:
         type_ok = True
   return type_ok

# ---------------------------------------------------------------------------


def make_var_with_poisson(infile, units, maskfile=None, gain=None, rdnoise=0.,
                          texp=None, origbkgd=0., epsilon=None, hext=0,
                          statcent=None, statsize=None, outfile=None,
                          outtype='var', outsnr=None, returnvar=True,
                          verbose=False):
   """
   Makes a variance file from an input image.  The variance will be based
   both on the variance in the pixel values in regions of blank sky and on
   the Poisson noise where there are objects in the image.

   Required inputs:
     infile - input image
     units  - units for the pixel values in the input image.  Allowed choices
              are:
                1. 'counts' - for input values in ADU
                2. 'cps'    - for input values in ADU/sec
                3. 'e-'     - for input values in electrons
                4. 'eps'    - for input values in electrons/sec

   Optional inputs:
     maskfile - file indicating bad pixels, where bad pixels have the value 0
                (or are larger than the epsilon parameter, if that parameter
                is set to something besides None), and good pixels have 
                non-zero values (or values > epsilon, if epsilon is not None).
                The mask file can be a weight file or bad pixel mask produced
                by standard image processing.
     gain     - gain (in e-/ADU).  Only used if units is 'counts' or 'cps'
     texp     - exposure time of the input image.  Only used if units is
                'cps' or 'eps'
   """

   """
   Start with some error checking based on the value chosen for the units
   """
   units_ok = False
   validtypes = [float, np.float32, np.float64, np.ndarray, int]

   if (units == 'counts' or units == 'cps'):
      units_ok = True
      if gain is None:
         print 'Error: Selected are counts or cps but the gain value is not set'
         raise TypeError
      if check_type(gain, validtypes) is not True:
         print 'Error: gain must either be a single number or a 2-D array'
         raise TypeError

   if (units == 'cps' or units == 'eps'):
      units_ok = True
      if texp is None:
         print 'Error: Selected units are cps or eps but texp is not set'
         raise TypeError
      if check_type(texp, validtypes) is not True:
         print 'Error: texp must either be a single number or a 2-D array'
         raise TypeError

   if units == 'e-':
      units_ok = True

   if units_ok is not True:
      print 'Error: Invalid value for units parameter'
      raise ValueError
   if verbose:
      print 'Finished unit and type check'

   """
   Load the input image and the mask image
   """
   im = imf.Image(infile, verbose=verbose)
   hdr = im.hdu[0].header.copy()
   data = im.hdu[hext].data
   if maskfile is not None:
      maskim = imf.Image(maskfile, verbose=verbose)
      mask = maskim.hdu[hext].data.copy()
      if epsilon is not None:
         mask[mask < epsilon] = 0
      mask[mask > 0] = 1
      mask = mask.astype(bool)
      del maskim
   else:
      mask = None

   """ Convert the texp to a 2-d array if it is not one already """
   if texp is not None:
      if type(texp) != np.ndarray:
         exptime = np.ones(data.shape) * texp
      else:
         exptime = texp

   """
   Determine the clipped RMS in the input image, which is a good proxy
   for the RMS in the blank sky pixels.  Use the mask to exclude bad
   pixels before doing the calculation.  Use the statcent and statsec
   parameters to define the region to use for the calculation, otherwise
   use the whole image.
   """
   im.sigma_clip(mask=mask, hext=hext, verbose=verbose)  # Put in the statsec later...
   if units == 'counts':
      rms = sqrt(im.rms_clip**2 - rdnoise**2)  # Think about the 'cps' case
   else:
      rms = im.rms_clip
   var = np.ones(data.shape) * rms**2

   """
   Now that the base level of the variance has been set, add the Poisson
   noise for pixels associated with objects.  Define these pixels as the
   ones where (after smoothing) the pixel value is more than 1 sigma above
   the clipped mean (which is a good proxy for the sky level).
   """
   bgsub = data - im.mean_clip
   snrsmo = filters.gaussian_filter(bgsub / rms, 1.)
   snrmask = (snrsmo > 1.) & (mask)

   """
   Determine the Poisson noise based on the units.  For the
   four cases, the additional amount to add for the relevant pixels is:
     1. units == 'counts':  gain * pixval
     2. [Needs to be double-checked]
     3. units == 'e-':      pixval
     4. units == 'eps':     pixval / texp
   """
   if verbose:
      print ''
      print 'Adding the Poisson noise'
      print '------------------------'

   if units == 'counts':
      var[snrmask] += gain * (data[snrmask] + origbkgd)
   elif units == 'cps':
      var[snrmask] += gain * (data[snrmask] + origbkgd) / exptime[snrmask]
   elif units == 'e-':
      var[snrmask] += (data[snrmask] + origbkgd)
   elif units == 'eps':
      var[snrmask] += (data[snrmask] + origbkgd) / exptime[snrmask]
   if verbose:
      print 'Done'

   """ Clean up """
   del im, mask, texp, gain

   """ Write to output files if requested """
   if outsnr is not None:
      pf.PrimaryHDU(bgsub / var**0.5).writeto(outsnr, clobber=True)
      del bgsub
   if outfile is not None:
      if outtype == 'rms':
         pf.PrimaryHDU(var**0.5).writeto(outfile, clobber=True)
      else:
         pf.PrimaryHDU(var).writeto(outfile, clobber=True)
   if returnvar:
      return var
   else:
      del var

# ---------------------------------------------------------------------------

def coadd_clean(infiles, medfile, outfile, whtsuff='wht', medwhtsuff=None):
   """
   Coadds the input data files after first masking out discrepant pixels
   """

   print ''
   print 'Input files'
   print '-------------------'
   for i in infiles:
      print i


   """ Check that images are same size """
   hdr0 = pf.getheader(infiles[0])
   x0 = hdr0["naxis1"]
   y0 = hdr0["naxis2"]
   print ''
   print 'Check that images are all the same size: %dx%d' % (x0, y0)
   print '-------------------------------------------------------'

   for i in infiles:
      hdr = pf.getheader(i)
      x = hdr["naxis1"]
      y = hdr["naxis2"]
      if x != x0 or y != y0:
         print ""
         print "Error: shape of %s does not match shape of %s" \
             % (i, infiles[0])
         print ""
         exit()
      else:
         print "%s: Shape matches" % i


   """ Get the median file data """
   med = imf.Image(medfile, verbose=False)
   meddata = med.hdu[0].data
   med.sigma_clip()
   noisemed = med.rms_clip

   """ initialize 3d arrays with zeros """
   insci  = np.zeros((len(infiles), y0, x0))
   wht    = np.zeros((len(infiles), y0, x0))
   newwht = np.zeros((len(infiles), y0, x0))

   """ Find and mask the bad pixels, and the area around them """
   print ''
   print 'Masking the bad pixels...'
   print '-------------------------'
   for i in range(len(infiles)):
      print '  %s' % infiles[i]
      '''Put data in arrays'''
      infile = infiles[i]
      whtfile = infile.replace('.fits', '_%s.fits', whtsuff)
      imi = imf.Image(infile, verbose=False)
      insci[i, :, :] = imi.hdu[0].data.copy()
      insci_i = insci[i, :, :].copy()
      whti = imf.Image(whtfile, verbose=False)
      wht[i, :, :] = whti.hdu[0].data.copy()
      tmpwht = np.ones(whti.hdu[0].data.shape)

      '''Mask pixels associated with weight of zero'''
      mask = whti.hdu[0].data != 0
      newwht[i, :, :][mask] = 1
      newwhti = newwht[i, :, :].copy()

      """
      Find bad pixels
      NOTE: right now this does not account for Poisson noise where there are
      real objects.  This needs to be fixed!
      """
      imi.sigma_clip(mask=mask)
      noisei = imi.rms_clip
      noise = sqrt(noisei**2 + noisemed**2)
      diffi = (insci_i - meddata) * newwhti
      mask2 = diffi > (3 * noise)

      ''' Grow the masked region around the bad pixels '''
      tmpwht[mask2] = 0
      tmp2 = filters.minimum_filter(tmpwht, size=3)
      # outfile = "tmp%d.fits" %i
      wht[i, :, :] *= tmp2
      # newwht[i,:,:] *= tmp2

      """ Do the weighted sum and also create the sum of the weights """
      print ''
      print 'Creating the weighted average image'
      datasum = (insci * wht).sum(axis=0)  # 2d array
      whtsum = wht.sum(axis=0)  # 2d array

      """ Avoid devision by zero errors """
      mask = whtsum == 0
      tmpwht = whtsum.copy()
      tmpwht[mask] = 1
      datasum[mask] = 0

      """' Finally, define the weighted average """
      whtedav = datasum / tmpwht
      print '  ...Done'

      """ Save the outputs """
      print 'Writing output to %s' % outfile
      pf.PrimaryHDU(whtedav).writeto(outfile, clobber=True)
      outwht = outfile.replace('.fits', '_%s.fits', whtsuff)
      pf.PrimaryHDU(whtsum).writeto(outwht, clobber=True)
      print ''
