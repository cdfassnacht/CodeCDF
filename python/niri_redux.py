"""
A series of functions to interact with the iraf/pyraf routines to reduce
imaging data from NIRI.  These are in the gemini package under iraf/pyraf
"""

import os
import numpy as n
import pyfits as pf
import astromatic as astrom
import astrom_simple as astsimp
from ccdredux import median_combine,sigma_clip,fixpix_wht
from scipy.ndimage import filters

#-----------------------------------------------------------------------

def make_flat(on_frames, off_frames, rawroot, band, rawdir='.'):
   """
   Makes an imaging flat-field file using the niflat task in pyraf
   """

   from pyraf import iraf
   iraf.gemini()
   iraf.niri()
   iraf.unlearn('nprepare')
   iraf.unlearn('niflat')

   """
   Set up the file lists that will be used in the various tasks
   """

   rawname = 'flats_%s_all_raw.list' % band
   npname  = 'flats_%s_all_nprepare.list' % band
   onname  = 'flats_%s_on.list' % band
   offname = 'flats_%s_off.list' % band
   outflat = 'Flat_%s.fits' % band

   """
   Create the lists
   """

   print ""
   print "Generating list of lamps-on files"
   print "---------------------------------"
   all_list = []
   f_in  = open(rawname,'w')
   f_np  = open(npname,'w')
   f_on  = open(onname,'w')
   f_off = open(offname,'w')
   for i in on_frames:
      infile = '%s%s.fits' % (rawroot,i)
      ofile  = 'np%s.fits' % i
      print " %s" % infile
      f_in.write('%s\n' % infile)
      f_np.write('%s\n' % ofile)
      f_on.write('%s\n' % ofile)
      all_list.append(ofile)
   print ""
   print "Generating list of lamps-off files"
   print "----------------------------------"
   for i in off_frames:
      infile = '%s%s.fits' % (rawroot,i)
      ofile  = 'np%s.fits' % i
      print " %s" % infile
      f_in.write('%s\n' % infile)
      f_np.write('%s\n' % ofile)
      f_off.write('%s\n' % ofile)
      all_list.append(ofile)
   f_in.close()
   f_np.close()
   f_on.close()
   f_off.close()

   """ 
   Run pyraf task nprepare to convert format of raw files to the one that is
   expected for subsequent tasks.
   """
   print ""
   print "Running nprepare on raw flat frames"
   print ""
   iraf.nprepare('@%s' %rawname,rawpath=rawdir,outimages='@%s'%npname)

   """
   Run pyraf task niflat to actually make the flat-field file
   """
   iraf.niflat('@%s'%onname,flatfile=outflat,lampsoff='@%s'%offname)


   """ Clean up """
   for i in all_list:
      os.remove(i)

   """ Give some information about the output file """
   print ""
   print "New flatfield file summary"
   print "------------------------------------------------------------------"
   flathdu = pf.open(outflat)
   flathdu.info()
   flathdu.close()

#-----------------------------------------------------------------------

def calib_1(sci_frames, rawroot, outroot, bpmfile, rawdir='.', obsdate=None, 
            do_nresid=False):
   """
   Does the non-linearity correction and then converts the linearized
    input files to the format expected by the Gemini pyraf tasks.  
   These steps are done through one Gemini and three pyraf tasks:
      nirlin    - to correct the raw files for nonlinearity
      nprepare  - to convert the linearized files to the appropriate format
      nresidual - [NOT IMPLEMENTED YET] to deal with persistence, if desired
      nisky     - to create the first-pass sky image
   """

   from pyraf import iraf
   iraf.gemini()
   iraf.niri()
   iraf.unlearn('nprepare')
   iraf.unlearn('nisky')
   import nirlin

   """
   Set up the file lists that will be used in the various tasks
   """

   rawname    = '%s_raw.list' % outroot
   lincorname = '%s_lincor.list' % outroot
   npname     = '%s_nprep.list' % outroot
   outname    = '%s_sky.fits' % outroot

   """
   Create the lists
   """

   print ""
   print "Generating list of science files"
   print "---------------------------------"
   raw_list = []
   lc_list = []
   #np_list = []
   f_in  = open(rawname,'w')
   f_lc  = open(lincorname,'w')
   f_np  = open(npname,'w')
   for i in sci_frames:
      infile = '%s%s.fits' % (rawroot,i)
      raw_list.append('%s/%s' % (rawdir,infile))
      if obsdate is not None:
         outframe = '%s_%s' % (obsdate,i)
      else:
         outframe = i
      lcfile = 'lincor%s.fits' % outframe
      ofile  = 'np%s.fits' % outframe
      print " %s" % infile
      f_in.write('%s\n' % infile)
      f_lc.write('%s\n' % lcfile)
      f_np.write('%s\n' % ofile)
      lc_list.append(lcfile)
      #np_list.append(ofile)
   f_in.close()
   f_lc.close()
   f_np.close()

   """
   Run nirlin to do the non-linearity correction
   """
   print ''
   print 'Running nirlin to correct for non-linearities'
   print '---------------------------------------------'
   print ''
   for i in range(len(raw_list)):
      nirlin.nirlin(raw_list[i],outputfile=lc_list[i])
   #   if i==0:
   #      npstring = np_list[0]
   #   else:
   #      npstring += ',%s' % np_list[i]

   """ 
   Run pyraf task nprepare to convert format of raw files to the one that is
   expected for subsequent tasks.
   """
   print ""
   print "Running nprepare on linearized raw frames"
   print ""
   print "%s %s %s %s" %(rawdir,lincorname,npname,bpmfile)
   #print "%s" % npstring
   tmplc = '@%s' % lincorname
   tmpnp = '@%s' % npname
   iraf.nprepare(tmplc,outimages=tmpnp,bpm=bpmfile)

   """
   Run pyraf task nisky to make the initial sky frame
   """
   iraf.nisky(tmpnp,outimage=outname)

#-----------------------------------------------------------------------

def reduce_sci(sci_frames, outroot, skyfile, flatfile, dodark=False, 
               darkfile='Dark.fits'):

   """
   Runs pyraf task nireduce to reduce the science frames using the 
   skyflat frame generated from these same science frames (using the
   calib_1 function) and the domeflat generated using the make_flat function.
   """

   from pyraf import iraf
   iraf.gemini()
   iraf.niri()
   iraf.unlearn('nprepare')
   iraf.unlearn('nireduce')
   """
   Set up the file lists that will be used in the various tasks
   """

   npname  = '%s_nprep.list' % outroot
   outname = '%s_ff.fits' % outroot

   """
   Create the lists
   """

   print ""
   print "Generating list of files"
   print "---------------------------------"
   f_out  = open(outname,'w')
   for i in sci_frames:
      ofile  = 'ff%s.fits' % i
      f_out.write('%s\n' % ofile)
   f_out.close()

   if dodark:
      iraf.nireduce('@%s'%npname,outimages='@%s'%outname,fl_sky=True,
                    skyimage=skyfile,fl_flat=True,flatimage=flatfile,
                    fl_dark=True,darkimage=darkfile)
   else:
      iraf.nireduce('@%s'%npname,outimages='@%s'%outname,fl_sky=True,
                    skyimage=skyfile,fl_flat=True,flatimage=flatfile,
                    fl_dark=False)

#-----------------------------------------------------------------------

def niri_bpm_from_sky(sky_frames, outbpm, inroot='ff', sigbpm=5.,
                      sigclip=3., outsky=None):
   """
   Given a list of flat-fielded frames (e.g., produced by the reduce_sci
   function), create a median sky from the science HDUs of the input files.
   Then create a bad pixel mask by selecting all pixels that are
   offset (in an absolute value sense) by more than sigbpm*rms from the 
   clipped mean, where rms is the standard deviation of the sigma-clipped 
   version of the median sky data.

   Inputs:
      sky_frames: list/array of input frame numbers
      outname:    name of output bad pixel mask file
      inroot:     prefix for input frames (default='ff'), so input files
                  will be named 'inroot''sky_frames[i]'.fits, e.g., ff123.fits
      sigbpm:     minimum number of sigma offset from the clipped mean used to 
                  define a bad pixel, where sigma is the rms of the clipped
                  data set (default=20.).  
                  So:  abs(data-clipped_mean) > sigbmp*clipped_rms
      sigclip:    clipping level for doing the sigma clip on the data set
   """

   """ Read in data (from HDU 1 in the files) and create temporary files """
   tmproot = '__tmp__'
   tmplist = []
   for i in sky_frames:
      infile = '%s%s.fits' % (inroot,i)
      data = pf.getdata(infile,1)
      tmpfile = infile.replace(inroot,tmproot)
      pf.PrimaryHDU(data).writeto(tmpfile,clobber=True)
      tmplist.append(tmpfile)
      del data

   """ 
   Take the median of the input data files, after subtracting each file's
   median value to take care of different sky levels
   """

   median_combine(tmplist,'__tmp__med.fits',zeromedian=True)

   """ Sigma-clip the median file and reject the bad pixels """
   data = pf.getdata('__tmp__med.fits')
   m,s = sigma_clip(data,sigclip)
   bpm = n.zeros(data.shape,dtype=int)
   mask = n.absolute(data - m)>sigbpm*s
   bpm[mask] = 1
   bpmhdu = pf.PrimaryHDU(bpm)
   bpmhdu.header.update('object','Bad Pixel Mask')
   bpmhdu.writeto(outbpm,clobber=True)

   """ Clean up """
   if outsky is None:
      os.remove('__tmp__med.fits')
   else:
      os.rename('__tmp__med.fits',outsky)
   for i in tmplist:
      os.remove(i)

#-----------------------------------------------------------------------

def niri_coadd(outroot):
   """
   Runs the pyraf task imcoadd.  
   *** NB: This works well for small dithers, but not so well for mosaics ***
   *** For mosaics use the tasks below ***
   """

   """
   Set up the input and output file names
   """

   ffname  = '%s_ff.fits' % outroot
   outname = '%s_coadd.fits' % outroot

   """
   Run imcoadd
   """

   iraf.imcoadd('@%s'%ffname,outimage=outname)

#------------------------------------------------------------------------------

def niri_fixpix(infiles, boxsize=11, verbose=True):
   """
   Does two passes of fixing the bad pixels in a set of NIRI images.  This
   two-pass process is necessary to deal with the large blocks of bad pixels
   along the edges of the images.
   """

   """ Loop over input files"""
   if verbose:
      print "Running fixpix"
      print "-------------------------------------------------------------"
   for i in infiles:
      """ Open input files """
      sciname = i
      whtname = sciname.replace('.fits','_wht.fits')
      scihdu = pf.open(sciname,mode='update')
      scidat = scihdu[0].data
      whthdu = pf.open(whtname)
      whtdat = whthdu[0].data
      mask = whtdat==0
      """
      First pass: replace all bad pixels with the median data value.
      """
      datmed = n.median(scidat)
      if verbose:
         print " %s - Median value is %7.2f" % (sciname,datmed)
      scidat[mask] = datmed
      scihdu.flush()
      whthdu.close()
      """
      Second pass: run fixpix_wht to replace the bad pixels with the local
      median value
      """
      fixpix_wht(sciname,whtname,boxsize=boxsize)
   print ""

#------------------------------------------------------------------------------

def split_and_fix_ff(inframes, inprefix='ff', outprefix='fc', badpixfile=None, 
                     fixpix=True, bpmgrow=3):
   """

   Splits the ff*fits files created by calib_2 into the format that will
   be used by all of the following tasks.  
   These ff*fits files are multi-extension files containing 3 images:
     science (HDU1)
     variance (HDU2)
     data-quality (HDU3)
   SExtractor needs these images as separate files, which are produced by
    this function.
   At the same time, fixes the bad pixels in all of the files if requested.

   Inputs:
      inframes:   list of input frames
      inprefix:   prefix for input files.  Default='ff' (flat-fielded)
      outprefix:  prefix for output files. Default='fc' (final calibration)
      badpixfile: previously generated bad pixel file.  If this file is
                  given (default=None), then it is combined with the bad
                  pixel list generated here.
      fixpix:     set to True (default) to correct the bad pixels
      bpmgrow:    amount to flag around the bad pixels. Default=3
   """

   """ Set up for running the task """
   if type(inframes) is str or type(inframes) is int:
      print ""
      print "Single input file"
      tmplist = [inframes,]
   elif type(inframes) is n.ndarray:
      print ""
      print "Input numpy array"
      tmplist = inframes
   elif type(inframes) is list:
      print ""
      print "Input file list"
      tmplist = inframes
   else:
      print ""
      print "Warning.  Input frames need to be either a list of frames "
      print " (python type==list) or a single input frame number."
      print ""
      return

   """ Loop through the files in the list """
   sci_files = []
   texp = []
   ncoadd = []
   for i in range(len(tmplist)):

      """ Open the input file """
      tmpname = '%s%d.fits' % (inprefix,tmplist[i])
      f = pf.open(tmpname)
      print ""
      print "Splitting input file in preparation for running SExtractor"
      print "------------------------------------------------------------"
      f.info()

      """ Split the input file data """
      sciname = '%s%d_sci.fits' % (outprefix,tmplist[i])
      #sciname  = tmplist[i].replace('.fits','_sci.fits')
      whtname  = sciname.replace('.fits','_wht.fits')
      flagname = sciname.replace('sci.fits','flag.fits')
      sci_files.append(sciname)
      scidat   = f[1].data.copy()
      vardat   = f[2].data.copy()
      flagdat  = f[3].data.copy()

      """
      Create the weight file
      For now, just set the weight=1., since there is a artifact in the
      variance images produced by the gemini pipeline
      """
      #vardat[vardat<1.0e-29] = 1.0e-29
      #whtdat   = 1./vardat
      whtdat = n.ones(scidat.shape)

      """ 
      Modify the weight file to reflect the flagged pixels, including
      growing the flagged regions a bit.
      If there is an external bad pixel file, include it before growing
      the flagged regions.
      """
      if badpixfile is not None:
         extflags = pf.getdata(badpixfile)
         allflags = flagdat + extflags
      else:
         allflags = flagdat
      flagdat = filters.maximum_filter(allflags,bpmgrow)
      whtdat[flagdat>0] = 0.

      """ Write out the 3 separate output files """
      hdr = f[0].header.copy()
      print "Writing %s" %sciname
      pf.PrimaryHDU(scidat,hdr).writeto(sciname,clobber=True)
      print "Writing %s" %whtname
      pf.PrimaryHDU(whtdat,hdr).writeto(whtname,clobber=True)
      print "Writing %s" %flagname
      pf.PrimaryHDU(flagdat,hdr).writeto(flagname,clobber=True)

      """ Save info for later steps """
      texp.append(hdr['coaddexp'])
      ncoadd.append(hdr['coadds'])

      """ Clean up """
      del scidat,vardat,whtdat,flagdat,allflags,hdr
      f.close()

   """ Run fixpix if requested """
   if fixpix:
      print ""
      niri_fixpix(sci_files)

#------------------------------------------------------------------------------

def niri_sextractor(infiles, catformat='ascii'):
   """
   Runs SExtractor on the files that were produced by the split_and_fix_ff
    function.
   After running SExtractor, make a ds9 region file from the catalog
   """

   """ Set up for running the task """
   if type(infiles) is str:
      print ""
      print "Single input file"
      sci_files = [infiles,]
   elif type(infiles) is list:
      print ""
      print "Input file list"
      sci_files = infiles
   else:
      print ""
      print "Warning.  Input files need to be either a list of files "
      print " (python type==list) or a single input file name."
      print ""
      return

   """ Run SExtractor on the input list """
   print ""
   for i in range(len(sci_files)):
      f = sci_files[i]
      sciname = f
      whtname = f.replace('.fits','_wht.fits')
      catname = f.replace('.fits','.cat')
      regname = f.replace('.fits','.reg')
      hdr = pf.getheader(f)
      texp = hdr['coaddexp']
      ncoadd = hdr['coadds']

      astrom.make_cat_niri(sciname,catname,regname,catformat=catformat,
                           texp=texp,ncoadd=ncoadd,weight_file=whtname,
                           weight_type='MAP_WEIGHT')

#------------------------------------------------------------------------------

def ccmap_tile(frames, astref, ccdbfile, bpmfile=None, swarplist=None):
   """
   Solves for the relative astrometry for all the frames in a single tile
   using the pyraf ccmap task.  In this process, the astrometry
   of the reference frame is assumed to be correct, and all of the other
   frames in the tile are aligned to it.  The steps are:
     1. Run SExtractor on each frame, using that frame's WCS data
     2. Create a matched catalog, containing the (x,y) coordinates for all
       objects in the frame that have (RA,Dec) coordinates from the reference
       frame.
     3. Use the matched catalog as input for the pyraf ccmap task.

   NOTE: A "tile" consists of a number of dithered positions, where the dithers 
   are small enough so that there is significant overlap between frames
   """

   try:
      os.remove(ccdbfile)
   except:
      pass

   """ 
   Run a new loop, so that the refcat is properly created from the previous
   loop, just in case that the reference frame is not the first frame in the
   list.
   """
   slist = []
   for i in frames:
      """ Set things up for astrometry, skipping the astrometric reference """
      fitscat = 'ff%d.cat' % i
      fitsfile  = fitscat.replace('.cat','_sci.fits')
      ccmapfile = fitscat.replace('.cat','.ccmap')
      slist.append(fitsfile)
      if fitscat == astref:
         print "Skipping over %s since it is the reference catalog" % fitscat
         continue

      """ 
      Match the objects in the input file to those in the astrometric reference 
      file 
      """
      astsimp.match_fits_to_ast(fitsfile,fitscat,astref,ccmapfile,max_offset=40.,
                                doplot=False)

      """ Run ccmap on the matched catalog """
      astsimp.rscale_ccmap(ccmapfile,ccdbfile,fitsfile,interactive=True)

   """ 
   Write an output file listing the *sci.fits files for use in swarp, if
   requested
   """

   if swarplist is not None:
      f = open(swarplist,'w')
      for i in slist:
         f.write('%s\n' % i)
      f.close()

#------------------------------------------------------------------------------

def niri_photom(infiles, photocat, magcol, photzp=30., fixzp=None, 
                mautocol=19, mapermaxcol=49, returnzp=False, 
                weightfiles='default', photconfig='sext_photom.config', 
                racol=1, deccol=2, xcol_fits=8, ycol_fits=9, max_offset=40.,
                verbose=True, summary=True, doplot=False):
   """ 
   Use a photometric catalog (probably from 2MASS) to set the proper zero-point 
   for the input fits file.  This function will run SExtractor on the input fits 
   fileand get the instrumental magnitude of any of the photometric objects in 
   the field of view of the image.  By comparison to the true magnitudes, the
   median zero-point for the image will be determined.
   """

   """ Set up for running the task """
   if type(infiles) is str:
      if verbose:
         print ""
         print "Single input file"
      tmplist = [infiles,]
      tmpwht = [weightfiles,]
   elif type(infiles) is list:
      if verbose:
         print ""
         print "Input file list"
      tmplist = infiles
      if type(weightfiles) is str:
         tmpwht = []
         for i in range(len(infiles)):
            tmpwht.append(weightfiles)
      else:
         tmpwht  = weightfiles
   else:
      print ""
      print "Warning.  Input files need to be either a list of files "
      print " (python type==list) or a single input file name."
      print ""
      return

   """ Get the data from the photometric catalog """
   photdat = n.loadtxt(photocat)

   """ Loop over the input file(s) """
   zpmedind  = n.zeros(len(tmplist))
   count = 0
   for f in tmplist:
      """ Open up the input file and extract information """
      hdu = pf.open(f,mode='update')
      hdr = hdu[0].header
      try:
         texp = hdr['coaddexp']
      except:
         texp = hdr['exptime']
      try:
         ncoadd = hdr['coadds']
      except:
         ncoadd = 1

      """ Set up the weight file and run SExtractor on the input file """
      if tmpwht[count] == 'default':
         weightfile = f.replace('.fits','_wht.fits')
      outcat = f.replace('.fits','_phot.cat')
      print ""
      print "Running SExtractor"
      print "-----------------------------------------------"
      print "Input file:     %s" % f
      print "Weight file:    %s" % weightfile
      print "Output catalog: %s" % outcat
      try:
         astrom.make_cat_niri(f,outcat,configfile=photconfig,texp=texp,
                              ncoadd=ncoadd,weight_file=weightfile,zeropt=photzp,
                              verbose=False)
      except:
         return
      count += 1


   """ Do the photometric comparison between the catalogs """
   if returnzp:
      zpmedind = astrom.do_photom(infiles,photocat,magcol,photzp,'default',
                                  fixzp,mautocol,mapermaxcol,returnzp,
                                  racol,deccol,xcol_fits,ycol_fits,
                                  max_offset,'2MASS',verbose,summary,doplot)
      return zpmedind
   else:
      astrom.do_photom(infiles,photocat,magcol,photzp,'default',
                       fixzp,mautocol,mapermaxcol,returnzp,
                       racol,deccol,xcol_fits,ycol_fits,
                       max_offset,'2MASS',verbose,summary,doplot)
