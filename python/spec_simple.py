"""
spec_simple.py - A library of functions to do various basic CCD spectroscopy
  processing operations

Classes (UNDER CONSTRUCTION, BUT GETTING THERE...)
  Spec1d
  Spec2d

NOTE: More and more of the previous stand-alone functionality has been moved
into either the Spec1d or  Spec2d class.  However, some stand-alone functions
will probably not be moved.  These include the list below.

Stand-alone functions:
   plot_model_sky_ir  - plots the NIR sky transmission as well as the night-sky
                        emission lines to aid in planning NIR spectroscopy
                        observations.
   xxxx               - more descriptions here

"""

import sys, os
import glob
from math import sqrt,pi
try:
   from astropy.io import fits as pf
except:
   import pyfits as pf
import imfuncs as imf
import numpy as np
import scipy as sp
from scipy import optimize,interpolate,ndimage
import matplotlib.pyplot as plt
import ccdredux as ccd
#import nirspec

#-----------------------------------------------------------------------

def clear_all():
   """
   Clears all of the open figures
   """

   for i in plt.get_fignums():
      plt.figure(i)
      plt.clf()

#-----------------------------------------------------------------------

class Spec1d():
   """
   A class to process and analyze 1-dimensional spectra.
   """

   def __init__(self, infile=None, informat='text', 
                wav=None, flux=None, var=None, sky=None, logwav=False):
      """
      Reads in the input 1-dimensional spectrum.
      This can be done in two mutually exclusive ways:

        1. By providing some 1-d arrays containing the wavelength vector,
           the flux vector, and, optionally, the variance vector.

                 ---- or ----

        2. By providing the name of a file that contains the spectrum.
            There are four possible input file formats:
              fits: A multi-extension fits file
                    Extension 1 is the wavelength vector
                    Extension 2 is the extracted spectrum (flux)
                    Extension 3 is the variance spectrum
              fitstab: A binary fits table, stored as a recarray, that has 
                    columns (or "fields" as they are referred to in a recarray
                    structure) corresponding to, at a mininum, wavelength and
                    flux.
                    NOTE. The table is assumed to be in Extension 1 and the
                    table is assumed to have wavelength in field 0 and flux in
                    field 1
              mwa:  A multi-extension fits file with wavelength info in
                    the fits header
                    Extension 1 is the extracted spectrum (flux)
                    Extension 3 is the variance spectrum
              text: An ascii text file with information in columns:
                 Column 1 is the wavelength
                 Column 2 is the extracted spectrum
                 Column 3 (optional) is the variance spectrum
                 Column 4 (optional) is the sky spectrum [NOT YET IMPLEMENTED]
                 
                 Therefore, an input text file could have one of three formats:
                    A.  wavelength flux
                    B.  wavelength flux variance
                    C.  wavelength flux variance sky


      Inputs:
         infile     - Name of the input file.  If infile is None, then the
                       spectrum must be provided via the wavelength and flux
                       vectors.
         informat   - format of input file ('fits', 'fitstab', 'mwa', or 'text').
                       Default = 'text'
         wav        - 1-dimensional array containing "wavelength" information,
                       either as actual wavelength or in pixels
         flux       - 1-dimensional array containing the flux information for
                       the extracted spectrum
         var        - [OPTIONAL] 1-dimensional array containing the variance
                       spectrum.  Remember: rms = sqrt(variance)
         sky        - [OPTIONAL] 1-dimensional array containing the sky spectrum
         logwav     - if True then input wavelength is logarithmic, i.e., the
                       numbers in the input wavelength vector are actually
                       log10(wavelength)
                       if False (the default), then the input wavelength vector
                       is linear.
      """

      """ Initialize some variables """
      self.wav  = None
      self.flux = None
      self.var  = None
      self.sky  = None
      self.varspec = None
      self.infile = None
      self.smoflux = None
      self.smovar  = None

      """ Check inputs """
      self.logwav = logwav
      if infile is not None:
         self.infile = infile
         try:
            self.read_from_file(informat)
         except:
            print ''
            print 'Could not read input file %s' % infile
            print ''
            return
      elif wav is not None and flux is not None:
         if self.logwav:
            self.wav = 10.**wav
         else:
            self.wav = wav
         self.flux = flux
         if var is not None:
            self.var = var
         if sky is not None:
            self.sky = sky
      else:
         print ''
         print 'ERROR: Must provide either:'
         print '  1. A name of an input file containing the spectrum'
         print '  2. At minimum, both of the following:'
         print '       A. a wavelength vector (wav)'
         print '       B. a flux vector (flux)'
         print '     and optionally one or both of the following'
         print '       C. a variance vector (var)'
         print '       D. a sky spectrum vector (sky)'
         print ''
         return

      """ Read in the list that may be used for marking spectral lines """
      self.load_linelist()

   #-----------------------------------------------------------------------

   def read_from_file(self, informat, verbose=True):
      if verbose:
         print ""
         print "Reading spectrum from %s" % self.infile
         print "Input file has format: %s" % informat

      """ Read in the input spectrum """
      if informat=='fits':
         hdu = pf.open(self.infile)
         print hdu[1].data
         if self.logwav:
            self.wav  = 10.**(hdu[1].data)
         else:
            self.wav  = hdu[1].data.copy()
         self.flux = hdu[2].data.copy()
         self.var  = hdu[3].data.copy()
         self.varspec = True
         del hdu
      elif informat=='fitstab':
         hdu = pf.open(self.infile)
         #hdu.info()
         tdat = hdu[1].data
         self.wav  = tdat.field(0)
         self.flux = tdat.field(1)
         if len(tdat[0]) > 2:
            self.var = tdat.field(3)
            self.varspec = True
         del hdu
      elif informat=='mwa':
         hdu = pf.open(self.infile)
         self.flux = hdu[1].data.copy()
         self.var  = hdu[3].data.copy()
         self.varspec = True
         self.wav = np.arange(self.flux.size)
         hdr1 = hdulist[1].header
         if self.logwav:
            self.wav = hdr1['crval1'] + 10.**(self.wav*hdr1['cd1_1'])
         else:
            self.wav = hdr1['crval1'] + self.wav*hdr1['cd1_1']
         del hdu
      else:
         spec = np.loadtxt(self.infile)
         if self.logwav:
            self.wav  = 10.**(spec[:,0])
         else:
            self.wav  = spec[:,0]
         self.flux = spec[:,1]
         if spec.shape[1] > 2:
            self.var = spec[:,2]
         if spec.shape[1] > 3:
            self.sky = spec[:,3]
         del spec

      """ Check for NaN's, which this code can't handle """
      if self.varspec:
         mask = (np.isnan(self.flux)) | (np.isnan(self.var))
         varmax = self.var[~mask].max()
         self.var[mask] = varmax * 5.
      else:
         mask = (np.isnan(self.flux))
      self.flux[mask] = 0
      if verbose:
         print " Spectrum Start: %8.2f" % self.wav[0]
         print " Spectrum End:   %8.2f" % self.wav[-1]
         print " Dispersion (1st pixel): %6.2f" % (self.wav[1]-self.wav[0])
         print " Dispersion (average):   %6.2f" % \
             ((self.wav[-1]-self.wav[0])/(self.wav.size-1))
         print ""

   #-----------------------------------------------------------------------

   def plot(self, xlabel='Wavelength (Angstroms)', ylabel='Relative Flux', 
            title='default', docolor=True, color='b', linestyle='', label=None,
            fontsize=12, rmscolor='r', rmsoffset=0, rmsls=None, 
            add_atm_trans=False, atmscale=1.05, atmfwhm=15., atmoffset=0., 
            atmls='-', usesmooth=False, verbose=True):
      """
      Plots the spectrum
      """

      """ Set the title """
      if title == 'default':
         if self.infile is None:
            title = 'Extracted Spectrum'
         else:
            title = 'Spectrum for %s' % self.infile

      """ Override color assignments if docolor is False"""
      if not docolor:
         color = 'k'
         rmscolor  = 'k'

      """ Draw the flux=0 line"""
      plt.axhline(color='k')

      """ Plot the spectrum """
      if usesmooth and self.smoflux is not None:
         flux = self.smoflux
         var  = self.smovar
      else:
         flux = self.flux
         var  = self.var
      ls = "steps%s" % linestyle
      if label is not None:
         plabel = label
      else:
         plabel = 'Flux'
      plt.plot(self.wav,flux,color,linestyle=ls,label=plabel)
      plt.tick_params(labelsize=fontsize)
      plt.xlabel(xlabel,fontsize=fontsize)

      """ Plot the RMS spectrum if the variance spectrum exists """
      if var is not None:
         rms = np.sqrt(var)+rmsoffset
         if rmsls is None:
            if docolor:
               rlinestyle = 'steps'
            else:
               rlinestyle = 'steps:'
         else:
            rlinestyle = 'steps%s' % rmsls
         if docolor:
            plt.plot(self.wav,rms,rmscolor,linestyle=rlinestyle,label='RMS')
         else:
            plt.plot(self.wav,rms,rmscolor,linestyle=rlinestyle,label='RMS',lw=2)

      """ More plot labels """
      plt.ylabel(ylabel,fontsize=fontsize)
      if(title):
         plt.title(title)
      if(self.wav[0] > self.wav[-1]):
         plt.xlim([self.wav[-1],self.wav[0]])
      else:
         plt.xlim([self.wav[0],self.wav[-1]])
      print self.wav[0],self.wav[-1]

      """ Plot the atmospheric transmission if requested """
      if add_atm_trans:
         plot_atm_trans(self.wav, atmfwhm, self.flux, scale=atmscale, 
                        offset=atmoffset, linestyle=atmls)

   #-----------------------------------------------------------------------

   def smooth_boxcar(self, filtwidth, doplot=True, outfile=None, 
                     color='b', title='default'):
      """
      Does a boxcar smooth of the spectrum.  
      The default is to do inverse variance weighting, using the variance 
       spectrum if it exists.
      The other default is not to write out an output file.  This can be
      changed by setting the outfile parameter.
      """

      """ Set the weighting """
      if self.var is not None:
         print 'Weighting by the inverse variance'
         wht = 1.0 / self.var
      else:
         print 'Uniform weighting'
         wht = 0.0 * self.flux + 1.0

      """ Smooth the spectrum and store results in smoflux and smovar """
      yin = wht * self.flux
      smowht = ndimage.filters.uniform_filter(wht,filtwidth)
      self.smoflux = ndimage.filters.uniform_filter(yin,filtwidth)
      self.smoflux /= smowht
      if self.var is not None:
         self.smovar = 1.0 / (filtwidth * smowht)

      """ Plot the smoothed spectrum if desired """
      if doplot:
         self.plot(usesmooth=True,title=title,color=color)

      """ Save the output file if desired """
      #if(outfile):
      #   print "Saving smoothed spectrum to %s" % outfile
      #   if varwt:
      #      save_spectrum(outfile,wavelength,outflux,outvar)
      #   else:
      #      save_spectrum(outfile,wavelength,outflux)
      #   print ""

   #-----------------------------------------------------------------------

   def load_linelist(self):

      linefmt = [('name','S10'),('wavelength',float),('label','S10'),\
                    ('dxlab',float),('dir',int),('plot',bool)]
      self.lineinfo = np.array([
            ("Ly-alpha",  1216.,    r"Ly $\alpha$",0.0,1,True),
            ("C IV",      1549.,    "C IV",        0.0,1,True),
            ("C III]",    1909.,    "C III]",      0.0,1,True),
            ("Mg II",     2800.,    "Mg II",       0.0,0,True),
            ("[O II]",    3726.03,  "[O II]",      0.0,1,True),
            ("[O II]",    3728.82,  "[O II]",      0.0,1,False),
            ("CN bandhd", 3883,     "CN",          0.0,-1,True),
            ("CaII K",    3933.667, "CaII K",      0.0,-1,True),
            ("CaII H",    3968.472, "CaII H",      0.0,-1,True),
            ("H-delta",   4101,     r"H$\delta$",  0.0,0,True),
            ("G-band",    4305,     "G-band",      0.0,-1,True),
            ("H-gamma",   4340,     r"H$\gamma$",  0.0,0,True),
            ("Fe4383",    4383,     "Fe4383",      0.0,-1,True),
            ("Ca4455",    4455,     "Ca4455",      0.0,-1,True),
            ("Fe4531",    4531,     "Fe4531",      0.0,-1,True),
            ("H-beta",    4861,     r"H$\beta$",   0.0,0,True),
            ("[O III]",   4962.,    "[O III]",     0.0,1,False),
            ("[O III]",   5007.,    "[O III]",     0.0,1,True),
            ("Mg I (b)",  5176,     "Mg b",        0.0,-1,True),
            ("[N I]",     5199.,    "[N I]",       0.0,1,True),
            ("HeI",       5876.,    "He I",        0.0,1,True),
            ("Na I (D)",  5893,     "Na D",        0.0,-1,True),
            ("[O I]",     6300.,    "[O I]",       0.0,1,True),
            ("[N II]",    6548.,    "[N II]",      0.0,1,False),
            ("H-alpha",   6562.8,   r"H$\alpha$",  0.0,0,True),
            ("[N II]",    6583.5,   "[N II]",      0.0,1,False),
            ("[S II]",    6716.4,   "[S II]",      0.0,1,False),
            ("[S II]",    6730.8,   "[S II]",      0.0,1,True),
            ("Ca triplet",8498.03,  "CaII",        0.0,-1,True),
            ("Ca triplet",8542.09,  "CaII",        0.0,-1,True),
            ("Ca triplet",8662.14,  "CaII",        0.0,-1,True),
            ("[S III]",   9069,     "[S III]",     0.0,1,True),
            ("[S III]",   9532,     "[S III]",     0.0,1,True),
            ("Pa-gamma", 10900.,    r"Pa $\gamma$",0.0,1,True),
            ("Pa-beta",  12800.,    r"Pa $\beta$", 0.0,1,True),
            ("Pa-alpha", 18700.,    r"Pa $\alpha$",0.0,1,True)\
            ],dtype=linefmt)

   #-----------------------------------------------------------------------

   def mark_speclines(self, linetype, z, usesmooth=False, marktype='tick', 
                      labww=20., labfs=12, tickfrac=0.05, tickfac=0.75,
                      showz=True, labloc='default', labcolor='k'):
      """
      A generic routine for marking spectral lines in the plotted spectrum.
      The required linetype parameter can be either 'abs' or 'em' and will 
       determine whether absorption or emission lines are marked.

      Inputs:
        linetype - Must be either 'abs' or 'em' to mark absorption or emission
                    lines, respectively
        z        - redshift to be marked
        labww    - width in pixels of the window used to set the vertical
                   location of the tickmark (location is set from the minimum
                   value within the window).
        labfs    - font size for labels, in points
        ticklen  - override of auto-determination of tick length if > 0
      """

      """ Check linetype """
      if linetype == 'abs':
         pm = -1.
         labva = 'top'
      elif linetype == 'em':
         pm = 1.
         labva = 'bottom'
      else:
         print ''
         print "ERROR: linetype must be either 'abs' or 'em'"
         print ''
         return

      """ Set the display limits """
      lammin,lammax = self.wav.min(),self.wav.max()
      if plt.xlim()[0] > lammin: lammin = plt.xlim()[0]
      if plt.xlim()[1] < lammax: lammax = plt.xlim()[1]
      dlam = self.wav[1] - self.wav[0]
      if usesmooth:
         ff = self.smoflux[(self.wav>=plt.xlim()[0]) & (self.wav<=plt.xlim()[1])]
      else:
         ff = self.flux[(self.wav>=plt.xlim()[0]) & (self.wav<=plt.xlim()[1])]
      fluxdiff = ff.max() - ff.min()
      x0,x1 = plt.xlim()
      y0,y1 = plt.ylim()
      xdiff = x1 - x0
      ydiff = y1 - y0
      dlocwin = labww / 2.

      """ Select lines within current display range """
      zlines = (z+1.0) * self.lineinfo['wavelength']
      zmask = np.logical_and(zlines>lammin,zlines<lammax)
      if linetype == 'em':
         tmask = self.lineinfo['dir'] > -1
      else:
         tmask = self.lineinfo['dir'] < 1
      mask = zmask & tmask
      tmplines = self.lineinfo[mask]
      zlines = (z+1.0) * tmplines['wavelength']
      print ""
      print "Line      lambda_obs"
      print "--------  -----------"
      for i in range(len(tmplines)):
         print "%-9s %8.2f" % (tmplines['name'][i],zlines[i])

      """ Set the length of the ticks """
      ticklen = tickfrac * ydiff

      """ Choose whether to use the smoothed flux or not """
      if usesmooth:
         flux = self.smoflux
      else:
         flux = self.flux

      print ''
      if (len(tmplines) > 0):
         xarr = np.zeros(len(tmplines))
         specflux = np.zeros(len(tmplines))
         for i in range(len(tmplines)):
            xarr[i] = tmplines['wavelength'][i]*(z+1.0)
            tmpmask = np.fabs(self.wav-xarr[i])<dlocwin
            if linetype == 'em':
               specflux[i] = flux[tmpmask].max()
            else:
               specflux[i] = flux[tmpmask].min()
         for i in range(len(tmplines)):
            info = tmplines[i]
            tickstart = specflux[i] + pm*tickfac*ticklen
            tickend = tickstart + pm*ticklen
            labstart  = tickstart + pm*1.5*ticklen
            plt.plot([xarr[i],xarr[i]],[tickstart,tickend],'k')
            if info['plot']:
               plt.text(xarr[i]+info['dxlab'],labstart,info['label'],
                        rotation='vertical',ha='center',va=labva,
                        color=labcolor,fontsize=labfs)
            #elif tmpfmin[i] > plt.ylim()[0]:
            #   plt.axvline(xarr[i],linestyle='--',color='k',lw=1)
            #   #print i['label'],tickstart,labstart,ticklen,tmpfmin

      """ Label the plot with the redshift, if requested """
      if showz:
         if labloc=='topright':
            labx = x0 + 0.95*xdiff
            laby = y0 + 0.95*ydiff
            ha = 'right'
         else:
            labx = x0 + 0.05*xdiff
            laby = y0 + 0.95*ydiff
            ha = 'left'
         print labx,laby
         plt.text(labx,laby,'z = %5.3f'%z,ha=ha,va='center',color=labcolor,
                  fontsize=labfs+4)

   #-----------------------------------------------------------------------

   def save(self, outfile, outformat='text'):
      """
      Saves a spectrum as a text file
      """

      if self.var is not None:
         if self.sky is not None:
            outdata = np.zeros((self.wav.shape[0],4))
            fmtstring = '%7.2f %9.3f %10.4f %9.3f'
            outdata[:,3] = self.sky
         else:
            outdata = np.zeros((self.wav.shape[0],3))
            fmtstring = '%7.2f %9.3f %10.4f'
         outdata[:,2] = self.var
      else:
         outdata = np.zeros((self.wav.shape[0],2))
         fmtstring = '%7.2f %9.3f'
      outdata[:,0] = self.wav
      outdata[:,1] = self.flux
      print ""
      print "Saving spectrum to file %s" % outfile
      np.savetxt(outfile,outdata,fmt=fmtstring)
      del outdata

   #-----------------------------------------------------------------------
   #-----------------------------------------------------------------------

#===========================================================================
#
# Start of Spec2d class
#
#===========================================================================

class Spec2d(imf.Image):
   """
   A class to process 2-dimensional spectra, i.e., the CCD data that
    comes out of a typical spectrograph.  
   The main purpose of this Spec2d class and its associated functions is to
    extract a 1-dimensional spectrum from a 2-dimensional spectrum.
    The extracted 1-dimensional spectrum will, in the end, be output into a file
    that can be analyzed using the Spec1d class.
   NOTE: Spec2d inherits the properties of the Image class that is defined
    in imfuncs.py
   """

   def __init__(self, inspec, hext=0, extvar=None,
                xtrim=None, ytrim=None, transpose=False, fixnans=True, 
                logwav=False, verbose=True):
      """
      Reads in the 2-dimensional spectrum from a input fits file (or the
      HDUList from a previously loaded fits file) and
      stores it in a Spec2d class container.

      Required inputs:
         inspec    - The input spectrum.  This can either be: 
                     1. a filename, the most common case
                          - or -
                     2. a HDU list.  It is sometimes the case, e.g., with 
                     multi-extension fits files such as those produced by ESI, 
                     that the fits file has already been loaded.  In this case 
                     it is more efficient to first read the HDU list from the
                     input file (external to this class) and then pass that
                     HDU list and a desired HDU (set by the hext parameter, see
                     below) to Spec2d instead.  For example, the Esi2d class
                     does this.

      Optional inputs:
         hext      - The header-data unit (HDU) that contains the 2-dimensional
                     spectroscopic data.  The default value (hdu=0) should work
                     for most fits files.
         extvar    - If the 2d variance spectrum has already been computed by
                     previous reduction steps and stored as a separate external
                     file, then it needs to be read in.  
                     Default value is None, implying no external variance file.
                     If set, this can either be a filename or a hdulist if the
                     file has already been opened.
         xtrim     - Change from the default value (None) if the input spectrum
                     needs to be trimmed along the x-axis.
                     Example format for trimming:  xtrim=[300,701]
         ytrim     - Change from the default value (None) if the input spectrum
                     needs to be trimmed along the y-axis.
                     Example format for trimming:  ytrim=[300,701]
         transpose - If transpose=True, transpose the x and y dimensions of the
                     input spectrum.  This is done, e.g., to change the 
                     dispersion axis from vertical to horizontal. 
                     NOTE: If transpose=True, the transpose action happens AFTER
                     any trimming that is done if xtrim and/or ytrim have
                     a value different from None
                     Default = False.
         verbose   - Set to True (the default) for information about the 
                     input file to be printed.
      """

      """ Initialize some variables """
      self.dispaxis  = 'x'
      self.specaxis  = 1
      self.spaceaxis = 0
      self.extvar    = None
      self.sky1d     = None
      self.sky2d     = None
      self.skysub    = None
      self.spec1d    = None
      self.fitrange  = None
      self.apmin     = -4.
      self.apmax     = 4.
      self.muorder   = 3
      self.sigorder  = 3
      self.p0        = None  # Parameters describing fit to the spatial profile
      self.logwav    = logwav

      """ 
      Read in the data and call the superclass initialization for useful 
      Image attributes (for now this superclass initialization only works if
      inpsec is a file name.
      """
      test = str(type(inspec))
      if test.rfind('hdu') > 0:
         self.hdu = inspec
      else:
         imf.Image.__init__(self,inspec,verbose=verbose)

      """ Read in the external variance file if there is one """
      if extvar is not None:
         test = str(type(extvar))
         if test.rfind('hdu') > 0:
            self.extvar = extvar
         else:
            self.extvar = pf.open(extvar)

      """ Set the portion of the input spectrum that should be used """
      self.hdr = self.hdu[hext].header
      nx = self.hdr['naxis1']
      ny = self.hdr['naxis2']
      trimmed = False
      if xtrim is not None:
         xmin = xtrim[0]
         xmax = xtrim[1]+1
         trimmed = True
      else:
         xmin = 0
         xmax = nx
      if ytrim is not None:
         ymin = ytrim[0]
         ymax = ytrim[1]+1
         trimmed = True
      else:
         ymin = 0
         ymax = ny
      
      """ Put the data in the appropriate container """
      if transpose:
         self.data = (self.hdu[hext].data[ymin:ymax,xmin:xmax]).transpose()
      else:
         self.data = self.hdu[hext].data[ymin:ymax,xmin:xmax]

      """ 
      Do the same thing for the external variance file, if there is one 
      ASSUMPTION: the external variance file (if it exists) is the same size as 
       the data file and thus should be trimmed, transposed, etc. in the same way
      """
      if self.extvar is not None:
         if transpose:
            self.vardata = \
                (self.extvar[hext].data[ymin:ymax,xmin:xmax]).transpose()
         else:
            self.vardata = self.extvar[hext].data[ymin:ymax,xmin:xmax]

      """ Set other variables and report results """
      self.xmin = xmin
      self.xmax = xmax
      self.ymin = ymin
      self.ymax = ymax
      if verbose:
         print ''
         print '----------------------------------------------------------------'
         print ''
         if self.infile is None:
            print 'Read in 2-dimensional spectrum from HDU=%d' % hext
         else:
            print 'Read in 2-dimensional spectrum from %s (HDU=%d)' % \
                (self.infile,hext)
         if trimmed:
            print 'The input dataset was trimmed'
            print ' xrange: %d:%d.  yrange: %d:%d' % (xmin,xmax,ymin,ymax)
         if transpose:
            print 'The input dataset was transposed'
         print 'Final data dimensions (x y): %d x %d' % \
             (self.data.shape[1],self.data.shape[0])
      self.get_dispaxis(verbose)

      """ 
      Check for NaN's within the spectrum and replace them if they are there
      """
      if fixnans:
         self.fix_nans(verbose=True)

   #-----------------------------------------------------------------------

   def get_dispaxis(self, verbose=True):
      """
      The dispersion axis is the axis corresponding to the spectral direction
      in a 2-dimensional spectrum.  get_dispaxis does the simple task of
      showing the current value of the dispaxis variable, either 'x' or 'y'
      """

      self.npix  = self.data.shape[self.specaxis]
      self.nspat = self.data.shape[self.spaceaxis]
      if verbose:
         print ''
         print 'Current value of dispaxis:              %s' % self.dispaxis
         print 'Number of pixels along dispersion axis: %d' % self.npix
         print ''

   #-----------------------------------------------------------------------

   def set_dispaxis(self, dispaxis):
      """
      The dispersion axis is the axis corresponding to the spectral direction
      in a 2-dimensional spectrum.  
      set_dispaxis is used to change the value of the dispaxis variable.

      For example, if the 2d spectrum was loaded as:
        myspec = Spec2d('myspec.fits')
      then to change the dispersion axis from x to y (the only two choices)
      type: 
        myspec.set_dispaxis('y')

      Required input:
         dispaxis - Dispersion axis: 'x' and 'y' are the only two possible
                    choices
      """

      oldval = self.dispaxis
      if dispaxis == 'x' or dispaxis == 'y':
         self.dispaxis = dispaxis
         if self.dispaxis == "y":
            self.specaxis  = 0
            self.spaceaxis = 1
         else:
            self.specaxis  = 1
            self.spaceaxis = 0
         print ''
         print 'Old value of dispaxis: %s' % oldval
         self.get_dispaxis()
         print ''
         return
      else:
         print ''
         print "ERROR: dispaxis must be either 'x' or 'y'"
         print '%s is not a valid value' % dispaxis
         print ''
         print 'Keeping current value for dispaxis:  %s' % self.dispaxis
         print ''
         return

   #-----------------------------------------------------------------------

   def get_wavelength(self, hext=None, verbose=False):
      """
      Gets a wavelength vector from the fits header, if it exists

      Inputs:
         hext - By default the wavelength information will be searched
                for in the the header that is already associated with
                the data, which was set when reading inp.  This
                process is followed for the default: hext=None 
                If hext is not None, search for the wavelength info 
                in the indicated HDU instead
      """

      if self.dispaxis == 'y':
         dim = 2
      else:
         dim = 1
      cdkey = 'cd%d_%d' % (dim,dim)
      crpix = 'crpix%d' % dim
      crval = 'crval%d' % dim
      if hext is not None:
         hdr = self.hdu[hext].header
      else:
         hdr = self.hdr
      #print cdkey,crpix,crval
      self.has_cdmatx = True
      try:
         dw = hdr[cdkey]
      except:
         self.has_cdmatx = False
         dw = 1
      try:
         wstart = hdr[crval]
      except:
         self.has_cdmatx = False
         wstart = 0
      try:
         wpix = hdr[crpix] - self.xmin - 1
      except:
         self.has_cdmatx = False
         wpix = 0

      """ Create the wavelength vector from the information above """
      self.wavelength = wstart + (np.arange(self.npix) - wpix) * dw
      if self.logwav:
         self.wavelength = 10.**self.wavelength
      if verbose:
         print self.wavelength

   #-----------------------------------------------------------------------

   def fix_nans(self, verbose=False):
      """
      Detects NaN's within the spectrum and replaces them with real numbers
      if they are there.
      """

      nanmask = np.isnan(self.data)
      nnan = nanmask.sum()
      if nnan>0:
         if verbose:
            print 'Found %d NaNs in the two-dimensional spectrum' % nnan

         """ First replace the NaNs with a temporary value """
         self.data[nanmask] = -999

         """ 
         Now find the median sky values by calling the subtract_sky_2d method
         """
         self.subtract_sky_2d()

         """ 
         Finally, replace the NaNs with the median sky for their row/column
         """
         self.data[nanmask] = self.sky2d[nanmask]

   #-----------------------------------------------------------------------

   def display_spec(self, show_skysub=True):
      """
      Displays the two-dimensional spectrum and also, by default, the
      same spectrum after a crude sky subtraction.  To show only the
      input spectrum without the additional sky-subtracted version,
      just set show_skysub=False

      Optional inputs:
         show_skysub - If this is True (the default) then make a second
                       plot, showing the 2-D spectrum after a crude
                       sky-subtraction has been performed.
      """

      if show_skysub:
         """ Subtract the sky if this has not already been done """
         if self.skysub is None:
            self.subtract_sky_2d()

         """ Plot the input spectrum """
         ### NOT DONE YET ###


   #-----------------------------------------------------------------------

   def subtract_sky_2d(self, outfile=None, outsky=None):
      """
      Given the input 2D spectrum, creates a median sky and then subtracts
      it from the input data.  Two outputs are saved: the 2D sky-subtracted
      data and a 1D sky spectrum.

      Optional inputs:
      data       - array containing the 2D spectrum
      outfile    - name for output fits file containing sky-subtracted spectrum
      outskyspec - name for output 1D sky spectrum
      """

      """ Set the dispersion axis direction """
      if self.dispaxis == 'y':
         spaceaxis = 1
      else:
         spaceaxis = 0
   
      """ Take the median along the spatial direction to estimate the sky """
      if self.data.ndim < 2:
         print ''
         print 'ERROR: subtract_sky needs a 2 dimensional data set'
         return
      else:
         self.sky1d = np.median(self.data,axis=spaceaxis)

      """ Turn the 1-dimension sky spectrum into a 2-dimensional form """
      self.sky2d = np.tile(self.sky1d,(self.data.shape[spaceaxis],1))

      """ Subtract the sky from the data """
      self.skysub = self.data - self.sky2d

      ### NOT DONE YET (needs possible saving of sky spectra) ###

   #-----------------------------------------------------------------------

   def spatial_profile(self, pixrange=None, showplot=True, do_subplot=False,
                       title='Spatial Profile', fit=None, normalize=False):
      """
      Compresses a 2d spectrum along the dispersion axis to create a spatial
       profile, and then plots it if requested
      """

      """ Set the data range in which to find the trace """
      if pixrange is not None:
         if self.data.ndim < 2:
            tmpdat = self.data[pixrange[0]:pixrange[1]]
         else:
            if self.specaxis == 0:
               tmpdat = self.data[pixrange[0]:pixrange[1],:]
            else:
               tmpdat = self.data[:,pixrange[0]:pixrange[1]]
         #print pixrange
      else:
         tmpdat = self.data.copy()
         

      """ Compress the data along the dispersion axis and find the max value """
      if self.data.ndim < 2:
         self.cdat = tmpdat
      else:
         self.cdat = np.median(tmpdat,axis=self.specaxis)
      if normalize:
         cmax = self.cdat.max()
         self.cdat /= cmax
         print cmax
      self.x = np.arange(self.cdat.shape[0])

      """ 
      Plot the compressed spectrum, showing the best-fit Gaussian if requested
      """
      if(showplot):
         plt.plot(self.x,self.cdat,linestyle='steps')
         plt.xlabel('Spatial direction (0-indexed)')
         if fit is not None:
            xmod = np.arange(1,self.cdat.shape[0]+1,0.1)
            p = np.atleast_1d(fit) # Make sure parameters are in a numpy array
            ngauss = int((p.size-1)/3)
            if p.size - (ngauss*3+1) != 0:
               print 'Fit not plotted since fit has wrong number of parameters'
               pass
            if normalize:
               p[0] /= cmax
               for j in range(ngauss):
                  p[j*3+3] /= cmax
            ymod = make_gauss(xmod,p)
            plt.plot(xmod,ymod)
            plt.axvline(fit[1]+self.apmin,color='k')
            plt.axvline(fit[1]+self.apmax,color='k')
         if title is not None:
            plt.title(title)

   #-----------------------------------------------------------------------

   def locate_trace(self, pixrange=None, init=None, fix=None, 
                    showplot=True, do_subplot=False, ngauss=1,
                    verbose=True):
      """
      Compresses a 2d spectrum along the dispersion axis so that
       the trace of the spectrum can be automatically located by fitting
       a gaussian + background to the spatial direction.  The function
       returns the parameters of the best-fit gaussian.
      The default dispersion axis is along the x direction.  To change this
       set the dispaxis to "y" with the set_dispaxis method in this Spec2d class.
      """

      """ Start by compressing the data, but don't show it yet """
      self.spatial_profile(pixrange,showplot=False)

      """ 
      Set up the container for the initial guesses for the parameter values
      """
      nparam = 3*ngauss + 1
      p_init = np.zeros(nparam)

      """
      Put default choices, which may be overridden, into p_init
      """
      p_init[0] = np.median(self.data,axis=None)
      for i in range(ngauss):
         """ 
         In this loop the parameters are set as follows:
           p_init[ind] is either mu (if i==0) or an offset from p_init[1]
           p_init[ind+1] is sigma
           p_init[ind+2] is amplitude
         """
         ind = 3*i + 1
         if i==0:
            tmp = self.cdat.argsort()
            p_init[ind] = 1.0 * tmp[tmp.shape[0]-1]
         else:
            p_init[ind] = 5. * i * (-1.)**(i+1)
         p_init[ind+1] = 3.
         p_init[ind+2] = self.cdat.max() - p_init[0]

      """
      Override the default values if init has been set.
      NOTE: the init parameter must be an array (or list) of length nparam
       otherwise the method will quit
      NOTE: A value of -999 in init means keep the default value
      """
      if init is not None:
         if len(init) != nparam:
            print ''
            print 'ERROR: locate_trace -- init parameter must have length %d' \
                % nparam
            print '  (since ngauss=%d ==> nparam= 3*%d +1)' % (ngauss,ngauss)
            print ''
            return np.nan
         for j in range(nparam):
            if init[j]>-998.:
               p_init[j] = init[j]

      """ 
      Set up which parameters are fixed based on the fix parameter
      The default value (fix=None) means that all of the parameters are varied
       in the fitting process
      """
      fixstr = np.zeros(nparam,dtype='S3')
      if fix is None:
         fixvec = np.zeros(nparam)
      else:
         fixvec = np.atleast_1d(fix)
         if fixvec.size != nparam:
            print ''
            print 'ERROR: locate_trace -- fix parameter must have length %d' \
                % nparam
            print '  (since ngauss=%d ==> nparam= 3*%d +1)' % (ngauss,ngauss)
            print ''
            return np.nan
         fixstr[fixvec==1] = 'Yes'
      fitmask = fixvec==0
      fitind  = np.arange(nparam)[fitmask]

      """ Fit a Gaussian plus a background to the compressed spectrum """
      mf=100000
      p = p_init[fitmask]
      p_fit,ier = optimize.leastsq(fit_gauss,p,
                                      (self.x,self.cdat,p_init,fitind),
                                      maxfev=mf)
      """ 
      Create the full parameter list for the fit by combining the fitted
      parameters and the fixed parameters
      """
      p_out = p_init.copy()
      p_out[fitind] = p_fit

      """ Give results """
      if(verbose):
         print ""
         print "Profile fit results"
         print "-------------------"
         print '                              Held'
         print 'Parameter        Init Value  fixed? Final Value'
         print '--------------   ----------  ------ -----------'
         print "background       %9.3f    %3s    %9.3f"     \
             % (p_init[0],fixstr[0],p_out[0])
         for i in range(ngauss):
            ind = 3*i + 1
            j=i+1
            if i==0:
               mustr = 'mu_1'
            else:
               mustr = 'offset_%d' % j
            print "%-9s        %9.3f    %3s    %9.3f" \
                % (mustr,p_init[ind],fixstr[ind],p_out[ind])
            print "sigma_%d          %9.3f    %3s    %9.3f" \
                % (j,p_init[ind+1],fixstr[ind+1],p_out[ind+1])
            print "amp_%d            %9.3f    %3s    %9.3f" \
                % (j,p_init[ind+2],fixstr[ind+2],p_out[ind+2])
         print ""

      """ Now plot the spatial profile, showing the best fit """
      if showplot:
         if(do_subplot):
            plt.subplot(221)
         else:
            plt.figure(1)
            plt.clf()
         self.spatial_profile(pixrange,showplot,do_subplot,fit=p_out)

      """ 
      Return the relevant parameters of the fit 
        p_out[0] = bkgd, the background level
        p_out[1] = mu, the location of the peak of the fit to the trace
        p_out[2] = sigma, the width of the fit to the trace
        p_out[3] = amp, the amplitude of the fit
      """
      return p_out

   #-----------------------------------------------------------------------

   def fit_poly_to_trace(self, x, data, fitorder, data0, fitrange=None,
                         doplot=True, markformat='bo', 
                         ylabel='Centroid of Trace (0-indexed)', 
                         title='Location of the Peak'):

      """ Do a sigma clipping to reject clear outliers """
      if fitrange is None:
         tmpfitdat = data
      else:
         fitmask = np.logical_and(x>=fitrange[0],x<fitrange[1])
         tmpfitdat = data[fitmask]
      dmu,dsig = ccd.sigma_clip(tmpfitdat,3.0)
      goodmask = np.absolute(data - dmu)<3.0*dsig
      badmask  = np.absolute(data - dmu)>=3.0*dsig
      dgood    = data[goodmask]
      dbad     = data[badmask]
      xgood    = x[goodmask]
      xbad     = x[badmask]

      """ Fit a polynomial to the trace """

      if fitrange is None:
         xpoly = xgood
         dpoly = dgood
      else:
         fitmask = np.logical_and(xgood>=fitrange[0],xgood<fitrange[1])
         xpoly  = xgood[fitmask]
         dpoly  = dgood[fitmask]

      if fitorder == -1:
         polyorder = 0
      else:
         polyorder = fitorder
      dpoly = np.polyfit(xpoly,dpoly,polyorder)

      if fitorder == -1:
         dpoly[0] = data0

      """ Calculate the fitted function """
      fitx  = np.arange(self.npix).astype(np.float32)
      fity  = 0.0 * fitx
      for i in range(dpoly.size):
         fity += dpoly[i] * fitx**(dpoly.size - 1 - i)

      """ Plot the results """
      ymin = dmu - 4.5*dsig
      ymax = dmu + 4.5*dsig
      if doplot:
         plt.plot(x,data,markformat)
         plt.xlabel("Pixel (dispersion direction)")
         plt.ylabel(ylabel)
         plt.title(title)

         """ Show the value from the compressed spatial profile """
         plt.axhline(data0,color='k',linestyle='--')

         """ Mark the bad points that were not included in the fit """
         plt.plot(xbad,dbad,"rx",markersize=10,markeredgewidth=2)

         """ Show the fitted function """
         plt.plot(fitx,fity,"r")

         """
         Show the range of points included in the fit, if fitrange was set
         """
         if fitrange is not None:
            plt.axvline(fitrange[0],color='k',linestyle=':')
            plt.axvline(fitrange[1],color='k',linestyle=':')
            xtmp = 0.5 * (fitrange[1] + fitrange[0])
            xerr = xtmp - fitrange[0]
            ytmp = fity.min() - 0.2 * fity.min()
            plt.errorbar(xtmp,ytmp,xerr=xerr,ecolor="g",capsize=10)
         plt.xlim(0,self.npix)
         plt.ylim(ymin, ymax)

      """ Return the parameters produced by the fit and the fitted function """
      print dpoly
      return dpoly,fity


   #-----------------------------------------------------------------------

   def trace_spectrum(self, ngauss=1, stepsize=25, muorder=3, sigorder=4, 
                      fitrange=None, doplot=True, do_subplot=False, 
                      verbose=True):
      """
      Fits a gaussian plus background to the spatial profile of the spectrum
       This is done in binned segments, because for nearly all cases the SNR
       in a single column (assuming that y is the spatial axis) is much too low
       to get a good fit to the profile.  The bin size, in pixels in the 
       dispersion direction, is set by the stepsize parameter (default is 25).
      The steps in this method are as follow:
       1. Obtain the parameters of the spatial fit in each bin and save them
          in the mustep and sigstep arrays
       2. Under the assumption that the parameters that describe the spatial 
          profile vary slowly with wavelength, fit a polynomial to the values
          in the mustep and sigstep arrays.  
          This polynomial will then provide the profile-description parameters
          for each individual column (or row, depending on the dispersion
          direction) in the spectrum.
          The order of the polynomials are set by the muorder and sigorder
          parameters. NOTE: values of -1 for these parameters mean to just
          take the values from the overall spatial profile and to not do
          this fitting exercise.
      """

      """
      As a first step, see if either muorder or sigorder are set to -1.
      If that is the case, then we can skip the fitting entirely for
      that parameter
      """
      fitmu  = True
      fitsig = True
      if muorder == -1:
         self.mu = np.ones(self.npix) * self.p0[1]
         fitmu = False
      if sigorder == -1:
         self.sig = np.ones(self.npix) * self.p0[2]
         fitsig = False

      """
      Define the slices through the 2D spectrum that will be used to find
       the centroid and width of the object spectrum as it is traced down 
       the chip
      """
      xstep = np.arange(0,self.npix-stepsize,stepsize)

      """ Set up containers for mu and sigma along the trace """
      mustep  = np.zeros((xstep.size,ngauss))
      sigstep = mustep.copy()
      nsteps  = np.arange(xstep.shape[0])


      if fitmu or fitsig:
         """ Step through the data """
         print ''
         print "Running fit_trace"
         print "--------------------------------------------------------------- "
         print "Finding the location and width of the trace at %d segments" % \
             nsteps.shape[0]
         print "   of the 2D spectrum..."
         for i in nsteps:
            pixrange = [xstep[i],xstep[i]+stepsize]
            p = self.locate_trace(pixrange=pixrange,showplot=False,verbose=False)
            for j in range(ngauss):
               mustep[i,j]  = p[j*3+1]
               sigstep[i,j] = p[j*3+2]
         print "   Done"

      """ Fit a polynomial to the location of the trace """
      if fitmu:
         if doplot:
            if(do_subplot):
               plt.subplot(222)
            else:
               plt.figure(2)
               plt.clf()
         print "Fitting a polynomial of order %d to the location of the trace" \
             % muorder
         self.mupoly,self.mu = \
             self.fit_poly_to_trace(xstep,mustep[:,0],muorder,self.p0[1],
                                    fitrange,doplot=doplot)

      """ Fit a polynomial to the width of the trace """
      if fitsig:
         if doplot:
            if(do_subplot):
               plt.subplot(223)
            else:
               plt.figure(3)
               plt.clf()
         print "Fitting a polynomial of order %d to the width of the trace" \
             % sigorder
         self.sigpoly,self.sig = \
             self.fit_poly_to_trace(xstep,sigstep[:,0],sigorder,self.p0[2],
                                    fitrange,markformat='go',
                                    title='Width of Peak',
                                    ylabel='Width of trace (Gaussian sigma)',
                                    doplot=doplot)
      
   #-----------------------------------------------------------------------

   def find_and_trace(self, ngauss=1, stepsize=25, muorder=3, sigorder=4, 
                      fitrange=None, doplot=True, do_subplot=True, verbose=True):

      """
      The first step in the spectroscopy reduction process.

      The find_and_trace function will:
        1. Locate roughly where the target object is in the spatial direction 
           (usually the y axis is the spatial direction) by taking the
           median in the spectral direction so the peak in the spatial
           direction stands out.  This step provides the initial guesses
           for the location (mu0) and width (sig0) of the peak that are
           then used in the second step.

           * This step is done by a call to the locate_trace method

        2. Once the rough location of the peak has been found, determines how
           its location and width change in the spectral direction.
           That is, this will find where the peak falls in each column.
           It returns the position (pos) and width (width) of the peak as
           a function of x location

           * This step is done by a call to the trace_spectrum method

      Inputs:
         stepsize
         muorder
         sigorder
      """

      self.p0 = self.locate_trace(showplot=doplot,do_subplot=do_subplot,
                                  ngauss=ngauss,pixrange=fitrange,
                                  verbose=verbose)

      self.trace_spectrum(ngauss,stepsize,muorder,sigorder,fitrange,doplot,
                          do_subplot,verbose=verbose)

   #-----------------------------------------------------------------------

   def extract_old(self, weight='gauss', sky=None, gain=1.0, rdnoise=0.0):
      """
      This method contains the old code that was used to do the extraction
      of a 1-D spectrum from the 2-D spectrum in the extract_spectrum
      method.  It has been split out from that while a new method that
      is more compact and may be faster is developed.  Once the new
      method has been written and is shown to produce the same results
      as this old code, the old code will be deleted.
      By splitting this code out in this way and allowing it to be
      called from the extract_spectrum method, we can make the changeover
      from the old code to the new code transparent to the user.
      """

      """ Set the wavelength axis """
      pix = np.arange(self.npix)
      
      """
      Set up the containers for the amplitude and variance of the spectrum
      along the trace
      """
      amp = 0.0 * pix
      var = 0.0 * pix
      skyspec = 0.0 * pix
      tmpy = 0.0 * pix
      tmpyy = 0.0 * pix

      """ Extract the spectrum """
      print ""
      print "Extracting the spectrum..."
      apmin = self.apmin
      apmax = self.apmax
      for i in pix:
         if self.specaxis == 0:
            tmpdata = self.data[i,:]
         else:
            tmpdata = self.data[:,i]
         if (sky == None):
            skyval = None
         else:
            skyval = sky[i]
      
         amp[i],var[i],skyspec[i] = \
             extract_wtsum_col(tmpdata,self.mu[i],apmin,apmax,sig=self.sig[i],
                               gain=gain,rdnoise=rdnoise,sky=skyval,
                               weight=weight)
      print "   Done"

      """ Get the wavelength/pixel vector """
      self.get_wavelength()

      """ Create a Spec1d container for the extracted spectrum """
      if sky is not None:
         skyspec = sky
      self.spec1d = Spec1d(wav=self.wavelength,flux=amp,var=var,sky=skyspec)

   #-----------------------------------------------------------------------

   def extract_new(self, gain=1.0, rdnoise=0.0):
      """

      STILL TO DO:
        - implement uniform weighting - OR - pass the method a previously
          generated profile instead of generating it here

      Extracts a 1d spectrum from the 2d spectrum by doing a weighted
      sum across the spectrum in the spatial direction at each wavelength.
      The extraction has now been set up to follow, at least partially,
      the "optimal extraction" scheme developed by Horne (1986, PASP, 98, 609).
      See more discussion of this scheme below the description of the
      weighting.

      There are three components to the weighting:
        1. The profile of the trace, P, i.e., aperture weighting (for now only 
           uniform or a single gaussian are implemented).  In future, this may 
           be provided in terms of an already constructed profile image rather 
           than calculated within this method.
        2. The aperture definition itself (stored in the apmin and apmax
           variables that are part of the Spec2d class definition).
           This weighting is, in fact, not really a weighting but just a mask 
           set up so that a pixel will get a weight of 1.0 if it is inside the 
           aperture and 0.0 if it is outside.
        3. The statistical errors associated with the detector, etc.,
           in the form of inverse variance weighting.
           The variance can either be provided as an external variance image,
            if the previous reduction / processing has provided this.  
           If no external variance spectrum is provided, then the variance image
            will be constructed from the data counts (including counts from
            a 2d sky spectrum if the sky has already been subtracted from
            the data) plus the gain and readnoise of the detector.

      According to the Horne paper, the optimal extraction of a spectrum
      that has a profile P and proper knowledge of the noise/variance associated
      with each pixel is as follows.  Below D represents the calibrated data,
      S is the sky, V is the pixel variance (based on counts and the detector
      gain and readnoise):

               Sum{ P * (D - S) / V}
           f = ---------------------
                  Sum{ P^2 / V}

      and the variance on the extracted flux, sigma_f^2, is

                      Sum{ P }
        sigma_f^2 = -------------
                    Sum{ P^2 / V}

      NOTE: P must be normalized for each wavelength
      """

      """ 
      Set up arrays for coordinate system 
      Remember, self.npix is the number of pixels in the spectral direction
       and self.nspat is the number of pixels in the spatial direction
      """
      x1d = np.arange(self.npix)
      y1d = np.arange(self.nspat)
      x,y = np.meshgrid(x1d,y1d)
      y = y.astype(float)

      """
      First weighting: Aperture profile
      ---------------------------------
      NOTE: in future, this may be moved into the find_and_trace code

      Make the 1d mu and sig polynomials into 2d polynomials that vary
      along the spectral direction but are identical along a given column.
      The transpose (.T) at the end is necessary because doing a np.repeat
      directly on the desired shape does not give the proper behavior
      (constant along columns in the right way).

      Code below is just for gaussian weighting.
      """
      self.mu2d = self.mu.repeat(self.nspat).reshape((self.npix,self.nspat)).T
      self.sig2d = self.sig.repeat(self.nspat).reshape((self.npix,self.nspat)).T
      ydiff = 1.0*y - self.mu2d
      P = (1./(self.sig2d * sqrt(2.*pi))) * np.exp(-0.5 * (ydiff/self.sig2d)**2)

      """ Make sure the profile is normalized in the spatial direction """
      Pnorm = (P.sum(axis=self.spaceaxis))
      Pnorm = Pnorm.repeat(self.nspat).reshape((self.npix,self.nspat)).T
      self.profile = P / Pnorm

      """
      Second "weighting"/mask: Aperture limits
      ----------------------------------------
      Put in the aperture limits, delimited by apmin and apmax
      """
      apmask = (ydiff>self.apmin-1) & (ydiff<self.apmax)
      bkgdmask = np.logical_not(apmask)

      """ 
      Third weighting: Inverse variance
      ---------------------------------
      Set up the variance based on the detector characteristics
       (gain and readnoise) if an external variance was not provided 
      """
      if self.extvar is not None:
         varspec = self.vardata
      else:
         varspec = (gain * self.data + rdnoise**2) / gain**2

      """ Check for NaNs """
      nansci  = np.isnan(self.data)
      nanvar  = np.isnan(varspec)
      nanmask = np.logical_or(np.isnan(self.data),np.isnan(varspec))
      nnans = nansci.sum()
      nnanv = nanvar.sum()
      nnan = nanmask.sum()

      """ 
      Set up a 2d background grid (think about doing this as part of a call
      to the sky subtraction routine in the future)
      """
      tmp = self.data.copy()
      tmp[apmask] = np.nan
      bkgd = np.nanmedian(tmp,axis=self.spaceaxis)
      bkgd2d = bkgd.repeat(self.nspat).reshape((self.npix,self.nspat)).T
      del tmp

      """ 
      Create the total weight array, combining (1) the aperture profile,
      (2) the aperture limits, and (3) the inverse variance weighting
      following the optimal extraction approach of Horne (1986) as described
      above.
      """
      self.extwt = np.zeros(self.data.shape) 
      vmask = varspec <= 0.
      varspec[vmask] = 1.e8
      varspec[nanvar] = 1.e8
      self.extwt[apmask]  = self.profile[apmask] / (varspec[apmask])
      self.extwt[nanmask] = 0.
      self.extwt[vmask]   = 0.
      wtdenom = (self.profile * self.extwt).sum(axis=self.spaceaxis)
      #wtdenom *= apmask.sum(axis=self.spaceaxis)

      """ Compute the weighted sum of the flux """
      data = self.data
      data[nansci] = 0.
      wtdenom[wtdenom==0] = 1.e9
      flux = \
          ((self.data - bkgd2d) * self.extwt).sum(axis=self.spaceaxis) / wtdenom

      """ 
      Compute the proper variance.
      """
      var = self.profile.sum(axis=self.spaceaxis) / wtdenom

      """ 
      Fix any remaining NaNs (there shouldn't be any, but put this in just to
       be on the safe side
      """
      nansci = np.isnan(flux)
      nanvar = np.isnan(var)
      flux[nansci] = 0.
      var[nansci] = 1.e9
      var[nanvar] = 1.e9

      """ Get the wavelength/pixel vector """
      self.get_wavelength()

      """
      Save the result as a Spec1d instance
      """
      #print '*** Number of nans: %d %d %d ***' % (nnans,nnanv,nnan)
      print ''
      self.spec1d = Spec1d(wav=self.wavelength,flux=flux,var=var,sky=bkgd)
      self.apmask = apmask

   #-----------------------------------------------------------------------

   def extract_spectrum(self, weight='gauss', sky=None, gain=1.0, rdnoise=0.0, 
                        doplot=True, do_subplot=True, outfile=None):
      """
      Second step in reduction process.

      This function extracts a 1D spectrum from the input 2D spectrum
      It uses the information about the trace profile that has been generated 
      by the trace_spectrum function and which is stored (for now) in the
      self.mu and self.sig arrays.
      """

      """ Extract the spectrum """
      #self.extract_old(weight,sky,gain,rdnoise)
      self.extract_new(gain,rdnoise)

      """ Plot the extracted spectrum if desired """
      if doplot:
         print ""
         print "Plotting the spectrum"
         if(do_subplot):
            plt.subplot(224)
         else:
            plt.figure(4)
            plt.clf()
         if self.has_cdmatx:
            xlab = 'Wavelength'
         else:
            xlab = 'Pixel number along the %s axis' % self.dispaxis
         self.spec1d.plot(xlabel=xlab,title='Extracted spectrum')

      """ Save the extracted spectrum to a file if requested """
      if outfile is not None:
         self.spec1d.save(outfile)

   #-----------------------------------------------------------------------


#===========================================================================
#
# *** End of Spec2d class definition ***
#
# Below here are:
#   1. a few stand-alone functions that do not fit particularly well into
#      either the Spec1d or Spec2d class and so will be maintained in
#      their stand-alone format
#   2. mostly old functions that will be kept for legacy purposes but will
#      slowly be modified to function in exactly the way but via calls to
#      the Spec2d and Spec1d classes.  This will be done in such a way to
#      be transparent to the user.
#   3. code that really should be incorporated into the Spec1d or Spec2d
#      classes.  This code will slowly disappear as it gets integrated
#      into the new classes.
#
#===========================================================================

#-----------------------------------------------------------------------

def load_2d_spectrum(filename, hdu=0):
   """
   Reads in a raw 2D spectrum from a fits file

   NOTE: This has now been replaced by the Spec2d class.
   """
   data = Spec2d(filename,hext=hdu)
   print ''
   print 'Loaded a 2-dimensional spectrum from %s' % filename
   print 'Data dimensions (x y): %dx%d' % (data.shape[1],data.shape[0])
   return data

#-----------------------------------------------------------------------

def zap_cosmic_rays(data, outfile, sigmax=5., boxsize=7, dispaxis="x"):
   """
   A method to reject cosmic rays from a 2D spectrum via the following steps:
     1. Creates the median sky from the spectrum
     2. Subtracts the sky from the spectrum
     3. Divides the subtracted spectrum by the square root of the sky, which
        gives it a constant rms
     4. Rejects pixels that exceed a certain number of sigma
     5. Adds the sky back in
   """

   from scipy.ndimage import filters

   """ Set the dispersion axis direction """
   if dispaxis == "y":
      spaceaxis = 1
   else:
      spaceaxis = 0
   
   """ Take the median along the spatial direction to estimate the sky """
   if data.ndim < 2:
      print ""
      print "ERROR: zap_cosmic_rays needs a 2 dimensional data set"
      return
   else:
      sky1d = np.median(data,axis=spaceaxis)
   sky = np.zeros(data.shape)
   for i in range(data.shape[spaceaxis]):
      if spaceaxis == 1:
         sky[:,i] = sky1d
      else:
         sky[i,:] = sky1d

   """ Subtract the sky  """
   skysub = data - sky

   """ 
   Divide the result by the square root of the sky to get a rms image
   """
   ssrms = skysub / np.sqrt(sky)
   m,s = ccd.sigma_clip(ssrms)
   
   """ Now subtract a median-filtered version of the spectrum """
   tmpsub = ssrms - filters.median_filter(ssrms,boxsize)

   """ 
   Make a bad pixel mask that contains pixels in tmpsub with values > sigmax*s
   """
   mask = tmpsub > sigmax * s
   tmpsub[mask] = m

   """ Replace the bad pixels in skysub with a median-filtered value """
   m2,s2 = ccd.sigma_clip(skysub)
   skysub[mask] = m2
   ssfilt = filters.median_filter(skysub,boxsize)
   skysub[mask] = ssfilt[mask]

   """ Add the sky back in and save the final result """
   szapped = skysub + sky
   pf.PrimaryHDU(szapped).writeto(outfile)
   print ' Wrote szapped data to %s' % outfile

   """ Clean up """
   del sky,skysub,ssrms,tmpsub,szapped

#-----------------------------------------------------------------------

def find_blank_columns(data,comp_axis=0,output_dims=1,findblank=False):
   """
   Takes 2-dimensional data and outputs indices of columns not entirely
   composed of zeros. If output_dims is 1, only the indices of the columns
   are give. If output_dims = 2, the indices for every point in any of these
   columns are given. By default, it is actually non-blank columns that are
   found. Setting findblank to True switches this.
   """
   if data.ndim != 2:
      sys.exit("find_blank_columns takes only 2-dimensional data")
   if output_dims==1:
      fbc_tmp = np.zeros(np.shape(data)[1-comp_axis])
      if comp_axis == 0:
         gprelim = np.where(data[int(np.shape(data)[comp_axis]/2),:] == 0)[0]
         for ifbc in range(0,len(gprelim)):
            if len(data[data[:,gprelim[ifbc]] != 0]) == 0:
               fbc_tmp[gprelim[ifbc]] = 1
      else:
         gprelim = np.where(data[:,int(np.shape(data)[comp_axis]/2)] == 0)[0]
         for ifbc in range(0,len(gprelim)):
            if len(data[data[gprelim[ifbc],:] != 0]) == 0:
               fbc_tmp[gprelim[ifbc]] = 1
      if findblank: fbc_tmp = 1-fbc_tmp
      gfbc = np.where(fbc_tmp == 0)[0]
   elif output_dims==2:
      fbc_tmp = np.zeros(np.shape(data))
      if comp_axis == 0:
         gprelim = np.where(data[int(np.shape(data)[comp_axis]/2),:] == 0)[0]
         for ifbc in range(0,len(gprelim)):
            if len(data[data[:,gprelim[ifbc]] != 0]) == 0:
               fbc_tmp[:,gprelim[ifbc]] = 1
      else:
         gprelim = np.where(data[:,int(np.shape(data)[comp_axis]/2)] == 0)[0]
         for ifbc in range(0,len(gprelim)):
            if len(data[data[gprelim[ifbc],:] != 0]) == 0:
               fbc_tmp[gprelim[ifbc],:] = 1
      if findblank: fbc_tmp = 1-fbc_tmp
      gfbc = np.where(fbc_tmp == 0)
   else:
      sys.exit("output_dims parameter for find_blank_columns must be either 1 or 2. Value was: " + str(output_dims))
   return gfbc


#-----------------------------------------------------------------------

def read_spectrum(filename, informat='text', varspec=True, verbose=True):
   """
   Reads in an extracted 1D spectrum and possibly the associated variance
   spectrum.  There are two possible input file formats:
      mwa:  A multi-extension fits file with wavelength info in the fits header
            Extension 1 is the extracted spectrum (flux)
            Extension 3 is the variance spectrum
      text: An ascii text file with information in columns:
            Column 1 is the wavelength
            Column 2 is the extracted spectrum
            Column 3 (optional) is the variance spectrum

   Inputs:
      filename - input file name
      informat - format of input file ("mwa" or "text")
      varspec  - if informat is text, then this sets whether to read in a 
                 variance spectrum.  Default: varspec = True
   """

   if verbose:
      print ""
      print "Reading spectrum from %s" % filename

   """ Set default for variance spectrum """
   var = None

   """ Read in the input spectrum """
   if informat=="mwa":
      hdulist = pf.open(filename)
      flux = hdulist[1].data.copy()
      var  = hdulist[3].data.copy()
      varspec = True
      wavelength = np.arange(flux.size)
      hdr1 = hdulist[1].header
      wavelength = hdr1['crval1'] + wavelength*hdr1['cd1_1']
      del hdulist
   else:
      spec = np.loadtxt(filename)
      wavelength = spec[:,0]
      flux       = spec[:,1]
      if varspec:
         if spec.shape[1] < 3:
            print ''
            print '**Warning: varspec=True but only 2 columns in input file'
            print ''
         else:
            var = spec[:,2]
      del spec

   if verbose:
      print " Spectrum Start: %8.2f" % wavelength[0]
      print " Spectrum End:   %8.2f" % wavelength[-1]
      print " Dispersion (1st pixel): %6.2f" % (wavelength[1]-wavelength[0])
      print " Dispersion (average):   %6.2f" % \
          ((wavelength[-1]-wavelength[0])/(wavelength.size-1))
      print ""

   return wavelength, flux, var

#-----------------------------------------------------------------------

def save_spectrum(filename,x,flux,var=None):
   """
   Saves a spectrum as a text file
   """

   if var is not None:
      outdata = np.zeros((x.shape[0],3))
      outdata[:,2] = var
      fmtstring = '%7.2f %9.3f %10.4f'
   else:
      outdata = np.zeros((x.shape[0],2))
      fmtstring = '%7.2f %9.3f'
   outdata[:,0] = x
   outdata[:,1] = flux
   print ""
   print "Saving spectrum to file %s" % filename
   np.savetxt(filename,outdata,fmt=fmtstring)
   del outdata

#-----------------------------------------------------------------------

def plot_spectrum_array(x, flux, var=None, xlabel="Wavelength (Angstroms)",
                        ylabel="Relative Flux", title='Extracted Spectrum', 
                        docolor=True, speccolor='b', rmscolor='r', 
                        rmsoffset=0, rmsls=None, fontsize=12):

   """
   Given two input arrays, plot a spectrum.
   """

   """ Override color assignments if docolor is False"""
   if not docolor:
      speccolor = 'k'
      rmscolor  = 'k'
   plt.axhline(color='k')
   plt.plot(x,flux,speccolor,linestyle='steps',label='Flux')
   plt.tick_params(labelsize=fontsize)
   plt.xlabel(xlabel,fontsize=fontsize)
   if var is not None:
      rms = np.sqrt(var)+rmsoffset
      if rmsls is None:
         if docolor:
            rlinestyle = 'steps'
         else:
            rlinestyle = 'steps:'
      else:
         rlinestyle = 'steps%s' % rmsls
      if docolor:
         plt.plot(x,rms,rmscolor,linestyle=rlinestyle,label='RMS')
      else:
         plt.plot(x,rms,rmscolor,linestyle=rlinestyle,label='RMS',lw=2)
   plt.ylabel(ylabel,fontsize=fontsize)
   if(title):
      plt.title(title)
   if(x[0] > x[-1]):
      plt.xlim([x[-1],x[0]])
   else:
      plt.xlim([x[0],x[-1]])

#-----------------------------------------------------------------------

def plot_spectrum(filename, varspec=True, informat="text", 
                  xlabel="Wavelength (Angstroms)", ylabel="Relative Flux",
                  title='default', fontsize=12, docolor=True, speccolor='b', 
                  rmscolor='r', rmsoffset=0, rmsls=None,
                  add_atm_trans=False, atmscale=1.05, atmfwhm=15.,
                  atmoffset=0., atmls='-', verbose=True):
   """
   Given an input file with spectroscopy information, plot a spectrum.  
   The input file can have one of two formats: text (the default) or mwa
   For the text format, the input file is an ascii text file, with 2 or 3 columns
      1. Wavelength (or pixel position in the dispersion direction)
      2. Flux (counts, etc.)
      3. Variance spectrum (OPTIONAL)
   The mwa format is created by the make_spec function in the spectools
    library (probably found in the mostools package).
    The format is a multi-extension FITS file:
      Extension 1: the science spectrum
      Extension 2: a smoothed version of the science spectrum
      Extension 3: the variance spectrum.
    NOTE: For this format, the wavelength information is contained in the
      FITS header cards, and is not stored as a separate array.
   """

   """ Read in the spectrum, using the appropriate input format """
   wavelength,flux,var = read_spectrum(filename, informat, varspec, verbose)

   """ Plot the spectrum """
   if title == 'default':
      title = 'Spectrum for %s' % filename
   if varspec:
      plot_spectrum_array(wavelength,flux,var=var,xlabel=xlabel,ylabel=ylabel,
                          title=title,docolor=docolor,speccolor=speccolor,
                          rmscolor=rmscolor,rmsoffset=rmsoffset,
                          rmsls=rmsls,fontsize=fontsize)
   else:
      plot_spectrum_array(wavelength,flux,xlabel=xlabel,ylabel=ylabel,
                          title=title,docolor=docolor,speccolor=speccolor,
                          rmscolor=rmscolor,rmsoffset=rmsoffset,
                          fontsize=fontsize)

   """ Plot the atmospheric transmission if requested """
   if add_atm_trans:
      plot_atm_trans(wavelength, atmfwhm, flux, scale=atmscale, 
                     offset=atmoffset, linestyle=atmls)

   del wavelength
   del flux
   if(varspec):
      del var

#-----------------------------------------------------------------------

def plot_blue_and_red(bluefile, redfile, outfile=None, smooth_width=7,
                      bscale=10., xlim=[3000.,9500.], title='default', 
                      z=None, mark_em=False, mark_abs=True):
   """
   Creates a single plot given files that contain data from the blue and red
   sides of a spectrograph.

   Required inputs:
      bluefile - file containing the blue-side spectral data, which (for now)
                 is expected to be a text file containing two or three columns:
                    wavelength flux [variance]
      redfile  - file containing the red-side spectral data.  Same format
                 as for bluefile
   Optional inputs:
      outfile      - Name of output file containing the plot.  The default value
                     (None) means no output file.
      smooth_width - The input spectra will be smoothed by a boxcar of 
                     width = smooth_width (default value = 7).  If no
                     smoothing is desired, then set smooth_width=None
      bscale       - Factor to multiply the blue-side flux by in order to
                     get it to match (roughly) the red-side flux where
                     the two spectra meet/overlap. Default value = 10.
      xlim         - x-axis range to be shown on the plot. Default=[3000.,9500.]
      z            - Redshift of the object.  Default (None) means that no
                     spectral lines are marked.  If z is not None, then
                     mark spectral lines, as set by the mark_em and mark_abs
                     parameters
      mark_em      - Mark emission lines on the spectrum. Default=False
      mark_abs     - Mark absorption lines on the spectrum. Default=True
   """

   """ Read in the spectra, including the scaling for the blue side """
   wb,fb,vb = read_spectrum(bluefile)
   wr,fr,vr = read_spectrum(redfile)
   fb *= bscale
   if vb is not None:
      vb *= bscale

   """ Smooth the data """
   if smooth_width is None:
      smb  = fb
      smvb = vb
      smr  = fr
      smvr = vr
   else:
      from scipy.ndimage.filters import uniform_filter as unif_filt
      if vb is None:
         tmpvb = np.ones(wb.size)
      else:
         tmpvb = vb
      if vr is None:
         tmpvr = np.ones(wr.size)
      else:
         tmpvr = vr
      wtb = 1./tmpvb
      wtr = 1./tmpvr
      tmpb = wtb * fb
      tmpr = wtr * fr
      smb = unif_filt(tmpb,smooth_width)
      smb /= unif_filt(wtb,smooth_width)
      smr = unif_filt(tmpr,smooth_width)
      smr /= unif_filt(wtr,smooth_width)
      if vb is None:
         smvb = None
      else:
         smvb = 1.0 / (smooth_width * unif_filt(wtb,smooth_width))
      if vr is None:
         smvr = None
      else:
         smvr = 1.0 / (smooth_width * unif_filt(wtr,smooth_width))
      

   """ Plot the spectra """
   plot_spectrum_array(wb,smb,smvb,rmscolor='k')
   plot_spectrum_array(wr,smr,smvr,speccolor='r',rmscolor='k')
   plt.xlim(xlim[0],xlim[1])

   """ Mark the lines if the redshift is given """
   if z is not None:
      is_z_shown = False
      w = np.concatenate((wb,wr))
      f = np.concatenate((smb,smr))
      if mark_em:
         mark_spec_emission(z,w,f)
         is_z_shown = True
      if mark_abs:
         if is_z_shown:
            mark_spec_absorption(z,w,f,showz=False)
         else:
            mark_spec_absorption(z,w,f)


#-----------------------------------------------------------------------

def subtract_sky(data, outfile, outskyspec, dispaxis='x', doplot=True):
   """
   Given the input 2D spectrum, creates a median sky and then subtracts
   it from the input data.  Two outputs are saved: the 2D sky-subtracted
   data and a 1D sky spectrum.

   Inputs:
      data       - array containing the 2D spectrum
      outfile    - name for output fits file containing sky-subtracted spectrum
      outskyspec - name for output 1D sky spectrum
   """
   """ Set the dispersion axis direction """
   if dispaxis == "y":
      spaceaxis = 1
   else:
      spaceaxis = 0
   
   """ Take the median along the spatial direction to estimate the sky """
   if data.ndim < 2:
      print ""
      print "ERROR: subtract_sky needs a 2 dimensional data set"
      return
   else:
      sky1d = np.median(data,axis=spaceaxis)
   sky = np.zeros(data.shape)
   for i in range(data.shape[spaceaxis]):
      if spaceaxis == 1:
         sky[:,i] = sky1d
      else:
         sky[i,:] = sky1d

   """ Plot the sky if desired """
   x = np.arange(sky1d.size)
   if doplot:
      if spaceaxis == 1:
         xlab = 'Row'
      else:
         xlab = 'Column'
      plot_spectrum_array(x,sky1d,xlabel=xlab,ylabel='Median Counts',
                          title='Median sky spectrum')

   """ Subtract the sky  """
   skysub = data - sky

   """ Save the sky-subtracted spectrum and the median sky """
   pf.PrimaryHDU(skysub).writeto(outfile)
   print ' Wrote sky-subtracted data to %s' % outfile
   save_spectrum(outskyspec,x,sky1d)
   print ' Wrote median sky spectrum to %s' % outskyspec
   print ''

   """ Clean up """
   del sky,sky1d,skysub


#-----------------------------------------------------------------------

#def make_gauss(x,mu,sigma,amp,bkgd):
def make_gauss(x,p):
   """
   Creates a model comprised of one or more Gaussian profiles plus a
    (for now) constant background.  
   NOTE: the only oddity is that, if there are more than one Gaussian in
    the profile, then the "mu" term for the subsequent Gaussians 
    (i.e., p[4], p[7], ..) are actually _offsets_ between the mean of the
    subsequent Gaussian and the mean of the first.  For example,
    mu_2 = p[0] + p[4]
    
   Inputs:
    x  - The independent variable that is used to generate the model, ymod(x)
    p  - The parameter values.  The length of this vector will be 1+3*n, where
          n>=1, for one constant background parameter (p[0]) plus one or more 
          Gaussian parameters, which come in sets of three.  Thus, p can
          be decomposed as follows:
          p[0] - background: required
          p[1] - mu_1: required
          p[2] - sigma_1: required
          p[3] - amplitude_1: required
          p[4] - offset between mu_2 and mu_1: optional
          p[5] - sigma_2: optional
          p[6] - amplitude_2: optional
          ... etc. for as many Gaussians are used to construct the profile
   """

   """ Calculate the number of Gaussians in the model """
   ngauss = int((p.size-1)/3)
   if p.size - (ngauss*3+1) != 0:
      print ''
      print 'ERROR: Gaussian model contains the incorrect number of parameters'
      print ''
      return np.nan

   """ Calculate y_mod using current parameter values """
   ymod = np.zeros(x.size) + p[0]
   for i in range(ngauss):
      ind = i*3+1
      if i == 0:
         mu = p[ind]
      else:
         mu = p[1] + p[ind]
         
      ymod += p[ind+2] * np.exp(-0.5 * ((x - mu)/p[ind+1])**2)

   return ymod

#-----------------------------------------------------------------------

def fit_gauss(p,x,y,p_init,fitind):
   """
   Compares the data to the model.  The model consists of at least one gaussian 
    plus a constant background and is created by a call to make_gauss.
    Thus the comparison is between ymod(x) and y, where the latter is the
    measured quantity.

   NOTE: the only oddity in the model is that, if there are more than one 
    Gaussian in the profile, then the "mu" term for the subsequent Gaussians 
    (i.e., p[4], p[7], ..) are actually _offsets_ between the mean of the
    subsequent Gaussian and the mean of the first.  For example,
    mu_2 = p[0] + p[4]
    
   Inputs:
    p  - The parameter values.  The length of this vector will be 1+3*n, where
          n>=1, for one constant background parameter (p[0]) plus one or more 
          Gaussian parameters, which come in sets of three.  Thus, p can
          be decomposed as follows:
          p[0] - background: required
          p[1] - mu_1: required
          p[2] - sigma_1: required
          p[3] - amplitude_1: required
          p[4] - offset between mu_2 and mu_1: optional
          p[5] - sigma_2: optional
          p[6] - amplitude_2: optional
          ... etc. for as many Gaussians are used to construct the profile
    x  - The independent variable that is used to generate the model, ymod(x)
    y  - The measured data, to be compared 
   """

   """
   Create the full list of model parameters by combining the fitted parameters
   and the fixed parameters
   """
   pfull = p_init.copy()
   pfull[fitind] = p

   """
   Compute the difference between model and real values
   """
   ymod = make_gauss(x,pfull)
   diff = y - ymod

   return diff

#-----------------------------------------------------------------------

def fit_gpb_fixmusig(p,x,y,mu,sigma):
   """
   Compares the data to the model.  The model is a gaussian plus a 
    constant background.  In the fit, mu and sigma are held fixed.
   The parameter values are:
    p[0] = background level
    p[1] = amplitude
   """

   """ 
   Put parameters into the form that is expected by make_gauss
   """
   ptmp = np.array([p[0],mu,sigma,p[1]])

   """
   Compute the difference between model and real values
   """

   ymod = make_gauss(x,ptmp)
   diff = y - ymod

   return diff

#-----------------------------------------------------------------------

def fit_gpb_fixmu(p,x,y,mu):
   """
   Compares the data to the model.  The model is a gaussian plus a 
    constant background.  In the fit, mu is held fixed.
   The parameter values are:
    p[0] = background level
    p[1] = amplitude
    p[2] = sigma
   """

   """ Unpack p """
   bkgd = p[0]
   amp  = p[1]
   sigma = p[2]
   if len(p) > 3:
      nps = (len(p)-1)/2
      for inpsf in range(1,nps):
         amp,sigma = np.append(amp,p[2*inpsf+1]),np.append(sigma,p[2*inpsf+2])

   """
   Compute the difference between model and real values
   """
   if np.shape(amp) != (): mu = np.ones(len(bkgd))*mu
   ymod = make_gauss(x,mu,sigma,amp,bkgd)
   diff = y - ymod

   return diff

#-----------------------------------------------------------------------

def plot_spatial_profile(infile, dispaxis="x"):
   """
   Given an input fits file with (by assumption) a 2d spectrum, this
   function will compress the spectrum along the dispersion direction
   to get an average spatial profile.

   Inputs:
      infile   - input file containing the 2d spectrum
      dispaxis - axis corresponding to the dispersion direction (i.e.,
                 the wavelength axis)
   """

   # Read the data
   data = load_raw_spectrum(infile)

   # Set the dispersion axis direction
   if dispaxis == "y":
      specaxis = 0
      spatlabel = "x"
   else:
      specaxis = 1
      spatlabel = "y"
   #print "specaxis = %d" % specaxis

   """ Compress the data along the dispersion axis and find the max value """
   if data.ndim < 2:
      print ""
      print "ERROR: plot_spatial_profile needs a 2 dimensional data set"
      del data
      return
   else:
      cdat = np.median(data,axis=specaxis)
   x = np.arange(1,cdat.shape[0]+1)

   """ Plot the spatial profile """
   plt.plot(x,cdat)
   plt.xlabel("Pixel in the %s direction" % spatlabel)
   plt.ylabel("Median counts")
   plt.title("Spatial profile for %s" % infile)

#-----------------------------------------------------------------------

def find_peak(data,dispaxis="x",mu0=None,sig0=None,fixmu=False,fixsig=False,
              showplot=True,do_subplot=False,verbose=True,apmin=-4.,apmax=4.):
   """
    Compresses a 2d spectrum along the dispersion axis so that
     the trace of the spectrum can be automatically located by fitting
     a gaussian + background to the spatial direction.  The function
     returns the parameters of the best-fit gaussian.
    The default dispersion axis is along the x direction.  To change this
     set the optional parameter dispaxis to "y"
   """

   # Set the dispersion axis direction
   if dispaxis == "y":
      specaxis = 0
   else:
      specaxis = 1
   #print "specaxis = %d" % specaxis

   """ Compress the data along the dispersion axis and find the max value """
   if data.ndim < 2:
      cdat = data.copy()
   else:
      cdat = np.median(data,axis=specaxis)
      cdat.shape
   x = np.arange(1,cdat.shape[0]+1)

   # Set initial guesses

   if fixmu:
      if mu0 is None:
         print ""
         print "ERROR: find_peak.  mu is fixed, but no value for mu0 given"
         return
      fixmunote = "**"
   else:
      if mu0 is None:
         i = cdat.argsort()
         mu0    = 1.0 * i[i.shape[0]-1]
      fixmunote = " "
   if fixsig:
      if sig0 is None:
         print ""
         print "ERROR: find_peak.  sigma is fixed, but no value for sig0 given"
         return
      fixsignote = "**"
   else:
      if sig0 is None:
         sig0 = 3.0
      fixsignote = " "
   amp0  = cdat.max()
   bkgd0 = np.median(data,axis=None)
   if(verbose):
      print ""
      print "Initial guesses for Gaussian plus background fit"
      print "------------------------------------------------"
      print " mu         = %7.2f%s"   % (mu0,fixmunote)
      print " sigma      =   %5.2f%s" % (sig0,fixsignote)
      print " amplitude  = %f"        % amp0
      print " background = %f"        % bkgd0
      print "Parameters marked with a ** are held fixed during the fit"
      print ""

   # Fit a Gaussian plus a background to the compressed spectrum
   mf=100000
   if fixmu and fixsig:
      p = [bkgd0,amp0]
      pt,ier = optimize.leastsq(fit_gpb_fixmusig,p,(x,cdat,mu0,sig0),maxfev=mf)
      p_out = [pt[0],mu0,sig0,pt[1]]
   #p_out,ier = optimize.leastsq(fit_gpb_fixmu,p,(x,cdat,mu0),maxfev=mf)
   #p_out,ier = optimize.leastsq(fit_gpb_fixsig,p,(x,cdat,sig0),maxfev=mf)
   else:
      p = [bkgd0,mu0,sig0,amp0]
      p_out,ier = optimize.leastsq(fit_gauss,p,(x,cdat),maxfev=mf)


   # Give results
   if(verbose):
      print "Fitted values for Gaussian plus background fit"
      print "----------------------------------------------"
      print " mu         = %7.2f%s"   % (p_out[1],fixmunote)
      print " sigma      =   %5.2f%s" % (p_out[2],fixsignote)
      print " amplitude  = %f"        % p_out[3]
      print " background = %f"        % p_out[0]
      print "Parameters marked with a ** are held fixed during the fit"
      print ""

   # Plot the compressed spectrum
   if(showplot):
      if(do_subplot):
         plt.subplot(221)
      else:
         plt.figure(1)
         plt.clf()
      plt.plot(x,cdat,linestyle='steps')
      xmod = np.arange(1,cdat.shape[0]+1,0.1)
      ymod = make_gauss(xmod,p_out[1],p_out[2],p_out[3],p_out[0])
      plt.plot(xmod,ymod)
      plt.axvline(p_out[1]+apmin,color='k')
      plt.axvline(p_out[1]+apmax,color='k')
      plt.xlabel('Pixel number in the spatial direction')
      plt.title('Compressed Spatial Plot')

   return p_out

#-----------------------------------------------------------------------

def extract_wtsum_col(spatialdat,mu,apmin,apmax,weight='gauss',sig=1.0,
                      gain=1.0,rdnoise=0.0,sky=None):
   """
   Extracts the spectrum from one row/column in the wavelength direction
   via a weighted sum.  The choices for the weighting are:
       'gauss'   - a Gaussian, where mu and sigma of the Gaussian are fixed.
       'uniform' - uniform weighting across the aperture, which is centered
                   at mu

   Inputs:
     spatialdat - a one-dimensional array, corresponding to a cut in
                  the spatial direction from the 2-d spectrum
     mu         - the fixed centroid of the trace
     apmin      - the lower bound, with respect to mu, of the aperture to be
                  extracted
     apmax      - the upper bound, with respect to mu, of the aperture to be
                  extracted
     weight     - the weighting scheme to be used.  Valid choices are:
                  'gauss'   (the default value)
                  'uniform' 
     sig        - the fixed sigma of the Gaussian fit to the trace
     gain       - CCD gain - used to compute the variance spectrum
     rdnoise    - CCD readnoise  - used to compute the variance spectrum
     sky        - sky value for this wavelength (default=None).  Used only
                  if the spectrum passed to this function has already been
                  background-subtracted.
   """

   """ Define aperture and background regions """
   y = np.arange(spatialdat.shape[0])
   ydiff = y - mu
   apmask = (ydiff>apmin-1) & (ydiff<apmax)
   bkgdmask = np.logical_not(apmask)
   #print bkgdmask.sum(), apmask.sum()

   """ Estimate the background """
   bkgd = np.median(spatialdat[bkgdmask],axis=None)
   #print "Background level is %7.2f" % bkgd

   """ Make the weight array """
   if(weight == 'uniform'):
      gweight = np.zeros(y.size)
      gweight[apmask] = 1.0
   else:
      gweight = make_gauss(y,mu,sig,1.0,0.0)

   """ Do the weighted sum """
   wtsum = ((spatialdat - bkgd)*gweight)[apmask].sum() / gweight[apmask].sum()

   """ Calculate the variance """
   if (sky == None):
      varspec = (gain * spatialdat + rdnoise**2)/gain**2
      var = (varspec * gweight)[apmask].sum() / gweight[apmask].sum()
   else:
      varspec = (gain * (spatialdat + sky) + rdnoise**2)/gain**2
      var = (varspec * gweight)[apmask].sum() / gweight[apmask].sum()

   return wtsum, var, bkgd

#-----------------------------------------------------------------------

def find_trace(data,dispaxis="x",apmin=-4.,apmax=4.,doplot=True,
               do_subplot=False):
   """
   First step in the reduction process.  
   Runs find_peak on the full 2d spectrum in order to set the initial 
    guess for trace_spectrum.
   """
   p = find_peak(data,dispaxis=dispaxis,apmin=apmin,apmax=apmax,
                 showplot=doplot,do_subplot=do_subplot)
   mu0  = p[1]
   sig0 = p[2]
   return mu0,sig0

#-----------------------------------------------------------------------

def fit_poly_to_trace(x, data, fitorder, data0, x_max, fitrange=None,
                      doplot=True, markformat='bo', ylabel='Centroid of Trace',
                      title='Location of the Peak'):

   # Do a sigma clipping to reject clear outliers
   if fitrange is None:
      tmpfitdat = data
   else:
      fitmask = np.logical_and(x>=fitrange[0],x<fitrange[1])
      tmpfitdat = data[fitmask]
   dmu,dsig = ccd.sigma_clip(tmpfitdat,3.0)
   goodmask = np.absolute(data - dmu)<3.0*dsig
   badmask  = np.absolute(data - dmu)>=3.0*dsig
   dgood    = data[goodmask]
   dbad     = data[badmask]
   xgood    = x[goodmask]
   xbad     = x[badmask]

   # Fit a polynomial to the trace

   if fitrange is None:
      xpoly = xgood
      dpoly = dgood
   else:
      fitmask = np.logical_and(xgood>=fitrange[0],xgood<fitrange[1])
      #print fitmask
      xpoly  = xgood[fitmask]
      dpoly  = dgood[fitmask]

   if fitorder == -1:
      polyorder = 0
   else:
      polyorder = fitorder
   dpoly = np.polyfit(xpoly,dpoly,polyorder)

   if fitorder == -1:
      dpoly[0] = data0

   # Plot the results

   ymin = dmu - 4.5*dsig
   ymax = dmu + 4.5*dsig
   if doplot:
      plt.plot(x,data,markformat)
      #plt.plot(xstep,mu,marker='o',mec='b',mfc='w',markersize=8,linestyle='')
      plt.xlabel("Pixel number in the dispersion direction")
      plt.ylabel(ylabel)
      plt.title(title)

      # Show the value from the compressed spatial profile
      plt.axhline(data0,color='k',linestyle='--')

      # Mark the bad points that were not included in the fit
      plt.plot(xbad,dbad,"rx",markersize=10,markeredgewidth=2)

      # Show the fitted function
      fitx = np.arange(0,x_max,0.1)
      fity = 0.0 * fitx
      for i in range(dpoly.size):
         fity += dpoly[i] * fitx**(dpoly.size - 1 - i)
      plt.plot(fitx,fity,"r")

      # Show the range of points included in the fit, if fitrange was set
      if fitrange is not None:
         plt.axvline(fitrange[0],color='k',linestyle=':')
         plt.axvline(fitrange[1],color='k',linestyle=':')
         xtmp = 0.5 * (fitrange[1] + fitrange[0])
         xerr = xtmp - fitrange[0]
         ytmp = fity.min() - 0.2 * fity.min()
         plt.errorbar(xtmp,ytmp,xerr=xerr,ecolor="g",capsize=10)
      plt.xlim(0,x_max)
      plt.ylim(ymin, ymax)

   # Return the parameters produced by the fit
   print dpoly
   return dpoly


#-----------------------------------------------------------------------

def trace_spectrum(data,mu0,sig0,dispaxis="x",stepsize=25,muorder=3,
   sigorder=4,fitrange=None,doplot=True,do_subplot=False):
   """
   Second step in the reduction process.
   Fits a gaussian plus background to portions of the spectrum separated
   by stepsize pixels (default is 25).
   """

   # Set the dispersion axis direction
   if dispaxis == "y":
      specaxis  = 0
      spaceaxis = 1
   else:
      specaxis  = 1
      spaceaxis = 0
   xlength   = data.shape[specaxis]

   # Define the slices through the 2D spectrum that will be used to find
   #  the centroid and width of the object spectrum as it is traced down 
   #  the chip
   xstep = np.arange(0,xlength-stepsize,stepsize)

   # Set up containers for mu and sigma along the trace
   mu = 0.0 * xstep
   sigma = 0.0 * xstep
   nsteps = np.arange(xstep.shape[0])

   # Step through the data
   print ""
   print "Running fit_trace"
   print "--------------------------------------------------------------- "
   print "Finding the location and width of the peak at %d segments of " % \
     nsteps.shape[0]
   print"   the 2D spectrum..."
   for i in nsteps:
      if(specaxis == 0):
         tmpdata = data[xstep[i]:xstep[i]+stepsize,:]
      else:
         tmpdata = data[:,xstep[i]:xstep[i]+stepsize]
      ptmp = find_peak(tmpdata,dispaxis=dispaxis,showplot=False,verbose=False)
      mu[i]    = ptmp[1]
      sigma[i] = ptmp[2]
   print "   Done"

   # Fit a polynomial to the trace
   if doplot:
      if(do_subplot):
         plt.subplot(222)
      else:
         plt.figure(2)
         plt.clf()
   print "Fitting a polynomial of order %d to the location of the trace" \
       % muorder
   mupoly = fit_poly_to_trace(xstep,mu,muorder,mu0,xlength,fitrange,
                              doplot=doplot)

   # Fit a polynomial to the width of the trace
   if(do_subplot):
      plt.subplot(223)
   else:
      plt.figure(3)
      plt.clf()
   print "Fitting a polynomial of order %d to the width of the trace" % sigorder
   sigpoly = fit_poly_to_trace(xstep,sigma,sigorder,sig0,xlength,fitrange,
                               markformat='go',title='Width of Peak',
                               ylabel='Width of trace (Gaussian sigma)',
                               doplot=doplot)

   # Return the fitted parameters
   return mupoly,sigpoly

#--------------------------------------------------------------------------

def find_and_trace(data, dispaxis="x", apmin=-4., apmax=4., stepsize=25,
                   muorder=3, sigorder=4, fitrange=None, doplot=True,
                   do_subplot=False):

   """
   Combines the first two steps in the spectroscopy reduction process.
   The find_and_trace function will:
      1. Locate roughly where the target object is in the spatial direction 
         (usually the y axis is the spatial direction) by taking the
         median in the spectral direction so the peak in the spatial
         direction stands out.  This step provides the initial guesses
         for the location (mu0) and width (sig0) of the peak that are
         then used in the second step.
      2. Once the rough location of the peak has been found, determines how
         its location and width change in the spectral direction.
         That is, this will find where the peak falls in each column.
         It returns the position (pos) and width (width) of the peak as
         a function of x location

   This function accomplishes these tasks by calling first the find_trace
   function and the the trace_spectrum function.
   """

   mu0,sig0 = find_trace(data,dispaxis,apmin,apmax,doplot,do_subplot)

   pos,width = trace_spectrum(data,mu0,sig0,dispaxis,stepsize,muorder,
                              sigorder,fitrange,doplot,do_subplot)

   return pos,width

#-----------------------------------------------------------------------

def extract_spectrum(data,mupoly,sigpoly,dispaxis="x",apmin=-4.,apmax=4.,
                     weight='gauss', sky=None, gain=1.0, rdnoise=0.0, 
                     doplot=True, do_subplot=False, outfile=None):
   """
   Third step in reduction process.
   This function extracts a 1D spectrum from the input 2D spectrum (data)
   It uses the information about the trace that has been generated by
   the trace_spectrum function.  In particular, it takes the two
   polynomials generated by trace_spectrum as the inputs mupoly and sigpoly.
   """

   # Set the dispersion axis direction
   if dispaxis == "y":
      specaxis  = 0
      spaceaxis = 1
   else:
      specaxis  = 1
      spaceaxis = 0

   # Set the wavelength axis
   pix = np.arange(data.shape[specaxis])

   # Set the fixed mu and sigma for the Gaussian fit at each point in the
   #  spectrum, using the input polynomials

   fitx = np.arange(data.shape[specaxis]).astype(np.float32)
   mu = 0.0 * fitx
   for i in range(mupoly.size):
      mu += mupoly[i] * fitx**(mupoly.size - 1 - i)
   sig = 0.0 * fitx
   for i in range(sigpoly.size):
      sig += sigpoly[i] * fitx**(sigpoly.size - 1 - i)

   # Set up the containers for the amplitude and variance of the spectrum
   #  along the trace
   amp = 0.0 * pix
   var = 0.0 * pix

   # Step through the data
   print ""
   print "Extracting the spectrum..."
   for i in pix:
      #print pix[i], mu[i], sig[i]
      if specaxis == 0:
         tmpdata = data[i,:]
      else:
         tmpdata = data[:,i]
      if (sky == None):
         skyval = None
      else:
         skyval = sky[i]
      amp[i],var[i] = extract_wtsum_col(tmpdata,mu[i],apmin,apmax,sig=sig[i],
                                        gain=gain,rdnoise=rdnoise,sky=skyval,
                                        weight=weight)
   print "   Done"

   # Plot the extracted spectrum
   if doplot:
      print ""
      print "Plotting the spectrum"
      if(do_subplot):
         plt.subplot(224)
      else:
         plt.figure(4)
         plt.clf()
      plot_spectrum_array(pix,amp,title='Extracted spectrum (pixels)',
                          xlabel='Pixel number in the dispersion direction')

   # Save the extracted spectrum
   if outfile is not None:
      save_spectrum(outfile,pix,amp,var)

   # Return the extracted spectrum
   return amp,var

#-----------------------------------------------------------------------

def resample_spec(w, spec, owave=None):
   """
   Given an input spectrum, represented by wavelength values (w) and fluxes
    (spec), resample onto a linearized wavelength grid.  The grid can either
    be defined by the input wavelength range itself (the default) or by
    a wavelength vector that is passed to the function.
   """

   if owave is None:
      w0 = w[0]
      w1 = w[-1]
      owave = np.linspace(w0,w1,w.size)

   specmod = interpolate.splrep(w,spec)
   outspec = interpolate.splev(owave,specmod)

   return owave,outspec

#-----------------------------------------------------------------------

def combine_spectra(txt_files,outfile):
   """
   Given the input spectra, stored as text files (e.g., "vega*txt"), reads in 
   the files and does an inverse-variance weighted combination of the counts 
   columns.
   """

   file_list = glob.glob(txt_files)

   """ Setup """
   tmpdat = np.loadtxt(file_list[0])
   wavelength = tmpdat[:,0].copy()
   wtflux = wavelength * 0.0
   wtsum  = wavelength * 0.0

   """ Create the weighted sum """
   print ""
   for f in file_list:
      print "Reading data from file %s" % f 
      wi,fi,vi = np.loadtxt(f,unpack=True,usecols=(0,1,2))
      wt = 1.0 / vi
      wt[vi==0] = 0.
      wtflux += wt * fi
      wtsum += wt

   """ 
   Normalize the flux, and calculate the variance of the coadded spectrum.
   Note that the equation below for the variance only works for the case
    of inverse variance weighting.
   """
   wtflux[wtsum==0] = 0
   wtsum[wtsum==0] = 1
   outspec = wtflux / wtsum
   outvar  = 1.0 / wtsum

   """ Plot the combined spectrum """
   plot_spectrum_array(wavelength,outspec,outvar,xlabel="Pixels",
                       title="Combined spectrum")

   print "Saving the combined spectrum"
   save_spectrum(outfile,wavelength,outspec,outvar)

#-----------------------------------------------------------------------

def plot_sky(infile):
   """
   Given an input 2-dimensional fits file for which the sky has NOT been 
   subtracted, this function will take the median along the spatial axis
   to produce a sky spectrum. 

   *** NB: Right now this ASSUMES that the dispersion is in the x direction ***
   """

   data = pf.getdata(infile)
   sky = np.median(data,axis=0)
   pix = np.arange(sky.size)
   plot_spectrum_array(pix,sky,xlabel='Pixels',title='Sky Spectrum')

#-----------------------------------------------------------------------

def make_sky_model(wavelength, smoothKernel=25., verbose=True):
   """
   Given an input wavelength vector, creates a smooth model of the
   night sky emission that matches the wavelength range and stepsize
   of the input vector.
   """
   import pickle

   # Get info from input wavelength vector
   wstart = wavelength.min()
   wend = wavelength.max()
   disp = wavelength[1] - wavelength[0]
   if verbose:
      print "Making model sky"
      print "--------------------------------------"
      print "Model starting wavelength: %f" % wstart
      print "Model ending wavelength:   %f" % wend
      print "Model dispersion:          %f" % disp

   # Use NIR sky model if starting wavelength is redder than 9000 Angstrom
   if wstart >= 9000.:
      # Read in skymodel, which is in a B-spline format
      skymodel_file = '/Users/cdf/git/KeckCDF/python/nirspec/nirspec_skymodel.dat'
      skymodel = np.load(skymodel_file)
   else:
      # Read in the sky model
      skymodel_file = '/Users/cdf/Code/python/LRISredux/data/uves_sky.model'
      f = open(skymodel_file)
      skymodel = pickle.load(f)
      f.close()

   # Resample and smooth the model spectrum
   wave = np.arange(wstart,wend,0.2)
   tmpskymod = interpolate.splev(wave,skymodel)
   tmpskymod = ndimage.gaussian_filter(tmpskymod,smoothKernel)

   # Create a B-spline representation of the smoothed curve for use in
   #  the wavecal optimization
   model = interpolate.splrep(wave,tmpskymod)

   # Finally use the initial guess for the dispersion and evaluate the
   #  model sky at those points, using the B-spline model
   skyflux = interpolate.splev(wavelength,model)

   """ Create a Spec1d instance containing the sky model """
   skymod = Spec1d(wav=wavelength,flux=skyflux)

   # Clean up and return
   del skymodel,tmpskymod,wave
   return skymod

#-----------------------------------------------------------------------

def apply_wavecal(infile, outfile, lambda0, dlambda, varspec=True):
   """
   Given an input file containing 2 columns (x and flux), and the y-intercept
   and slope of the x-lambda transformation, convert x to wavelength units
   and write the output as outfile.
   """

   """ Read the input file """
   x,flux,var = read_spectrum(infile,varspec=varspec)

   """ Convert x from pixels to wavelength units """
   wavelength = lambda0 + dlambda * x

   """ Plot and save the results """
   if varspec == True:
      plot_spectrum_array(wavelength,flux,var=var,
                          title="Wavelength-calibrated spectrum")
      save_spectrum(outfile,wavelength,flux,var)
   else:
      plot_spectrum_array(wavelength,flux,title="Wavelength-calibrated spectrum")
      save_spectrum(outfile,wavelength,flux)


#-----------------------------------------------------------------------

def check_wavecal(infile, informat='text', modsmoothkernel=25.):
   """
   MOVE TO SPEC1D CLASS!

   Plots the wavelength-calibrated sky information from the input file
   on top of a smoothed model of the night sky emission so that
   the quality of the wavelength calibration can be evaluated.

   The input file can either be a fits file (fileformat='fits') that has
   been produced by the niredux function or a text spectrum containing
   three columns (wavelength, flux, variance).  The variance spectrum in
   this case should contain clear sky line features if the exposure length
   was long enough.
   """

   """ Read in the observed sky spectrum """
   if informat=='fits':
      hdulist = pf.open(infile)
      varspec = hdulist[1].data.copy()
      skyobs  = np.sqrt(np.median(varspec,axis=0))
      skylab  = "RMS Spectrum"
      hdr = hdulist[0].header
      crval1  = hdr['crval1']
      crpix1  = hdr['crpix1']
      cd11    = hdr['cd1_1']
      waveobs = 1.0* np.arange(varspec.shape[1])
      waveobs *= cd11
      waveobs += crval1
   elif informat=='fitsold':
      hdulist = pf.open(infile)
      waveobs = hdulist[1].data.copy()
      skyobs  = hdulist[2].data.copy()
      skylab  = "Observed Sky"
   else:
      try:
         waveobs,varspec = np.loadtxt(infile,unpack=True,usecols=(0,2))
         skyobs = np.sqrt(varspec)
      except:
         print ""
         print "Cannot get variance spectrum from input text file %s" % infile
         print ""
         return
      skylab = "RMS Spectrum"

   """ Create the sky spectrum, with the appropriate smoothing """
   print ""
   skymod = make_sky_model(waveobs,modsmoothkernel)

   """ 
   Scale the sky spectrum to roughly be 75% of the amplitude of the observed
   spectrum
   """

   deltamod = skymod.flux.max() - skymod.flux.min()
   deltaobs = skyobs.max() - skyobs.min()
   skymod.flux *= 0.75*deltaobs/deltamod
   skymod.flux += skyobs.mean() - skymod.flux.mean()

   """ Make the plot """
   wrange = waveobs.max() - waveobs.min()
   xmin = waveobs.min() - 0.05*wrange
   xmax = waveobs.max() + 0.05*wrange
   plt.plot(waveobs,skyobs,'k',ls='steps', label=skylab)
   skymod.plot(color='r',label='Model sky')
   plt.legend()
   plt.xlim(xmin,xmax)

   del waveobs
   del skyobs
   del skymod

#-----------------------------------------------------------------------

def planck_spec(wavelength, T=1.0e4, waveunit='Angstrom'):
   """
   Given a wavelength vector and an input temperture, generates a thermal
    spectrum over the input wavelength range.
   The input spectrum is B_lambda(T) and NOT B_nu(T)
   """

   # Define the constants in the Planck function in SI units
   c = 3.0e8
   h = 6.626e-34
   k = 1.38e-23

   # Convert the wavelength into meters (default assumption is that the
   #  input wavelength is in Angstroms
   wtmp = wavelength.copy()
   if waveunit[0:6].lower()=='micron':
      print "Converting wavelength from microns to meters"
      wtmp *= 1.0e-6
   elif waveunit[0:5].lower()=='meter':
      wtmp *= 1.0
   else:
      print "Converting wavelength from Angstroms to meters"
      wtmp *= 1.0e-10

   # Generate the Planck function, and then scale it so that its mean matches
   #  the mean of the observed spectrum, just for ease in plotting.
   from math import e
   denom = e**(h * c /(wtmp * k * T)) - 1.0
   B_lam = 2.0 * h * c**2 / (wtmp**5 * denom)

   return B_lam

#-----------------------------------------------------------------------

def response_ir(infile,outfile,order=6,fitrange=None,filtwidth=9):
   """
   Calculates the response function for an IR spectral setup using observations
   of a standard star.  The assumption is that the standard is a hot
   star (A or B), and therefore its spectrum is just a power law in the NIR.
   The steps are:
     (1) Divide the model spectrum by the observed spectrum
     (2) Fit a polynomial to the result
     (3) Write out the result to the output file
   The response function in the output file can then be multiplied by other
   spectra to do a response correction (an approximation of flux calibration).

   Inputs:
      infile:    Input file containing observed spectrum of the hot star
                  This file should have 3 columns (wavelength, flux, variance)
      outfile:   Output file that will contain the response function
      order:     Order of polynomial fit (default = 6)
      fitrange:  A list of 2-element lists, where each of the smaller lists
                  contains the starting and ending values for a range of
                  good data to include in the fit.
                  E.g., fitrange=[[20150.,21500.],[22000.,25000.]]
                  The default (fitrange=None) uses the full wavelength range
                  in the fit.
      filtwidth: Width of box used in the maximum filtering step, which is
                  used to minimize the number of absorption lines in the
                  spectrum before fitting a low-order polynomial to the result
                  (default = 9)
   """

   # Read the input spectrum
   wave,fluxobs,var = read_spectrum(infile)
   rms = np.sqrt(var)

   # Generate the thermal spectrum and normalize it
   B_lam = planck_spec(wave)
   B_lam *= fluxobs.mean() / B_lam.mean()

   # The features in the observed spectrum that deviate from a thermal
   #  spectrum should only be absorption lines.  Therefore, run a maximum
   #  filter

   flux = ndimage.filters.maximum_filter(fluxobs,filtwidth)

   # Calculate the observed response function
   respobs = B_lam / flux

   # Show some plots
   plt.figure(1)
   plt.clf()
   plt.plot(wave,fluxobs)
   plt.plot(wave,B_lam)
   plt.plot(wave,flux)
   plt.plot(wave,rms)
   plt.figure(2)
   plt.clf()
   plt.plot(wave,respobs)

   # Define the spectral range to be included in the fit
   if fitrange is not None:
      mask = np.zeros(respobs.size,dtype=np.bool)
      fitr = np.atleast_2d(np.asarray(fitrange))
      for i in range(fitr.shape[0]):
         wmask = np.logical_and(wave>fitr[i,0],wave<fitr[i,1])
         mask[wmask] = True
      wavegood = wave[mask]
      respgood = respobs[mask]
   else:
      wavegood = wave
      respgood = respobs

   # Fit a polynomial to the observed response function
   fpoly = np.polyfit(wavegood,respgood,order)
   print ""
   print "Fit a polynomial of order %d to curve in Figure 2." % order
   print "Resulting coefficients:"
   print "-----------------------"
   print fpoly

   # Convert polynomial into a smooth response function
   p = np.poly1d(fpoly)
   resp = p(wave)

   # Add the smooth response to the plot and show corrected curve
   plt.plot(wave,resp,'r')
   fc = fluxobs * resp
   plt.figure(3)
   plt.clf()
   plt.plot(wave,fc)
   plt.plot(wave,B_lam)

   # Write smooth response to output file
   out = np.zeros((wave.size,2))
   out[:,0] = wave
   out[:,1] = resp
   np.savetxt(outfile,out,'%8.3f  %.18e')

#-----------------------------------------------------------------------

def response_correct(infile, respfile, outfile):
   """
   Applies a response correction, calculated previously by response_ir
   or another function, to the input file.  The output is placed in
   outfile.

   Inputs:
      infile:   Input spectrum
      respfile: Response correction spectrum
      outfile:  Output spectrum
   """

   # Read input files
   try:
      w,f,v = np.loadtxt(infile,unpack=True)
   except:
      print ""
      print "ERROR: response_correct.  Unable to read input spectrum %s" \
          % infile
      return
   try:
      wr,resp = np.loadtxt(respfile,unpack=True)
   except:
      print ""
      print "ERROR: response_correct.  Unable to read response spectrum %s" \
          % respfile
      return

   # Apply the response correction and save the spectrum
   f *= resp
   v *= resp**2
   save_spectrum(outfile,w,f,v)

#-----------------------------------------------------------------------

def normalize(infile,outfile,order=6,fitrange=None,filtwidth=11):
   """
   Normalizes a spectrum by fitting to the continuum and then dividing the
    input spectrum by the fit.

   Inputs:
      infile:    File containing the input spectrum
                  This file should have 3 columns (wavelength, flux, variance)
      outfile:   Output file that will contain the normalized spectrum
      order:     Order of polynomial fit (default = 6)
      fitrange:  A list of 2-element lists, where each of the smaller lists
                  contains the starting and ending values for a range of
                  good data to include in the fit.
                  E.g., fitrange=[[20150.,21500.],[22000.,25000.]]
                  The default (fitrange=None) uses the full wavelength range
                  in the fit.
      filtwidth: Width of box used in the boxcar smoothing step, which is
                  used to minimize the number of outlier points in the input
                  spectrum before fitting the polynomial to the spectrum
                  (default = 9)
   """

   # Read the input spectrum
   wave,fluxobs,var = read_spectrum(infile)
   rms = np.sqrt(var)

   # Try to minimize outliers due to both emission and absorption
   #  lines and to cosmetic features (cosmic rays, bad sky-line subtraction).
   #  Do this by doing a inverse-variance weighted boxcar smoothing.

   wt = 1.0 / var
   yin = wt * fluxobs
   flux = ndimage.filters.uniform_filter(yin,filtwidth)
   flux /= ndimage.filters.uniform_filter(wt,filtwidth)

   # Show some plots
   plt.figure(1)
   plt.clf()
   plt.plot(wave,fluxobs)
   plt.plot(wave,flux,'r')

   # Define the spectral range to be included in the fit
   if fitrange is not None:
      mask = np.zeros(flux.size,dtype=np.bool)
      fitr = np.atleast_2d(np.asarray(fitrange))
      for i in range(fitr.shape[0]):
         wmask = np.logical_and(wave>fitr[i,0],wave<fitr[i,1])
         mask[wmask] = True
      wavegood = wave[mask]
      fluxgood = flux[mask]
   else:
      wavegood = wave
      fluxgood = flux

   # Fit a polynomial to the observed response function
   fpoly = np.polyfit(wavegood,fluxgood,order)
   print ""
   print "Fit a polynomial of order %d to the red curve in Figure 1." % order
   print "Resulting coefficients:"
   print "-----------------------"
   print fpoly

   # Convert polynomial into a smooth response function
   p = np.poly1d(fpoly)
   cfit = p(wave)

   # Add the smooth response to the plot and show corrected curve
   plt.plot(wave,cfit,'k')
   fc = fluxobs / cfit
   vc = var / cfit**2
   plt.figure(2)
   plt.clf()
   plt.plot(wave,fc)

   # Write normalized spectrum to output file
   save_spectrum(outfile,wave,fc,vc)

#-----------------------------------------------------------------------

def atm_trans(w, fwhm=15., flux=None, scale=1.05, offset=0.0, modfile='default'):
   """
   Creates a Spec1d instance (i.e., a 1-dimensional spectrum) containing the 
   transmission of the Earth's atmosphere as a function of wavelength in
   the near-infrared (NIR) part of the spectrum.  The returned spectrum
   is for the wavelength range specified by the required w parameter, which
   is a wavelength vector.

   Inputs:
      w       - wavelength vector whose min and max values set the wavelength
                range of the returned atmospheric transmission spectrum
      fwhm    - smoothing parameter for the output spectrum
      modfile - the full path+name of the file containing the atmospheric
                transmission data.  The default location is in the Data
                subdirectory contained within the directory in which this
                code is found.
   """

   """ Read in the atmospheric transmission data"""
   if modfile != 'default':
      infile = modfile
   else:
      infile = '%s/Data/atm_trans_maunakea.fits' \
          % __file__.split('spec_simple')[0]
   print "Loading atmospheric data from %s" % infile
   atm0 = Spec1d(infile,informat='fitstab')
   atm0.wav *= 1.0e4

   """ Only use the relevant part of the atmospheric transmission spectrum"""
   mask = np.where((atm0.wav>=w.min())&(atm0.wav<=w.max()))
   watm = atm0.wav[mask]
   trans = atm0.flux[mask]
   del atm0
   
   """ Smooth the spectrum """
   trans = ndimage.gaussian_filter(trans,fwhm)
   
   """ Resample the smoothed spectrum """
   watm,trans = resample_spec(watm,trans,w)
   
   """ If an input spectrum has been given, then rescale the trans spectrum """
   if flux is not None:
      trans *= scale * flux.max()
   else:
      trans *= scale
   
   """ Add any requested vertical offset """
   trans += offset

   """ Save as a Spec1d instance and return """
   atm = Spec1d(wav=watm,flux=trans)
   del watm,trans
   return atm

#-----------------------------------------------------------------------

def plot_atm_trans(w, fwhm=15., flux=None, scale=1.05, offset=0.0,
                   color='g', linestyle='-', return_atm=False,
                   modfile='default', title=None):
   """
   Given an input spectrum, represented by the wavelength (w) and flux (spec)
   vectors, and a rough fwhm (in Angstrom), smooths and resamples the 
   atmospheric transmission spectrum for the NIR and plots it.
   """

   """ Get the atmospheric spectrum that matches the input wavelength array """
   atm = atm_trans(w,flux=flux,scale=scale,offset=offset,modfile=modfile)

   """ Plot the results """
   atm.plot(color=color,linestyle=linestyle,title=title)

   if return_atm:
      return atm
   else:
      del atm

#-----------------------------------------------------------------------

def plot_model_sky_ir():
   """
   Calls plot_atm_trans and make_sky_model to make a combined plot for
   the NIR sky that can be used to judge whether expected spectral lines
   will fall in good parts of the sky
   """

   wsky = np.arange(10000.,26000.)
   atm = atm_trans(wsky)
   atm.plot(color='g',title=None)
   skymod = make_sky_model(wsky)
   skymod /= skymod.max()
   plt.plot(wsky,skymod)

#-----------------------------------------------------------------------

def smooth_boxcar(infile, filtwidth, outfile=None, varwt=True):
   """
   Does a boxcar smooth of an input spectrum.  The default is to do
   inverse variance weighting, using the variance encoded in the third column
   of the input spectrum file.
   The other default is not to write out an output file.  This can be
   changed by setting the outfile parameter.
   """

   """ Read the input spectrum """
   inspec = np.loadtxt(infile)
   wavelength = inspec[:,0]
   influx = inspec[:,1]
   if(varwt):
      if(inspec.shape[1] < 3):
         print ""
         print "ERROR: Inverse variance weighting requested, but input file"
         print "       has fewer than 3 columns (ncol = %d)" % inspec.shape[1]
         return
      else:
         wt = 1.0 / inspec[:,2]
   else:
      wt = 0.0 * influx + 1.0

   """ Smooth spectrum """
   yin = wt * influx
   outflux = ndimage.filters.uniform_filter(yin,filtwidth)
   outflux /= ndimage.filters.uniform_filter(wt,filtwidth)
   if varwt:
      outvar = 1.0 / (filtwidth * ndimage.filters.uniform_filter(wt,filtwidth))

   """ Plot the smoothed spectrum """
   if varwt:
      plot_spectrum_array(wavelength,outflux,outvar,title=None)
   else:
      plot_spectrum_array(wavelength,outflux,title=None)

   """ Save the output file if desired """
   if(outfile):
      print "Saving smoothed spectrum to %s" % outfile
      if varwt:
         save_spectrum(outfile,wavelength,outflux,outvar)
      else:
         save_spectrum(outfile,wavelength,outflux)
      print ""

#-----------------------------------------------------------------------

def calc_lineflux(wavelength,flux,bluemin,bluemax,redmin,redmax,var=None,
                  showsub=False):
   """
   Given vectors of flux and wavelength, interactively calculates the integrated
   flux in an emission line.  The user enters the wavelength ranges to use
   for the continuum on both the blue (bluemin and bluemax) and red (redmin
   and redmax) sides of the line.  The function will do a first order fit 
   (i.e., a line) to the continuum using these ranges, subtract the continuum
   from the data, and then numerically integrate the flux/counts in the line.
   """

   """ Plot the data over this spectral range """
   mask = (wavelength>bluemin) & (wavelength<=redmax)
   tmplamb = wavelength[mask].copy()
   tmpflux = flux[mask].copy()
   if(var):
      tmpvar = var[mask].copy()
      plot_spectrum_array(tmplamb,tmpflux,tmpvar)
   else:
      plot_spectrum_array(tmplamb,tmpflux)

   """ Find a linear fit to the background regions """
   bkgdmask = ((tmplamb>bluemin) & (tmplamb<bluemax)) | \
       ((tmplamb>redmin) & (tmplamb<redmax))
   bkgdwave = tmplamb[bkgdmask].copy()
   bkgdflux = tmpflux[bkgdmask].copy()
   bkgdpoly = np.polyfit(bkgdwave,bkgdflux,1)
   continuum = tmplamb*bkgdpoly[0] + bkgdpoly[1]
   plt.plot(tmplamb,continuum,'r')
   plt.axvline(bluemin,color='k')
   plt.axvline(bluemax,color='k')
   plt.axvline(redmin,color='k')
   plt.axvline(redmax,color='k')
   plt.xlim(tmplamb[0],tmplamb[tmplamb.size - 1])

   """ Calculate the subtracted spectrum, and plot it if desired """
   subflux = tmpflux - continuum

   if(showsub):
      plt.figure()
      plt.clf()
      plt.plot(tmplamb,subflux)
      plt.xlim(tmplamb[0],tmplamb[tmplamb.size - 1])

   """ Numerically integrate the flux/counts in the line region """
   linemask = np.logical_not(bkgdmask)
   linewave = tmplamb[linemask].copy()
   lineflux = subflux[linemask].copy()
   # Assume that the wavelength scale is linear
   delwave = linewave[1] - linewave[0]
   print delwave
   intflux = (lineflux * delwave).sum()
   print intflux

#===========================================================================
#
# Rumbaugh code that may or may not get discarded at a later time
# (has now been commented out but not yet actually discarded)
#
#===========================================================================

##-----------------------------------------------------------------------
#
#def plot_multiple_peaks(cdat,tp,theight,apmin=-4.,apmax=4.,maxpeaks=2,fig=4,clearfig=True,plot_fits=True,apertures=None):
#   plt.figure(fig)
#   if clearfig: plt.clf()
#   plt.plot(np.arange(1,theight+1),cdat,linestyle='steps',color='black')
#   xmod = np.arange(1,theight+1,0.1)
#   tcolors = np.array(['red','cyan','magenta','green','blue','yellow'])
#   for ipg in range(0,maxpeaks):
#      ymod = make_gauss(xmod,tp[ipg][1],tp[ipg][2],tp[ipg][3],tp[0][0])
#      if plot_fits: plt.plot(xmod,ymod,color=tcolors[ipg])
#      if apertures == None:
#         plt.axvline(tp[ipg][1]+apmin,color=tcolors[ipg])
#         plt.axvline(tp[ipg][1]+apmax,color=tcolors[ipg])
#      else:
#         plt.axvline(tp[ipg][1]+apertures[ipg],color=tcolors[ipg])
#         plt.axvline(tp[ipg][1]-apertures[ipg],color=tcolors[ipg])
#      if tp[ipg][3]*1.05 > 2*np.max(cdat):
#         plt.text(tp[ipg][1],1.8*np.max(cdat),str(ipg+1),color=tcolors[ipg])
#      elif ((tp[ipg][3]*1.05 < 2*np.min(cdat)) & (tp[ipg][3]*1.05 < -2*np.max(cdat))):
#         plt.text(tp[ipg][1],np.min([1.8*np.min(cdat),-1.8*np.max(cdat)]),str(ipg+1),color=tcolors[ipg])
#      else:
#         plt.text(tp[ipg][1],tp[ipg][3]*1.05,str(ipg+1),color=tcolors[ipg])
#   if plt.ylim()[1] > 2*np.max(cdat):
#      plt.ylim(plt.ylim()[0],2*np.max(cdat))
#   if ((plt.ylim()[0] < 2*np.min(cdat)) & (plt.ylim()[0] < -2*np.max(cdat))):
#      plt.ylim(np.min([2*np.min(cdat),-2*np.max(cdat)]),plt.ylim()[1])
#   plt.xlabel('Pixel number in the spatial direction')
#   plt.title('Compressed Spatial Plot with Potential Peaks')
#
##-----------------------------------------------------------------------
#
#def find_multiple_peaks(data,dispaxis="x",apmin=-4.,apmax=4.,maxpeaks=2,output_plot=None,output_plot_dir=None,check_aps=False):
#   tdata = data.copy()
#   gfbc = find_blank_columns(tdata)
#   if dispaxis == 'x':
#      data[:,gfbc]
#   tp = np.zeros((maxpeaks,4,))
#   p_prelim = find_peak(tdata,dispaxis=dispaxis,apmin=apmin,apmax=apmax,showplot=False,do_subplot=False,nofit=True)
#   tp[0] = p_prelim
#   for ifmp in range(1,maxpeaks):
#      if dispaxis == 'x':
#         theight = np.shape(tdata[:,gfbc])[0]
#         tlength = np.shape(tdata[:,gfbc])[1]
#         if ifmp == 1: x = np.arange(1,theight+1)
#         gx = np.where((x < p_prelim[1]+2*apmin) | (x > p_prelim[1]+2*apmax))[0]
#         if ifmp != 1: gx = np.intersect1d(gx,gxprev)
#         p_prelim = find_peak(tdata[gx,:],dispaxis=dispaxis,apmin=apmin,apmax=apmax,showplot=False,do_subplot=False,nofit=True)
#         tp[ifmp] = p_prelim
#         tp[ifmp][1] = x[gx[int(tp[ifmp][1])-1]]
#         gxprev = gx.copy()
#      else:
#         tlength = np.shape(tdata[gfbc,:])[0]
#         theight = np.shape(tdata[gfbc,:])[1]
#         if ifmp == 1: x = np.arange(1,theight+1)
#         gx = np.where((x < p_prelim[1]+2*apmin) | (x > p_prelim[1]+2*apmax))[0]
#         if ifmp != 1: gx = np.intersect1d(gx,gxprev)
#         p_prelim = find_peak(tdata[:,gxprev],dispaxis=dispaxis,apmin=apmin,apmax=apmax,showplot=False,do_subplot=False,nofit=True)
#         tp[ifmp] = p_prelim
#         tp[ifmp][1] = x[gx[int(tp[ifmp][1])-1]]
#         gxprev = gx.copy()
#   tp = find_peak(tdata,dispaxis=dispaxis,apmin=apmin,apmax=apmax,showplot=False,do_subplot=False,mu0=tp[:,1])
#   tp[0] = np.ones(len(tp[1]))*tp[0]
#   tp = np.transpose(tp)
#   if dispaxis == 'x':
#      cdat = np.median(data[:,gfbc],axis=1)
#   else:
#      cdat = np.median(data[gfbc,:],axis=0)
#   plot_multiple_peaks(cdat,tp,theight,apmin=apmin,apmax=apmax,maxpeaks=maxpeaks)
#   print 'Plotting %i highest peaks found\n'%maxpeaks
#   tflag,fitmp,fixmu = False,False,False
#   while not tflag:
#      inp_fitmp = raw_input('Reduce secondary peaks? (y/n)\n')
#      if ((inp_fitmp == 'y') | (inp_fitmp == 'Y')):
#         tflag,fitmp = True,True
#      elif ((inp_fitmp == 'n') | (inp_fitmp == 'N')):
#         tflag = True
#      elif inp_fitmp == 'fixmu':
#         tflag,fixmu = True,True
#      else:
#         print 'Invalid input\n'
#   fitpeaks = np.zeros(maxpeaks,dtype='bool')
#   fitpeaks[0] = True
#   bounds_arr = np.array([0,np.min(np.shape(data))])
#   if fitmp:
#      tflag = False
#      while not tflag:
#         inp_chp1 = raw_input("Is peak 1 okay? (y/n)\n")
#         if ((inp_chp1 == 'y') | (inp_chp1 == 'Y')):
#            tflag = True
#         elif((inp_chp1 == 'n') | (inp_chp1 == 'N')):
#            mflag = False
#            while not mflag:
#               inp_newp = raw_input("Current mu for peak 1 is %f. Is this acceptable? Enter 'y' or new value for mu.\n"%(tp[0][1]))
#               if ((inp_newp == 'y') | (inp_newp == 'Y')):
#                  mflag,tflag = True,True
#               else:
#                  try:
#                     tp[0][1] = float(inp_newp)
#                     plot_multiple_peaks(cdat,tp,theight,apmin=apmin,apmax=apmax,maxpeaks=maxpeaks)
#                  except ValueError:
#                     print 'Invalid input\n'
#         else:
#            print 'Invalid input\n'
#      for iwp in range(1,maxpeaks):
#         tflag = False
#         while not tflag:
#            inp_whichp = raw_input("Reduce peak %i? (y/n/manual) Enter 'manual' to manually set peak position\n"%(iwp+1))
#            if ((inp_whichp == 'y') | (inp_whichp == 'Y')):
#               tflag,fitpeaks[iwp] = True,True
#            elif((inp_whichp == 'n') | (inp_whichp == 'N')):
#               tflag = True
#            elif ((inp_whichp == 'manual') | (inp_whichp == 'm')):
#               mflag = False
#               while not mflag:
#                  inp_newp = raw_input("Current mu for peak %i is %f. Is this acceptable? Enter 'y' or new value for mu.\n"%(iwp+1,tp[iwp][1]))
#                  if ((inp_newp == 'y') | (inp_newp == 'Y')):
#                     mflag = True
#                     tflag,fitpeaks[iwp] = True,True
#                  else:
#                     try:
#                        tp[iwp][1] = float(inp_newp)
#                        plot_multiple_peaks(cdat,tp,theight,apmin=apmin,apmax=apmax,maxpeaks=maxpeaks)
#                     except ValueError:
#                        print 'Invalid input\n'
#            else:
#               print 'Invalid input\n'
#      num_peaks = len(fitpeaks[fitpeaks])
#      mp_out = np.zeros((4,num_peaks))
#      for impo in range(0,maxpeaks): 
#         if fitpeaks[impo]: 
#            inow = len(fitpeaks[0:impo+1][fitpeaks[0:impo+1]])
#            mp_out[:,inow-1] = tp[inow-1]
#      mus_tmp = mp_out[1]
#      sort_mus = np.sort(mus_tmp)
#      argsort_mus = np.argsort(mus_tmp)
#      tbounds_arr = np.zeros(2*num_peaks)
#      tbounds_arr[2*num_peaks-1] = np.min(np.shape(data))
#      for il in range(0,num_peaks-1): tbounds_arr[2*il+1:2*il+3] = np.mean(sort_mus[il:il+2])
#      bounds_arr = np.zeros(2*num_peaks)
#      aa_mus = np.argsort(argsort_mus)
#      for il in range(0,num_peaks): bounds_arr[2*il:2*il+2] = tbounds_arr[2*aa_mus[il]:2*aa_mus[il]+2]
#      plot_multiple_peaks(cdat,tp,theight,apmin=apmin,apmax=apmax,maxpeaks=num_peaks)
#      for il in range(0,2*num_peaks): plt.axvline(bounds_arr[il],color='k')
#      
#      tflag = False
#      inp_aps = raw_input("Are these bounds okay? (y/n)\n")
#      while not tflag:
#         if ((inp_aps == 'y') | (inp_aps == 'Y')):
#            tflag = True
#         elif ((inp_aps == 'n') | (inp_aps == 'N')):
#            for ilf in range(0,num_peaks):
#               tflag2 = False
#               inp_aps2 = raw_input("Are the bounds for peak %i okay? (y/n)\n"%(ilf+1))
#               while not tflag2:
#                  if ((inp_aps2 == 'y') | (inp_aps2 == 'Y')):
#                     tflag2 = True
#                  elif ((inp_aps2 == 'n') | (inp_aps2 == 'N')):
#                     tflag3 = False
#                     nlb = raw_input("Current bounds for peak %i are (%.1f,%.1f). Enter new lower bound:\n"%(ilf+1,bounds_arr[2*ilf],bounds_arr[2*ilf+1]))
#                     nub = raw_input('Enter new upper bound:\n')
#                     while not tflag3:
#                        try: 
#                           nlb,nub = float(nlb),float(nub)
#                           if ((nlb < 0) | (nub > np.min(np.shape(data))) | (nub <= nlb)): raise ValueError
#                           bounds_arr[2*ilf],bounds_arr[2*ilf+1] = nlb,nub
#                           plot_multiple_peaks(cdat,tp,theight,apmin=apmin,apmax=apmax,maxpeaks=num_peaks)
#                           for ilt in range(0,2*num_peaks): plt.axvline(bounds_arr[ilt],color='k')
#                           tflag3 = True
#                           inp_aps2 = raw_input("New bounds for peak %i are: (%.1f,%.1f). Are these okay? (y/n)\n"%(ilf+1,bounds_arr[2*ilf],bounds_arr[2*ilf+1]))
#                        except ValueError:
#                           print 'Invalid input. Bounds must be floats between 0 and %.1f.\n'%(np.min(np.shape(data)))
#                           nlb = raw_input("Current bounds for peak %i are (%.1f,%.1f). Enter new lower bound:\n"%(ilf+1,bounds_arr[2*ilf],bounds_arr[2*ilf+1]))
#                           nub = raw_input('Enter new upper bound:\n')
#                  else:
#                     print 'Invalid input\n'
#                     inp_aps2 = raw_input("Are the bounds for peak %i okay? (y/n)\n"%(ilf+1))
#            tflag = True
#         else:
#            print 'Invalid input\n'
#            inp_aps = raw_input("Are these bounds okay? (y/n)\n")
#   try:
#      bnds_bool = (bounds_arr == np.array([0,np.min(np.shape(data))])).all()
#   except AttributeError:
#      bnds_bool = (bounds_arr == np.array([0,np.min(np.shape(data))]))
#   if ((fitpeaks[0]) & (len(fitpeaks[fitpeaks]) == 1) & bnds_bool): 
#      fitmp = False
#      print 'No secondary peaks selected. Reverting to normal analysis.'
#   num_peaks = len(fitpeaks[fitpeaks])
#   mp_out = np.zeros((4,num_peaks))
#   for impo in range(0,maxpeaks): 
#      if fitpeaks[impo]: 
#         inow = len(fitpeaks[0:impo+1][fitpeaks[0:impo+1]])
#         mp_out[:,inow-1] = tp[inow-1]
#   aflag,change_aps = False,False
#   if check_aps:
#      while not aflag:
#         inp_aps = raw_input("Change apertures? (y/n)\n")
#         if ((inp_aps == 'y') | (inp_aps == 'Y')):
#            aflag,change_aps = True,True
#         elif ((inp_aps == 'n') | (inp_aps == 'N')):
#            aflag = True
#         else:
#            print 'Invalid input.\n'
#   apertures = 4.*np.ones(num_peaks)
#   if change_aps:
#      for iaps in range(0,num_peaks):
#         aflag = False
#         while not aflag:
#            inp_aps = raw_input("Change apertures for peak %i? (y/n)\n"%(iaps+1))
#            if ((inp_aps == 'y') | (inp_aps == 'Y')):
#               aflag2 = False
#               while not aflag2:
#                  inp_aps2 = raw_input("Aperture for peak %i is +%.1f,-%.1f. Is this okay? Enter 'y' or new width.\n"%(iaps+1,apertures[iaps],apertures[iaps]))
#                  if ((inp_aps2 == 'y') | (inp_aps2 == 'Y')):
#                     aflag2 = True
#                  else:
#                     try:
#                        if inp_aps > 0: 
#                           apertures[iaps] = inp_aps2
#                           plot_multiple_peaks(cdat,tp,theight,apmin=-1*apertures[iaps],apmax=apertures[iaps],maxpeaks=num_peaks,apertures=apertures)
#                        else:
#                           print 'Input value must be greater than zero.'
#                     except:
#                        print 'Invalid input'
#               aflag = True
#            elif ((inp_aps == 'n') | (inp_aps == 'N')):
#               aflag = True
#            else:
#               print 'Invalid input.\n'
#   if output_plot != None:
#      outplotname = 'bounds.%s'%output_plot
#      if output_plot_dir != None: outplotname = '%s/%s'%(output_plot_dir,outplotname)
#      plot_multiple_peaks(cdat,np.transpose(mp_out),theight,apmin=apmin,apmax=apmax,maxpeaks=num_peaks,plot_fits=False,apertures=apertures)
#      for il in range(0,2*num_peaks): plt.axvline(bounds_arr[il],color='k')
#      plt.title('Compressed Spatial Plot with Extraction Regions')
#      plt.savefig(outplotname)
#   if check_aps:
#      if num_peaks == 1:
#         return False,fixmu,tp[0],bounds_arr,apertures
#      else:
#         return fitmp,fixmu,mp_out,bounds_arr,apertures
#   else:
#      if num_peaks == 1:
#         return False,fixmu,tp[0],bounds_arr
#      else:
#         return fitmp,fixmu,mp_out,bounds_arr

