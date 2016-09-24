"""
imfuncs.py - A library of functions to do various basic image processing
             operations

NB: Some of these functions are (slowly) being incorporated into the
Image class, the code for which is at the beginning of this file.

Methods (functions) in the Image class
--------------------------------------

Stand-alone functions
---------------------
   open_fits         - opens a fits file, incorporating the possibility of
                       the "missing end" problem that affects some of the Keck
                       NIR data
   imcopy            - copies a portion of a fits file given the corners (in 
                       pixels)
   poststamp         - copies a portion of a fits file given the center and size
                       (in pixels)
   image_cutout      - copies a portion of a fits file given the center and
                       size (center in RA,Dec)
   make_snr_image    - calculates the rms noise in an image and divides by
                       that value to create a SNR image, which can be written
                       as an output fits file.
   overlay_contours  - plots a greyscale image and overlays contours from a
                       second image
   calc_sky_from_seg - calculates sky level using SExtractor segmentation map
   display_image     - displays an image
   plot_cat          - given a fits image and a object catalog, marks the
                       positions of the catalog objects.
"""

import astropy
try:
   from astropy.io import fits as pf
except:
   import pyfits as pf
from astropy import wcs
if astropy.__version__[:3] == '0.3':
   from astropy.coordinates import ICRS as SkyCoord
else:
   from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as n
from scipy import ndimage
import matplotlib.pyplot as plt
from math import log,log10,sqrt,pi,fabs
from math import cos as mcos,sin as msin
import wcs as wcsmwa
import ccdredux as ccd

#-----------------------------------------------------------------------

class Image:

   def __init__(self, infile, mode='copyonwrite', verbose=True):
      """
      This method gets called when the user types something like
         im = Image(infile)

      Inputs:
         infile - input file containing the fits image
         mode   - mode in which to open the fits image.  The default value,
                  'copyonwrite' is the standard read-only mode.  To enable
                  modification in place, use mode='update'
      """

      """
      Start by loading the hdu information
      """
      if verbose:
         print ""
         print "Loading file %s" % infile
         print "-----------------------------------------------"
      try:
         self.hdu = open_fits(infile, mode)
      except:
         print "  ERROR. Problem in loading file %s" % infile
         print "  Check to make sure filename matches an existing file."
         print "  If it does, there may be something wrong with the fits header."
         print ""
         return
      if verbose:
         self.hdu.info()

      """ Set parameters related to image properties """
      self.infile = infile
      self.fitsmode = mode

      """ Do an initial import of the WCS information from the header """
      self.found_wcs = False
      self.radec = None
      self.ra  = None
      self.dec = None
      self.pixscale = None
      try:
         self.get_wcs()
      except:
         print 'Could not load WCS info for %s' % self.infile
         self.found_wcs = False

      """ Initialize figures """
      self.fig1 = None
      self.fig2 = None

      """ 
      Initialize default display parameters 

       - The scale for the display (i.e., the data values that correspond
         to full black and full white on a greyscale display) are (by default) 
         set in terms of the "clipped mean" and "clipped rms".  Those values
         are the mean and rms of the data after a sigma clipping algorithm
         has been applied to reject outliers.
       - The display min and max values are stored as self.fmin and self.fmax
       - For more information see the set_display_limits method
      """
      self.found_rms = False    # Have the clipped rms and mean been calculated?
      self.mean_clip = 0.0      # Value of the clipped mean
      self.rms_clip  = 0.0      # Value of the clipped rms
      self.fmin      = None     # Lower flux limit used in image display
      self.fmax      = None     # Upper flux limit used in image display
      self.fscale    = 'linear' # Flux scaling for display
      self.statsize  = 2048     # Region size for statistics if image is too big
      self.zoomsize  = 31       # Size of postage-stamp zoom
      self.dispunits = 'radec'  # Default display units are arcsec offsets
      self.extval    = None     # Just label the axes by pixels
      self.cmap      = plt.cm.YlOrBr_r # This corresponds to the 'gaia' cmap

      """ Initialize contouring parameters """
      self.contbase   = sqrt(3.)
      self.clevs      = None
      self.overlay_im = None  # Not currently used

      """ Initialize other parameters """
      self.reset_subim()
      self.reset_imex()

   #-----------------------------------------------------------------------

   def reset_subim(self):
      """
      Returns the sub-image variables to their initial, unset, state
      """

      self.subim = None
      self.subx1 = None
      self.subx2 = None
      self.suby1 = None
      self.suby2 = None
      self.subsizex = None
      self.subsizey = None
      self.subcentx = None
      self.subcenty = None

   #-----------------------------------------------------------------------

   def reset_imex(self):
      """
      Resets the parameters that are associated with the imexam-like
       processing
      """

      """ Initialize coordinates and data """
      self.imex_x = None
      self.imex_y = None
      self.imex_data = None

      """ Initialize moment calculations """
      self.imex_mux   = None
      self.imex_muy   = None
      self.imex_sigxx = None
      self.imex_sigyy = None
      self.imex_sigxy = None

   #-----------------------------------------------------------------------

   def close(self):
      """
      Closes the image
      """

      self.hdu.close()
      return

   #-----------------------------------------------------------------------

   def sigma_clip(self, nsig=3., mask=None, hext=0, verbose=False):
      """
      Runs a sigma-clipping on image data.  The code iterates over
      the following steps until it has converged:
         1. Compute mean and rms noise in the (clipped) data set
         2. Reject (clip) data points that are more than nsig * rms from
            the newly computed mean (using the newly computed rms)
         3. Repeat until no new points are rejected
      Once convergence has been reached, the final clipped mean and clipped
      rms are stored in the mean_clip and rms_clip variables

      Optional inputs:
         nsig    - Number of sigma from the mean beyond which points are
                    rejected.  Default=3.
         mask    - If some of the input data are known to be bad, they can
                    be flagged before the inputs are computed by including
                    a mask.  Clearly his mask must be set such that True
                    indicates good data and False indicates bad data
         hext    - Image HDU containing the data.  The default (hext=0) should
                    work for all single-extension fits files and may 
                    work for some multi-extension files.
                    NOTE: hext is ignored if the subim variable is already set
                    by, e.g., def_subim_xy or def_subim_radec
         verbose - If False (the default) no information is printed
      """

      """ Determine what the input data set is """
      if self.subim is not None:
         if mask:
            data = self.subim[mask]
         else:
            data = self.subim
      else:
         if mask:
            data = self.hdu[hext].data[mask]
         else:
            data = self.hdu[hext].data

      """ Report the number of valid data points """
      size = data[n.isfinite(data)].size 
      if verbose:
         print " sigma_clip: Full size of data       = %d" % data.size
         print " sigma_clip: Number of finite values = %d" % size

      """ Reject the non-finite data points and compute the initial values """
      d = data[n.isfinite(data)].flatten()
      mu = d.mean()
      sig = d.std()
      mu0 = d.mean()
      sig0 = d.mean()
      if verbose:
         print ''
         print 'npix = %11d. mean = %f. sigma = %f' % (size,mu,sig)
   
      """ Iterate until convergence """
      delta = 1
      clipError = False
      while delta:
         size = d.size
         if sig == 0.:
            clipError = True
            break
         d = d[abs(d-mu)<nsig*sig]
         mu = d.mean()
         sig = d.std()
         if verbose:
            print 'npix = %11d. mean = %f. sigma = %f' % (size,mu,sig)
         delta = size-d.size
      if clipError:
         print ''
         print 'ERROR: sigma clipping failed, perhaps with rms=0.'
         print 'Setting mu and sig to their original, unclipped, values'
         print ''
         mu = mu0
         sig = sig0

      """ Store the results and clean up """
      del data,d
      self.mean_clip = mu
      self.rms_clip  = sig
      return

   #-----------------------------------------------------------------------

   def start_interactive(self):
      self.xmark = None
      self.ymark = None
      self.cid_mouse = self.fig1.canvas.mpl_connect('button_press_event',
                                                    self.onclick)
      self.cid_keypress = self.fig1.canvas.mpl_connect('key_press_event',
                                                       self.keypress)
      self.keypress_info()
      return

   #-----------------------------------------------------------------------

   def keypress_info(self):
      """
      Prints useful information about what the key presses do
      """
      print ''
      print 'Actions available by pressing a key in the Figure 1 window'
      print '----------------------------------------------------------'
      print 'Key      Action'
      print '-------  ---------------------------------'
      print '[click]  Report (x,y) position, and (RA,Dec) if file has WCS'
      print '   m     Mark the position of an object'
      print '   z     Zoom in at the position of the cursor'
      print '   q     Quit and close the window'
      print ''

   #-----------------------------------------------------------------------

   def onclick(self, event):
      """
      Actions taken if a mouse button is clicked.  In this case the
      following are done:
        (1) Store and print (x,y) value of cursor
        (2) If the image has wcs info (i.e., if found_wcs is True) then
            store and print the (RA,Dec) value associated with the (x,y)
      """
      self.xclick = event.xdata
      self.yclick = event.ydata
      print ''
      print 'Mouse click x, y:   %7.1f %7.1f' % (self.xclick,self.yclick)

      """ 
      Also show the (RA,Dec) of the clicked position if the input file has WCS
      NOTE: This needs to be handled differently if the displayed image has
       axes in pixels or in arcsec offsets
      """
      if self.found_wcs:
         if self.dispunits == 'pixels':
            pix = n.zeros((1,self.wcsinfo.naxis))
            pix[0,0] = self.xclick
            pix[0,1] = self.yclick
            radec = self.wcsinfo.wcs_pix2world(pix,0)
            self.raclick  = radec[0,0]
            self.decclick = radec[0,1]
         else:
            """ For now use small-angle formula """
            cosdec = mcos(self.radec.dec.radian)
            self.raclick = self.radec.ra.degree + \
                (self.xclick + self.zeropos[0]) / (3600. * cosdec)
            self.decclick = self.radec.dec.degree + self.yclick/3600. + \
                self.zeropos[1]
         print 'Mouse click ra, dec: %11.7f %+11.7f' % \
             (self.raclick,self.decclick)
      return

   #-----------------------------------------------------------------------

   def keypress(self, event):
      """
      Actions taken if a key on the keyboard is pressed
      """

      if event.key == 'f':
         """
         Change the display range
         """
         print ''
         self.set_display_limits(fmax=None,funits='abs')
         
      if event.key == 'm':
         """
         Mark an object.  Hitting 'm' saves the (x,y) position into
         the xmark and ymark variables
         """
         global xmark, ymark
         print ''
         print 'Marking position %8.2f %8.2f' % (event.xdata,event.ydata)
         print ''
         self.xmark = event.xdata
         self.ymark = event.ydata
         subimsize = (self.zoomsize,self.zoomsize)
         subimcent = (self.xmark,self.ymark)
         self.display(subimcent=subimcent,subimsize=subimsize,show_xyproj=True)

      if event.key == 'z':
         """
         Zoom in by a factor of two at the location of the cursor
         """
         xzoom,yzoom = event.xdata,event.ydata
         xl1,xl2 = self.ax1.get_xlim()
         yl1,yl2 = self.ax1.get_ylim()
         dx = (xl2 - xl1)/4.
         dy = (yl2 - yl1)/4.
         xz1 = min((max(xl1,(xzoom-dx))),(xzoom-1.))
         xz2 = max((min(xl2,(xzoom+dx))),(xzoom+1.))
         yz1 = min((max(yl1,(yzoom-dy))),(yzoom-1.))
         yz2 = max((min(yl2,(yzoom+dy))),(yzoom+1.))
         self.ax1.set_xlim(xz1,xz2)
         self.ax1.set_ylim(yz1,yz2)
         self.fig1.show()
         return

      if event.key == 'q':
         print ''
         print 'Closing down'
         print ''
         if self.fig1:
            self.fig1.canvas.mpl_disconnect(self.cid_mouse)
            self.fig1.canvas.mpl_disconnect(self.cid_keypress)
         if self.fig2:
            self.fig2.canvas.mpl_disconnect(self.cid_keypress2)
         for ii in plt.get_fignums():
            plt.close(ii)
         return

      self.keypress_info()
      return

   #-----------------------------------------------------------------------

   def radec_to_skycoord(self, ra, dec):
      """
      Converts a (RA,Dec) pair into the astropy.coordinates SkyCoord format

      Required inputs:
        ra   - RA in one of three formats:
                Decimal degrees:  ddd.ddddddd  (as many significant figures
                  as desired)
                Sexigesimal:      hh mm ss.sss (as many significant figures
                  as desired)
                Sexigesimal:      hh:mm:ss.sss (as many significant figures
                  as desired)
        dec  - Dec in one of three formats:
                Decimal degrees:  sddd.ddddddd, where "s" is + or -
                Sexigesimal:      sdd mm ss.sss (as many significant figures
                  as desired)
                Sexigesimal:      sdd:mm:ss.sss (as many significant figures
                  as desired)
      """

      """ Get RA format """
      if type(ra)==float or type(ra)==n.float32 or type(ra)==n.float64:
         rafmt = u.deg
      else:
         rafmt = u.hourangle

      """ Dec format is always degrees, even if in Sexigesimal format """
      decfmt = u.deg

      """ Do the conversion """
      self.radec = SkyCoord(ra,dec,unit=(rafmt,decfmt))

   #-----------------------------------------------------------------------

   def get_wcs(self, hext=0):
      """
      Reads in WCS information from the header and stores it in 
      new wcsinfo (see below) and pixscale variables.
      NOTE: This used to use Matt Auger's wcs library, but it has
       been converted to the astropy wcs module

      Inputs:
         hext - Header extension that contains the WCS info.  Default=0

      """

      """ Read in the header and use it to set the WCS information"""
      hdr = self.hdu[hext].header
      try:
         self.wcsinfo = wcs.WCS(hdr)
      except:
         print 'get_wcs: No WCS information in file header'
         self.found_wcs = False
         return

      """ 
      Make sure that the WCS information is actually WCS-like and not,
      for example, pixel-based
      """

      if self.wcsinfo.wcs.ctype[0][0:2] != 'RA':
         print 'get_wcs: No valid WCS information in file header'
         self.found_wcs = False
         return

      """ Get the RA and Dec of the center of the image """
      xcent = hdr['naxis1'] / 2.
      ycent = hdr['naxis2'] / 2.
      imcent = n.ones((1,hdr['naxis']))
      imcent[0,0] = xcent
      imcent[0,1] = ycent
      imcentradec = self.wcsinfo.wcs_pix2world(imcent,1)
      self.radec_to_skycoord(imcentradec[0,0],imcentradec[0,1])

      """ Calculate the pixel scale """
      if self.wcsinfo.wcs.ctype[0][0:2].upper() == 'RA':
         try:
            self.pixscale = sqrt(self.wcsinfo.wcs.cd[0,0]**2 + 
                                 self.wcsinfo.wcs.cd[1,0]**2)*3600.
         except:
            try:
               self.pixscale = sqrt(self.wcsinfo.wcs.pc[0,0]**2 + 
                                    self.wcsinfo.wcs.pc[1,0]**2) * \
                                    self.wcsinfo.cdelt[0] * 3600.
            except:
               self.pixscale = abs(self.wcsinfo.wcs.cdelt[0]) * 3600.
         print 'Pixel scale: %7.3f arcsec/pix' % self.pixscale
         self.found_wcs = True
      else:
         print 'Warning: no WCS info in header %d' % hext
         self.found_wcs = False

   #-----------------------------------------------------------------------

   def set_pixscale(self):
      """
      Interactively set the pixel scale
      """

      print ''
      self.pixscale = \
          float(raw_input('Enter the pixel scale for the image in arcsec/pix'))

   #-----------------------------------------------------------------------

   def set_wcsextent(self, hext=0, zeropos=None):
      """
      For making plots with WCS information, it is necessary to define
      the boundaries in terms of RA and Dec offsets from the center, in
      arcsec.  For this purpose, the imshow and contour methods in 
      matplotlib.pyplot have an 'extent' parameter.  
      
      This set_wcsextent method will use the WCS information in the fits
      header to properly set the extent parameter values and return them.
      These are put into the "extval" container, which is part of the Image
      class.  extval is a four-element tuple containing the coordinates of the
      lower left and upper right corners, in terms of RA and Dec offsets.

      Optional inputs:
         hext    - HDU containing the WCS info.  Default=0
         zeropos - By default, which happens when zeropos=None, the (0,0)
                      point on the output image, as designated by the image
                      axis labels, will be at the center of the image.  However,
                      you can shift the (0,0) point to be somewhere else by
                      setting zeropos.  For example, zeropos=(0.5,0.3) will
                      shift the origin to the point that would have been
                      (0.5,0.3) if the origin were at the center of the image
      """

      self.get_wcs(hext)
      coords = n.indices(self.subim.shape).astype(n.float32)
      pltc = n.zeros(coords.shape)
      pltc[0] = (coords[0] - self.subim.shape[0]/2.)*self.pixscale
      pltc[1] = (coords[1] - self.subim.shape[1]/2.)*self.pixscale
      pltc[1] *= -1.
      maxi = n.atleast_1d(self.subim.shape) - 1
      extx1 = pltc[1][0,0]
      exty1 = pltc[0][0,0]
      extx2 = pltc[1][maxi[1],maxi[1]]-self.pixscale
      # *** IS THIS A BUG? SHOULD THE FF USE maxi[0]??? ***
      exty2 = pltc[0][maxi[1],maxi[1]]+self.pixscale

      if zeropos is not None:
         dx = zeropos[0]
         dy = zeropos[1]
      else:
         dx = 0.
         dy = 0.
      extx1 -= dx
      extx2 -= dx
      exty1 -= dy
      exty2 -= dy

      """ Set the extval values, and also record the zerpos values used """
      self.extval  = (extx1,extx2,exty1,exty2)
      self.zeropos = (dx,dy)

   #-----------------------------------------------------------------------

   def im_moments(self, x0, y0, rmax=10., detect_thresh=3., skytype='global', 
                  hext=0, verbose=False):
      """
      Given an initial guess of a centroid position, calculates the
      flux-weighted first and second moments within a square centered
      on the initial guess point and with side length of 2*rmax + 1.
      The moments will be estimates of the centroid and sigma of the
      light distribution within the square.

      Inputs:
        x0      - initial guess for x centroid
        y0      - initial guess for y centroid
        rmax    - used to set size of image region probed, which will be a
                   square with side = 2*rmax + 1. Default=10
        skytype - set how the sky/background level is set.  Three options:
                    'global' - Use the clipped mean as determined by the
                               sigma_clip method.  This is the default.
                    'local'  - Use a region that surrounds the source
                               NOT YET IMPLEMENTED
                    None     - Don't do any sky/background subtraction
        hext    - HDU extension that contains the data.  Default = 0
      """

      """ Define the data and the coordinate arrays """
      data = self.hdu[hext].data
      y,x = n.indices(data.shape)

      """ 
      Select the data within the square of interest
      """
      x1,x2 = x0-rmax-1,x0+rmax+1
      y1,y2 = y0-rmax-1,y0+rmax+1
      pixmask = (x>x1)&(x<x2)&(y>y1)&(y<y2)
      if skytype is None:
         f   = data[pixmask]
      else:
         if self.found_rms == False:
            self.sigma_clip(verbose=verbose)
            self.found_rms = True
         print self.mean_clip, self.rms_clip
         f   = data[pixmask] - self.mean_clip
      self.imex_x = x[pixmask]
      self.imex_y = y[pixmask]

      """ 
      Calculate the flux-weighted moments 
       NOTE: Do the moment calculations relative to (x1,y1) -- and then add 
        x1 and y1 back at the end -- in order to avoid rounding errors (see
        SExtractor user manual)
      """
      objmask = f > self.mean_clip + detect_thresh * self.rms_clip
      fgood = f[objmask]
      """ 
      """
      xgood = self.imex_x[objmask] - x1 
      ygood = self.imex_y[objmask] - y1
      fsum = fgood.sum()
      mux = (fgood * xgood).sum() / fsum
      muy = (fgood * ygood).sum() / fsum
      self.imex_mux = mux + x1
      self.imex_muy = mux + y1
      self.imex_sigxx = (fgood * xgood**2).sum() / fsum - mux**2
      self.imex_sigyy = (fgood * ygood**2).sum() / fsum - muy**2
      self.imex_sigxy = (fgood * xgood*ygood).sum() / fsum - mux*muy
      print self.imex_mux,self.imex_muy
      print sqrt(self.imex_sigxx),sqrt(self.imex_sigyy),self.imex_sigxy

   #-----------------------------------------------------------------------

   def radplot(self, x0, y0, rmax, skylevel=0. , zp=None, runit='pixel', 
               logr=False, hext=0):
      """
      Given a position in the image file (the x0 and y0 parameters), makes
       a plot of image flux as a function of distance from that (x,y) 
       position, out to a maximum distance of rmax.
      The default is to make the plot in flux units (or ADU or counts or
       counts/sec).  However, if the zero point parameter (zp) is set
       then the values in the data array will be converted into surface
       brightness in magnitudes, via the usual formula:
         mu = -2.5 log10(data) + zp
         
      Required inputs:
        x0   - x coordinate 
        y0   - y coordinate
        rmax - maximum radius, in pixels, for the plot
      Optional inputs:
        skylevel - If the sky has not been subtracted from the data, then
                    the integrated counts, surface brightness in mag/arcsec**2,
                    and integrated magnitude will all be wrong.  Set this
                    parameter to the rough sky level to address these issues.
                    The default (skylevel=0) is appropriate if the sky _has_
                    been subtracted.
        zp       - zero point.  If this parameter is set, then the output plot
                    will be in magnitude units (i.e., surface brightness) rather
                    than the default flux-like units (ADU, counts, 
                    counts/sec, etc.)
        runit    - units for the x-axis of the plot.  The only options are
                    'pixel' (the default) or 'arcsec'
        logr     - If False (the default) then x-axis is linear. If true, then 
                    it is in log
        hext     - HDU extension that contains the data.  Default = 0
      """

      """ Define the data and the coordinate arrays """
      data = self.hdu[hext].data
      y,x = n.indices(data.shape)

      """ 
      Find the offsets of each pixel in the data array from the requested
       central location (x0,y0).
      For better speed, only actually do the computations for pixels that
       might be in the correct area
      """
      x1,x2 = x0-rmax-1,x0+rmax+1
      y1,y2 = y0-rmax-1,y0+rmax+1
      pixmask = (x>x1)&(x<x2)&(y>y1)&(y<y2)
      dx = x[pixmask] - x0
      dy = y[pixmask] - y0
      r = n.sqrt(dx**2 + dy**2)

      """ 
      Get the pixel scale if needed, which it will be if either runit==arcsec
       or if zp is set (and thus the first plot is mag/arcsec**2).
      """
      if zp or (runit=='arcsec'):
         if self.pixscale is None:
            self.set_pixscale()
         print 'Using pixel scale of %6.3f arcsec/pix' % self.pixscale

      """ Select the points within rmax and convert to mags if desired """
      ii = n.argsort(r)
      if runit=='arcsec':
         rr = r[ii] * self.pixscale
         xlab = 'r (arcsec)'
         rmax *= self.pixscale
      else:
         rr = r[ii]
         xlab = 'r (pixels)'
      rflux = (data[pixmask])[ii] - skylevel
      if zp:
         domega = self.pixscale**2
         mu = -2.5 * n.log10(rflux/domega) + zp
         ftype = 'Surface Brightness'
         ttype = 'Magnitude'
         flab = 'Surface rrightness (mag/arcsec**2)'
         tlab = 'Integrated magnitude within r (mag)'
      else:
         ftype = 'Counts'
         ttype = 'Counts'
         flab = 'Counts / pixel'
         tlab = 'Integrated counts within r'

      """ Plot the surface brightness """
      ax1 = plt.subplot(211)
      #plt.setp(ax1.get_xticklabels(), visible=False)
      #plt.figure(1)
      if zp:
         if logr:
            plt.semilogx(rr,mu,'+')
         else:
            plt.plot(rr,mu,'+')
         yl1,yl2 = plt.ylim()
         plt.ylim(yl2,yl1)
      else:
         if logr:
            plt.semilogx(rr,rflux,'+')
         else:
            plt.plot(rr,rflux,'+')
      plt.xlim(0,rmax)
      plt.title('%s Profile centered at (%6.1f,%6.1f)'%(ftype,x0,y0))
      plt.xlabel(xlab)
      plt.ylabel(flab)

      """ Plot the integrated flux/mag """
      ax2 = plt.subplot(212,sharex=ax1)
      #plt.figure(2)
      ftot = n.cumsum(rflux)
      if zp:
         m = -2.5 * n.log10(ftot) + zp
         if logr:
            plt.semilogx(rr,m,'+')
         else:
            plt.plot(rr,m,'+')
         yl1,yl2 = plt.ylim()
         plt.ylim(yl2,yl1)
      else:
         if logr:
            plt.semilogx(rr,ftot,'+')
         else:
            plt.plot(rr,ftot,'+')
      plt.xlim(0,rmax)
      plt.title('Integrated %s centered at (%6.1f,%6.1f)'%(ttype,x0,y0))
      plt.xlabel(xlab)
      plt.ylabel(tlab)

   #-----------------------------------------------------------------------

   def set_contours(self, rms=None, hext=0, verbose=True):
      """
      Sets the contouring levels for an image.  If a subimage (i.e., cutout)
      has already been defined, then its properties are used.  Otherwise,
      the full image is used.

      The levels are set in terms of an rms, which is either passed
      explicitly via the optional rms parameter, or is determined from
      the properties of the data themselves (if rms=None).  The contours
      are multiples of (1) the rms, and (2) the contour base level (contbase), 
      which has a default value of sqrt(3).  Thus:

         clev = [-contbase**2, contbase**2, contbase**3, contbase**4,...] * rms

      Optional inputs:
         rms     - If rms is None (the default), then use the data to determine
                    the rms.  If it is not None, then use the passed value.
         hext    - Image HDU containing the data.  The default (hext=0) should
                    work for all single-extension fits files and may 
                    work for some multi-extension files.
                    NOTE: hext is ignored if the subim variable is already set
                    by, e.g., def_subim_xy or def_subim_radec
         verbose - Report contour levels if True (the default)
      """

      """ 
      Set the portion of the data to be used.  This may already have been
      set before calling set_contours.  If it has already been set, then
      self.subim will contain the data and the hext parameter will be ignored.
      If it has not been set, i.e., if self.subim is None, then set the
      data to be the full image.
      """
      if self.subim is None:
         self.subcentx = None
         self.subcenty = None
         self.get_subim(None)

      """ If no rms value has been requested, calculate the rms from the data """
      if rms is None:
         self.sigma_clip()
         rms = self.rms_clip

      """ Set the contours based on the rms and the contour base """
      maxcont = int(log((self.subim.max()/rms),self.contbase))
      if maxcont < 3:
         self.clevs = n.array([-3.,3.,self.contbase**3])
      else:
         poslevs = n.logspace(2.,maxcont,maxcont-1,base=self.contbase)
         self.clevs = n.concatenate(([-self.contbase**2],poslevs))
                                     
      if verbose:
         print "Contour levels: %f *" % rms
         print self.clevs
      self.clevs *= rms

   #-----------------------------------------------------------------------

   def read_overlay_image(self, file2name):
      """
      Reads in a second fits image in order to do an contour overlay
      on the main image.  This process is separated from the overlay_contour
      method so that a possibly large second image only needs to be read
      in once.

      Inputs:
         file2name - name of second fits file
      """

      """ Read in the image """
      try:
         self.overlay_hdu = open_fits(file2name)
      except:
         print ''
         print 'ERROR: Could not read in fits file %s' % file2name
         print ''
         return

      """ Set the overlay_im parameter """
      self.overlay_im = file2name

   #-----------------------------------------------------------------------

   def add_rough_wcs(self, ra, dec, pixscale, hext=0):
      """
      Uses the given RA, Dec, and pixel scale to create an initial guess
      of the WCS for the data.  This is a simplistic guess, with the RA
      and Dec assigned to the central pixel and the orientation assumed
      to be exactly north up, east left.

      Inputs:
         ra       - RA of the image center
         dec      - Dec of the image center
         pixscale - pixel scale in arcsec/pix
         hext     - HDU extension to assign the WCS to.  The default, 0, is
                    good for most imaging
      """

      """ First make sure that the fits file can be updated """
      if self.fitsmode != 'update':
         print ''
         print 'ERROR: Cannot update WCS unless image is loaded with'
         print '   mode="update".  The mode for this image is %s' % self.fitsmode
         print ''

      else:
         """ Create the header """
         #newwcs = wcs.ma
         print 'Not functioning yet'

   #-----------------------------------------------------------------------

   def get_subim_bounds(self, subimsize, hext=0):
      """
      Defines a subimage based on a subimage size and center
      """

      """ Get the full size of the image """
      hdr = self.hdu[hext].header
      nx = hdr['naxis1']
      ny = hdr['naxis2']

      """ 
      If the subimage center is not already set, then define it as the
      center of the full data set
      """
      if self.subcentx is None:
         self.subcentx = int((nx+1.)/2.)
      if self.subcenty is None:
         self.subcenty = int((ny+1.)/2.)

      """ 
      Define limits of subimage 
      For now does not deal with regions partially outside the input file
      """
      if subimsize is not None:
         subxy = n.atleast_1d(subimsize) # Converts subimsize to a numpy array
         subx = subxy[0]
         if subxy.size>1:
            suby = subxy[1]
         else:
            suby = subxy[0]
         halfx = int(subx/2.0)
         halfy = int(suby/2.0)
         self.subx1 = self.subcentx - halfx
         self.subx2 = self.subx1 + subx
         self.suby1 = self.subcenty - halfy
         self.suby2 = self.suby1 + suby
         self.subsizex = int(subx)
         self.subsizey = int(suby)
      else:
         self.subx1 = 0
         self.subx2 = nx
         self.suby1 = 0
         self.suby2 = ny
         self.subsizex = nx
         self.subsizey = ny

   #-----------------------------------------------------------------------

   def def_subim_xy(self, hext=0, verbose=True):
      """

      Selects the data in the subimage defined by the bounds x1, x2, y1, y2.
      These bounds are all contained within the Image class itself, and
      were either set directly (e.g., by a call to imcopy) or by the
      get_subim_bounds function (which takes a subimage center and size)

      Inputs:
         hext    - Image HDU number that contains the full image
         verbose - Print out useful information if True (the default)
      """

      """ 
      Cut out the subimage based on the bounds.
      Note that radio images often have 4 dimensions (x,y,freq,stokes)
       so for those just take the x and y data
      """
      hdr = self.hdu[hext].header
      if  hdr['naxis'] == 4:
         self.subim = self.hdu[hext].data[0,0,self.suby1:self.suby2,
                                          self.subx1:self.subx2].copy()
      else:
         self.subim = self.hdu[hext].data[self.suby1:self.suby2,
                                          self.subx1:self.subx2].copy()
      self.subim[~n.isfinite(self.subim)] = 0.
      self.subimhdr = self.hdu[hext].header.copy()

      """ Print out useful information """
      if verbose:
         print ''
         print "Cutout image center (x,y): (%d, %d)" % \
             (self.subcentx,self.subcenty)
         print "Cutout image size (x y): %dx%d" % \
             (self.subsizex,self.subsizey)

      """ 
      Update the header info, including updating the CRPIXn values if they
      are present.
      """
      self.subimhdr['ORIG_IM'] = 'Copied from %s' % self.infile
      self.subimhdr['ORIG_REG'] = 'Region in original image: [%d:%d,%d:%d]' % \
                    (self.subx1,self.subx2,self.suby1,self.suby2)

      """ Update the headers to reflect the cutout center"""
      try:
         self.subimhdr['CRPIX1'] -= self.subx1
      except:
         print 'Warning: CRPIX1 header not found in %s' % self.infile
         pass
      try:
         self.subimhdr['CRPIX2'] -= self.suby1
      except:
         print 'Warning: CRPIX2 header not found in %s' % self.infile
         pass

   #-----------------------------------------------------------------------

   def def_subim_radec(self, ra, dec, xsize, ysize=None, outscale=None, 
                       docdmatx=True, hext=0, dext=0, verbose=True):
      """
      Selects the data in the subimage defined by ra, dec, xsize, and ysize.

      The vast majority of the code is Matt Auger's (his image_cutout in 
      imagelib.py).
      Some modifications have been made by Chris Fassnacht.

      Inputs:
         ra       - Central right ascension in decimal degrees
         dec      - Central declination in decimal degrees
         xsize    - Output image x size in arcsec
         ysize    - Output image y size in arcsec
                    If ysize is None (the default) then use the same size for
                    y as is being used for x (i.e., ysize=xsize)
         outscale - Output image pixel scale, in arcsec/pix.  
                    If outscale is None (the default) then the output image
                    scale will be the same as the input image scale
         docdmatx - If set to True (the default), then put the output image scale
	            in terms of a CD matrix.  If False, then use the
		    CDELT and PC matrix formalism instead.
         hext     - Input file HDU number that contains the WCS info (default 0)
         dext     - Input file HDU number that contains the image data (default 
                    0)
         verbose - Print out informational statements if True (default=True)
      """

      """ First check to make sure that a subimage is even requested """
      if ra is None or dec is None or xsize is None:
         self.subim    = self.hdu[hext].data.copy()
         self.subimhdr = self.hdu[hext].header.copy()
         self.subsizex = self.hdu[hext].data.shape[1]
         self.subsizey = self.hdu[hext].data.shape[0]
         return

      """ Convert ra and dec into astropy.coordinates SkyCoord format """
      self.radec_to_skycoord(ra,dec)

      """ Calculate the (x,y) that is associated with the requested center"""
      self.subimhdr = self.hdu[hext].header.copy()
      x,y = wcsmwa.sky2pix(self.subimhdr,self.radec.ra.degree,
                           self.radec.dec.degree)

      """ 
      Get rough image size in pixels for the segment of input image, since the 
      pixel scale for the output image does not necessarily match that of the 
      input image.
      """
      if ysize is None:
         ysize = xsize
      inpixxsize = int(xsize / self.pixscale)
      inpixysize = int(ysize / self.pixscale)
      if outscale is None:
         outscale = self.pixscale
      self.subsizex = int(xsize / outscale)
      self.subsizey = int(ysize / outscale)

      """ Summarize the request """
      if verbose:
         print ''
         rastr = '%02d %02d %06.3f' % \
             (self.radec.ra.hms.h,self.radec.ra.hms.m,self.radec.ra.hms.s)
         decstr = '%+03d %02d %05.2f' % \
             (self.radec.dec.hms.h,self.radec.dec.hms.m,self.radec.dec.hms.s)
         print " Requested center (RA,Dec): %11.7f   %+10.6f" % \
             (self.radec.ra.deg,self.radec.dec.deg)
         print " Requested center (RA,Dec):  %s %s" % (rastr,decstr)
         print " Requested center (x,y):    %8.2f %8.2f" % (x,y)
         print " Requested image size (arcsec): %6.2f %6.2f" % \
             (xsize,ysize)
         print " Requested size in input pixels: %d %d" % (inpixxsize,inpixysize)

      """
      In order to account for rotations, etc., when cutting out the
      desired image section, start with a region that is larger
      (by a factor of 2, if the image is large enough).
      """
      x0 = max(0,int(x-inpixxsize))
      x1 = min(self.subimhdr['naxis1'],int(x+inpixxsize))
      y0 = max(0,int(y-inpixysize))
      y1 = min(self.subimhdr['naxis2'],int(y+inpixysize))
      if verbose:
         print " Cutting out image with x=%d--%d, y=%d--%d" % (x0,x1,y0,y1)

      """ Actually get the data in the large region """
      if  self.subimhdr['naxis'] == 4:
         data = self.hdu[dext].data[0,0,y0:y1,x0:x1].copy()
      else:
         data = self.hdu[dext].data[y0:y1,x0:x1].copy()
      data[~n.isfinite(data)] = 0.

      """ Update the headers to reflect the cutout center"""
      try:
         self.subimhdr['CRPIX1'] -= x0
      except:
         print 'Warning: CRPIX1 header not found in %s' % self.infile
         pass
      try:
         self.subimhdr['CRPIX2'] -= y0
      except:
         print 'Warning: CRPIX2 header not found in %s' % self.infile
         pass

      """ 
      Set up the output header and do the coordinate transform preparation 
      """
      outheader = wcsmwa.make_header(self.radec.ra.degree,self.radec.dec.degree,
                                     self.subsizex,self.subsizey,outscale,
                                     docdmatx=docdmatx)
      coords = n.indices((self.subsizey,self.subsizex)).astype(n.float32)
      skycoords = wcsmwa.pix2sky(outheader,coords[1],coords[0])
      ccdcoords = wcsmwa.sky2pix(self.subimhdr,skycoords[0],skycoords[1])
      coords[1] = ccdcoords[0]
      coords[0] = ccdcoords[1]
      self.coords = coords.copy()

      """ Transform the coordinates """
      self.subim = ndimage.map_coordinates(data,coords,output=n.float64,order=5)
      self.subimhdr = outheader.copy()

      """ Clean up """
      del data,outheader,coords,skycoords,ccdcoords

   #-----------------------------------------------------------------------

   def poststamp_xy(self, centx, centy, imsize, outfile, hext=0):
      """
      Creates a new fits file that is a cutout of the original image.  For
      this method, the image center is defined by its (x,y) coordinate
      rather than (ra,dec).

      Inputs:
         centx   - x coordinate of cutout center
         centy   - y coordinate of cutout center
         imsize  - size of cutout (postage stamp) image, in pixels
                   imsize can take any of the following formats:
                     1. A single number (which will produce a square image)
                     2. A 2-element numpy array
                     3. A 2-element list:  [xsize,ysize] 
                     4. A 2-element tuple: (xsize,ysize)
         outfile - name of output file
         hext    - HDU containing the image data in the input image (default=0)
      """

      print ""
      print "Input file:  %s" % self.infile
      print "Output file: %s" % outfile

      """ Read in relevant data """
      self.subcentx = centx
      self.subcenty = centy

      """ Make the cutout """
      self.get_subim_bounds(imsize,hext)
      self.def_subim_xy(hext)

      """ Write to the output file and clean up"""
      pf.PrimaryHDU(self.subim,self.subimhdr).writeto(outfile,clobber=True)
      print "Wrote postage stamp cutout to %s" % outfile

   #-----------------------------------------------------------------------

   def poststamp_radec(self, ra, dec, xsize, ysize, scale, outfile, 
                       docdmatx=True, hext=0, dext=0, verbose=True):
      """
      Given a central coordinate (RA,Dec), a size in pixels, and a pixel scale,
      creates an output cutout image

      The majority of the code is Matt Auger's (his image_cutout in 
      imagelib.py).
      Some modifications have been made by Chris Fassnacht.

      Inputs:
         ra       - Central right ascension in decimal degrees
         dec      - Central declination in decimal degrees
         xsize    - Output image x size in arcsec
         ysize    - Output image y size in arcsec
         scale    - Output image pixel scale, in arcsec/pix
         outfile  - Output file name
         docdmatx - If set to True (the default), then put the output image scale
	            in terms of a CD matrix.  If False, then use the
		    CDELT and PC matrix formalism instead.
         hext     - Input file HDU number that contains the WCS info (default 0)
         dext     - Input file HDU number that contains the image data 
                    (default 0)
      """

      """ Create the postage stamp data """
      self.def_subim_radec(ra,dec,xsize,ysize,scale,docdmatx,hext,dext,verbose)

      """ 
      Put the new WCS information into the original header, along with some
      additional info.
      """
      newhdr = self.hdu[hext].header.copy()
      wcskeys = \
          ['ra','dec','ctype1','ctype2','crval1','crpix1','crval2','crpix2']
      if docdmatx:
         for i in ('cd1_1','cd1_2','cd2_1','cd2_2'):
            wcskeys.append(i)
      else:
         for i in ('cdelt1', 'cdelt2', 'pc1_1','pc1_2','pc2_1','pd2_2'):
            wcskeys.append(i)

      for i in wcskeys:
         newhdr.update(i,self.subimhdr[i])
      newhdr.update('ORIG_IM',self.infile)

      """ Write the postage stamp to the output file """
      pf.PrimaryHDU(self.subim,newhdr).writeto(outfile,clobber=True)
      print "Wrote postage stamp cutout to %s" % outfile

   #-----------------------------------------------------------------------

   def imcopy(self, x1, x2, y1, y2, outfile, hext=0):
      """ 
      Description: Given the x and y coordinates of 
      the lower left corner and the upper right corner, creates a new 
      fits file that is a cutout of the original image.

      Inputs:
        x1:      x coordinate of the lower left corner of desired region
        x2:      x coordinate of the upper right corner of desired region
        y1:      y coordinate of the lower left corner of desired region
        y2:      y coordinate of the upper right corner of desired region
        outfile: file name of output image
        hext:    HDU containing the image data in the input image (default=0)
      """

      """ Get info about input image """
      inhdr = self.hdu[hext].header.copy()
      xmax = inhdr["NAXIS1"]
      ymax = inhdr["NAXIS2"]
      print ""
      print "imcopy: Input image %s has dimensions %d x %d" % \
          (self.infile,xmax,ymax)

      """Check to make sure that requested corners are inside the image"""

      """ Make sure that everything is in integer format """
      x1 = int(x1)
      y1 = int(y1)
      x2 = int(x2)
      y2 = int(y2)

      """ 
      Cut the file, and then update the CRPIXn header cards if they're there 
      """
      print "imcopy: Cutting out region between (%d,%d) and (%d,%d)" % \
          (x1,y1,x2,y2)
      outdat = self.hdu[hext].data[y1:y2,x1:x2].copy()
      inhdr.update('ORIG_IM','Copied from %s with region[%d:%d,%d:%d]' % \
                      (self.infile,x1,x2,y1,y2))
      print ""
      print "Updating CRPIXn header cards if they exist"
      print "------------------------------------------"
      try:
         crpix1 = inhdr['crpix1']
      except:
         print "   No CRPIX1 header found"
         crpix1 = n.nan
      try:
         crpix2 = inhdr['crpix2']
      except:
         print "   No CRPIX2 header found"
         crpix2 = n.nan
      if n.isnan(crpix1)==False:
         inhdr['crpix1'] -= x1
         print "   Updating CRPIX1:  %8.2f --> %8.2f" % (crpix1,inhdr['crpix1'])
      if n.isnan(crpix2)==False:
         inhdr['crpix2'] -= y1
         print "   Updating CRPIX2:  %8.2f --> %8.2f" % (crpix2,inhdr['crpix2'])
      

      """ Write to output file and clean up """
      outhdu = pf.PrimaryHDU(data=outdat,header=inhdr)
      outhdu.verify('fix')
      print "imcopy: Writing to output file %s" % outfile
      outhdu.writeto(outfile,clobber=True)
      del outdat

   #-----------------------------------------------------------------------

   def mark_fov(self, ra, dec, size, pa=0.0, color='g', lw=1):
      """
      Draws a rectangle on the currently displayed image data.  This rectangle
      can represent the FOV of a camera or, for example, the slit for a
      spectrograph.

      Required inputs:
        ra    - RA of the center of the rectangle
        dec   - Dec of the center of the rectangle
        size  - size of the rectangle IN ARCSEC.  This can be in one of the
                 following formats:
                  1. A single number (which will produce a square image)
                  2. A 2-element numpy array
                  3. A 2-element list:  [xsize,ysize] 
                  4. A 2-element tuple: (xsize,ysize)
                NOTE: The convention for a slit is for the narrow dimension
                 (i.e., the slit width) to be given as the xsize and the
                 long dimension (the slit length) to be given as the ysize

      Optional inputs: 
        pa    - Position angle of the FOV, in units of degrees E of N
                Default value of 0.0 will produce a vertical slit
        color - Line color for drawing the rectangle.  Default='g'
        lw    - Line width for drawing the rectangle.  Default=1
      """

      """ 
      This function is meaningless if the input image does not have WCS
      information in it.  Check on this before proceeding
      """

      if self.found_wcs==False:
         print ''
         print 'ERROR: Requested a FOV plot, but input image (%s)' % self.infile
         print ' does not have WCS information in it.'
         print ''
         exit()

      """ Set the rectangle size """
      imsize = n.atleast_1d(size) # Converts size to a numpy array
      xsize  = imsize[0]
      if imsize.size>1:
         ysize = imsize[1]
      else:
         ysize = xsize

      """ Set the original vertices of the FOV marker, in terms of dx and dy """
      dw = 1. * xsize / 2.
      dh = 1. * ysize / 2.
      dx0 = n.array([dw,dw,-dw,-dw,dw])
      dy0 = n.array([-dh,dh,dh,-dh,-dh])

      """ 
      Rotate the vertices.
      NOTE: With the standard convention of RA increasing to the left (i.e.,
       north up, east left) and the PA defined as north through east, we
       have to set the user's PA to its negative to get what the user wanted
       in this astronomical convention.
      """
      parad = -1. * pa * pi / 180.
      cpa = mcos(parad)
      spa = msin(parad)
      dx = dx0 * cpa - dy0 * spa
      dy = dx0 * spa + dy0 * cpa

      """ 
      Find the center point of the FOV.
      For now assume that it is close enough to the center point of the
       image that we can use the small-angle approximation to calculate
       the offsets.
      Note that we have to include the zeropos offset to get the alignment
       to be correct, since the origin of the axes may not be at the center
       pixel.
      """
      cosdec = mcos(self.radec.dec.radian)
      fovx0 = 3600. * cosdec * (ra - self.radec.ra.deg) - self.zeropos[0]
      fovy0 = 3600. * (dec - self.radec.dec.deg) - self.zeropos[1]
      fovx = fovx0 + dx
      fovy = fovy0 + dy

      """ Plot the FOV """
      xlim = plt.xlim()
      ylim = plt.ylim()
      plt.plot(fovx,fovy,color=color,lw=lw)
      plt.xlim(xlim)
      plt.ylim(ylim)

   #-----------------------------------------------------------------------

   def set_display_limits(self, fmin=-1., fmax=10., funits='sigma',
                          verbose=False):
      """

      The method used to set the flux limits for the image display.  The
       two numbers that are generated by this method will be used for the
       vmin and vmax values when the actual call to imshow (from
       matplotlib.pyplot) is made.  The two values will be stored within the
       Image class as fmin and fmax.

      Inputs:
        fmin      - Value that is used to set the minimum of the displayed flux 
                     range, where the actual value depends on the
                     value of the funits paramters (see below).
                    NOTE: If fmin is None then switch to interactive mode
        fmax      - Value that is used to set the maximum of the displayed flux 
                     range, where the actual value depends on the
                     value of the funits paramters (see below).
                    NOTE: If fmin is None then switch to interactive mode
        funits    - Either 'sigma' (the default) or 'abs'. Used to determine
                     the method of setting fmin and fmax.
                    If funits is 'abs' then the two numbers in the disprange
                     list just get stored as fmin and fmax.
                    If funits is 'sigma' (the default) then the two numbers 
                     in disprange represent the numbers of clipped standard 
                     devitations relative to the clipped mean.  In that case, 
                     the method will first calculate the clipped mean and 
                     standarddeviations and then multiply them by the passed 
                     values.

      """

      """ 
      If funits is 'abs', then just set self.fmin and self.fmax directly from
       the disprange values if those are set. Otherwise, query the user for the 
       values.
      """
      if funits == 'abs':

         """ If disprange was set, then just transfer the values """
         if fmin is not None and fmax is not None:
            self.fmin = fmin
            self.fmax = fmax

         else: # Otherwise, query the user
            """ 
            Set some default values if there aren't already some in the fmin
             and fmax containers
            """
            if self.fmin is None or self.fmax is None:
               if self.found_rms == False:
                  self.sigma_clip(verbose=verbose)
                  self.found_rms = True
               self.fmin = self.mean_clip - 1.*self.rms_clip
               self.fmax = self.mean_clip + 10.*self.rms_clip
            """ Query the user for new values """
            tmpmin = self.fmin
            tmpmax = self.fmax
            tmp = raw_input('Enter minimum flux value for display [%f]: ' \
                               % tmpmin)
            if len(tmp)>0:
               self.fmin = float(tmp)
            tmp = raw_input('Enter maximum flux value for display [%f]: ' \
                               % tmpmax)
            if len(tmp)>0:
               self.fmax = float(tmp)
         print 'fmin:  %f' % self.fmin
         print 'fmax:  %f' % self.fmax

      else:
         """
         If funits is not 'abs', then it must be 'sigma', which is the only other
         possibility, and the default value for funits.  In that case, set
         the display limits in terms of the clipped mean and sigma
         """
         
         """ Start by calculating the clipped statistics if needed """
         if self.found_rms == False:
            print "Calculating display limits"
            print "--------------------------"
            self.sigma_clip(verbose=verbose)
            self.found_rms = True

         """ If disprange is not set, then query the user for the range """
         if fmin is None or fmax is None:
            fmin = -1.
            fmax = 10.
            tmp = raw_input(
               'Enter min flux for display in terms of sigma from mean [%f]: ' \
                  % fmin)
            if len(tmp)>0:
               fmin = float(tmp)
            tmp = raw_input(
               'Enter max flux for display in terms of sigma from mean [%f]: ' \
                  % fmax)
            if len(tmp)>0:
               fmax = float(tmp)
               
         """ Set fmin and fmax in terms of clipped mean and sigma"""
         self.fmin = self.mean_clip + fmin*self.rms_clip
         self.fmax = self.mean_clip + fmax*self.rms_clip
         print " Clipped mean: %f" % self.mean_clip
         print " Clipped rms:  %f" % self.rms_clip
         s1='-' if fmin<0. else '+'
         s2='-' if fmax<0. else '+'
         print " fmin (mean %s %3d sigma):  %f" % (s1,fabs(fmin),self.fmin)
         print " fmax (mean %s %3d sigma):  %f" % (s2,fabs(fmax),self.fmax)

   #-----------------------------------------------------------------------

   def set_cmap(self, cmap='gaia'):
      """
      
      Sets the color map for the image display.

      Inputs:
       cmap - name of the color map to use.  There are only a limited
               number of choices:
               ---
               None  
               'gaia' (default)
               'gray' or 'grey'
               'gray_inv' or 'grey_inv'
               'heat' or 'hot'
               'jet'
      """

      if cmap == 'gray' or cmap == 'grey':
         self.cmap = plt.cm.gray
      elif cmap == 'gray_inv' or cmap == 'grey_inv':
         self.cmap = plt.cm.gray_r
      elif cmap == 'heat' or cmap == 'hot':
         self.cmap = plt.cm.hot
      elif cmap == 'Yl_Or_Br' or cmap == 'gaia':
         self.cmap = plt.cm.YlOrBr_r
      elif cmap == 'jet':
         self.cmap = plt.cm.jet
      else:
         print ' WARNING - Requested unknown color map.  Using gaia colors'
         self.cmap = plt.cm.YlOrBr_r


   #-----------------------------------------------------------------------

   def display_setup(self, hext=0, cmap='gaia', fmin=-1., fmax=10.,
                     funits='sigma', statsize=2048, title=None, 
                     subimdef='xy', subimcent=None, subimsize=None, 
                     dispunits='pixels', zeropos=None, axlabel=True, 
                     mask = None, show_xyproj=False, verbose=False):
      """
      Sets parameters within the Image class that will be used to actually
       display the image or the requested part of it.
      NOTE: This method is usually called from the display method, and is
       not meant to be used in a stand-alone manner
      For more information about the parameters, etc., please see the
       help information for the display method.
      """

      """ Set the region of the image to be displayed """
      if subimdef == 'radec':
         """ 
         If requesting a (RA,Dec) cutout, make the display units arcsec 
         """
         self.dispunits = 'radec'

         """ Set the display center"""
         if subimcent == None:
            ra = None
            dec = None
         else:
            ra = subimcent[0]
            dec = subimcent[1]

         """ Set the display size """
         if subimsize == None:
            xsize = None
            ysize = None
         else:
            xsize = subimsize[0]
            ysize = subimsize[1]
         self.def_subim_radec(ra,dec,xsize,ysize,hext=hext)

      else:
         """
         If not requesting a (RA,Dec) cutout, the code is simpler
         """
         self.dispunits = dispunits
         if subimcent is None:
            self.subcentx = None
            self.subcenty = None
         else:
            self.subcentx = int(subimcent[0])
            self.subcenty = int(subimcent[1])
         self.get_subim_bounds(subimsize,hext)
         self.def_subim_xy(hext)
         print "Display image center (x,y): (%d, %d)" % \
             (self.subcentx,self.subcenty)
      print "Displayed image size (x y): %dx%d" % \
          (self.subsizex,self.subsizey)
      print ''

      """ Set the image display limits """
      self.set_display_limits(fmin,fmax,funits)

      """ Set the color map """
      self.set_cmap(cmap)

      """ Set the displayed axes to be in WCS offsets, if requested """
      if self.dispunits == 'radec':
         if not self.found_wcs:
            print ''
            print "WARNING: dispunits='radec' but no WCS info in image header"
            print 'Using pixels instead'
            print ''
            self.dispunits = 'pixels'
            self.extval = None
         else:
            self.set_wcsextent(hext,zeropos)
      else:
         self.extval = None

      """ Set other display parameters """
      self.title = title

   #-----------------------------------------------------------------------

   def display_implot(self, show_xyproj=False, axlabel=True):
      """

      NOTE: DO NOT USE this routine/method unless you know exactly what
      you are doing.  It is meant to be called from the display() routine/method,
      as well as in a few other specialized cases, and is NOT meant to 
      have stand-alone functionality.

      Please see the help for the display method for more information.
      """

      """ 
      Set up for displaying the image data
       - If show_xyproj is False (the default), then just show self.subim
       - If show_xyproj is True, then make a three panel plot, with
          Panel 1: self.subim (i.e., what you would see in the default behavior)
          Panel 2: Projection of data in self.subim onto the x-axis
          Panel 3: Projection of data in self.subim onto the x-axis
        - Setting show_xyproj=True is most useful when evaluating, e.g., a 
          star in the image data.  The projections along the two axes of the
          cutout can be useful for evaluating whether the object is a star and/or
          whether it is saturated
      """
      if show_xyproj:
         self.fig2 = plt.figure(figsize=(10,3))
         self.fig2.add_subplot(131)
      else:
         self.fig1 = plt.gcf()
         self.ax1 = plt.gca()

      """ Choose the scaling for the display """
      fdiff = fabs(self.fmax - self.fmin)
      if self.fscale == 'log':
         data = self.subim - self.subim.min() + 1.
         data = n.log10(data)
         vmin = log10(self.fmin - self.subim.min() + 1.)
         vmax = log10(self.fmax - self.subim.min() + 1.)
      else:
         """ Linear scaling is the default """
         data = self.subim
         vmin = self.fmin
         vmax = self.fmax
         

      """ Display the image data """
      plt.imshow(data,origin='bottom',cmap=self.cmap,vmin=vmin,
                 vmax=vmax,interpolation='nearest',extent=self.extval)
      if axlabel is True:
         if self.dispunits == 'radec':
            plt.xlabel(r"$\Delta \alpha$ (arcsec)")
            plt.ylabel(r"$\Delta \delta$ (arcsec)")
         else:
            plt.xlabel('x (pix)')
            plt.ylabel('y (pix)')
      if self.title is not None:
         plt.title(title)

      """ 
      Now add the x and y projections if requested (i.e., if show_xyproj is True
      """
      if show_xyproj:
         self.fig2.add_subplot(132)
         xsum = self.subim.sum(axis=0)
         plt.plot(xsum)
         plt.xlabel('Relative x Coord')
         self.fig2.add_subplot(133)
         ysum = self.subim.sum(axis=1)
         plt.plot(ysum)
         plt.xlabel('Relative y Coord')
         self.cid_keypress2 = self.fig2.canvas.mpl_connect('key_press_event',
                                                           self.keypress)
         self.fig2.show()

   #-----------------------------------------------------------------------

   def display(self, hext=0, cmap='gaia', fmin=-1., fmax=10., funits='sigma',
               statsize=2048, title=None, 
               subimdef='xy', subimcent=None, subimsize=None, 
               dispunits='pixels', zeropos=None, axlabel=True, 
               mask = None, show_xyproj=False, verbose=False):
      """
      The main way to display the image data contained in the Image class.
      The default is to display the entire image, but it is possible to display
      cutouts (subimages), which can be defined either by (RA,Dec) or (x,y)

      Optional inputs:
         subimsize - size of the subimage to be displayed, either in pixels
                      (the default) or arcsec (if subimdef='radec').  The
                      default, designated by subimsize=None, is to display
                      the entire image.  The subimsize parameter can take
                      any of the following formats:
                         1. A single number (which will produce a square image)
                         2. A 2-element numpy array
                         3. A 2-element list:  [xsize,ysize] 
                         4. A 2-element tuple: (xsize,ysize)
         zeropos   - NOTE: Only used if dispunits='radec'
                      By default, which happens when zeropos=None, the (0,0)
                      point on the output image, as designated by the image
                      axis labels, will be at the center of the image.  However,
                      you can shift the (0,0) point to be somewhere else by
                      setting zeropos.  For example, zeropos=(0.5,0.3) will
                      shift the origin to the point that would have been
                      (0.5,0.3) if the origin were at the center of the image

      """
      print ""
      print "Input file:  %s" % self.infile

      """ Set up the parameters that will be needed to display the image """
      self.display_setup(hext=hext,cmap=cmap,fmin=fmin,fmax=fmax,funits=funits,
                         statsize=statsize,title=title,subimdef=subimdef,
                         subimcent=subimcent,subimsize=subimsize, 
                         dispunits=dispunits,zeropos=zeropos,axlabel=axlabel,
                         mask=mask,show_xyproj=show_xyproj,verbose=verbose)


      """ Now display the data """
      self.display_implot(show_xyproj, axlabel)

      #del data

   #-----------------------------------------------------------------------

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

def make_cutout(infile, ra, dec, imsize, scale, outfile, whtsuff=None,
                makerms=False, rmssuff='_rms', hext=0, dext=0, verbose=True):
   """
   Makes a cutout from an input image, based on a requested (RA,Dec) center
   and an image size in arcsec.
   Additional, optional functionality:
      - Makes the same-sized cutout for the associated weight file.  Done if
        whtsuff is not None
      - Makes an RMS image following Matt Auger's recipe for the NIRC2 data
        This only happens if BOTH whtsuff is not None AND makerms is True

   Inputs:
      infile  - input file
      ra      - RA of cutout center, in decimal degrees
      dec     - Dec of cutout center, in decimal degrees
      imsize  - output image size, in arcsec
      scale   - pixel scale of output image, in arcsec/pix
      outfile - output file name
      whtsuff - suffix for input weight file, if a cutout of the weight file
                is also desired.  If whtsuff is None (the default) then no
                weight-file cutout is made.
                Example: whtsuff='_wht' means that for infile='foo.fits' the
                  weight file is called 'foo_wht.fits'
      makerms - Set to True to make, in addition, an output rms file following
                Matt Auger's prescription for the NIRC2 data.  Default is False
                NB: Both whtsuff being something other than None and 
                 makerms=True are required for an output rms file to be created.
      rmssuff - Suffix for output rms file.  Default='_rms' means that for
                infile='foo.fits', the output file will be 'foo_rms.fits'
      hext    - Input file HDU number that contains the WCS info (default 0)
      dext    - Input file HDU number that contains the image data (default 0)
   """

   """ Make the input file cutout """
   infits = Image(infile)
   infits.poststamp_radec(ra,dec,imsize,imsize,scale,outfile,hext=hext,
                          dext=dext,verbose=verbose)

   """ Make the weight file cutout, if requested """
   if whtsuff is not None:
      whtfile = infile.replace('.fits','%s.fits' % whtsuff)
      outwht  = outfile.replace('.fits','%s.fits' % whtsuff)
      whtfits = Image(whtfile)
      whtfits.poststamp_radec(ra,dec,imsize,imsize,scale,outwht,hext=hext,
                              dext=dext,verbose=verbose)

   """ Make output RMS file, if requested """
   # CODE STILL TO COME

   """ Clean up """
   infits.close()
   if whtsuff is not None:
      whtfits.close()
   

#-----------------------------------------------------------------------

def open_fits(infile, mode='copyonwrite'):
   """
   Opens a fits file, allowing for the possibility of the missing end that
   plagues some of the NIR instruments on Keck.

   Inputs:
      infile - input file name
      mode   - [OPTIONAL] mode of opening the file.  Note that the default
               value ('copyonwrite') is the pyfits default value.  Look at
               the help information for pyfits open for other options.
   """

   try:
      hdulist = pf.open(infile,mode=mode)
   except:
      try:
         """ Try to get rid of read-in warnings """
         import warnings
         warnings.filterwarnings('ignore')
         hdulist = pf.open(infile,mode=mode,ignore_missing_end=True)
      except:
         print ""
         print "ERROR. Could not open fits file %s" % infile
         return None

   return hdulist

#---------------------------------------------------------------------------

def imcopy(infile,x1,x2,y1,y2,outfile):
   """ 
   Description: Given a fits file name, the x and y coordinates of 
     the lower left corner and the upper right corner, creates a new 
     fits file that is a cutout of the original image.

   Inputs:
     infile:  file name of input image
     x1:      x coordinate of the lower left corner of desired region
     x2:      x coordinate of the upper right corner of desired region
     y1:      y coordinate of the lower left corner of desired region
     y2:      y coordinate of the upper right corner of desired region
     outfile: file name of output image
   """

   hdu_list = open_fits(infile)
   if hdu_list is None:
      return

   """ Get info about input image """
   inhdr = hdu_list[0].header
   xmax = inhdr["NAXIS1"]
   ymax = inhdr["NAXIS2"]
   print ""
   print "imcopy: Input image %s has dimensions %d x %d" % (infile,xmax,ymax)

   """Check to make sure that requested corners are inside the image"""

   """ Make sure that everything is in integer format """
   x1 = int(x1)
   y1 = int(y1)
   x2 = int(x2)
   y2 = int(y2)

   """ Cut the file, and then update the CRPIXn header cards if they're there """
   print "imcopy: Cutting out region between (%d,%d) and (%d,%d)" % \
      (x1,y1,x2,y2)
   outdat = hdu_list[0].data[y1:y2,x1:x2]
   print ""
   print "Updating CRPIXn header cards if they exist"
   print "------------------------------------------"
   try:
      crpix1 = inhdr['crpix1']
   except:
      print "   No CRPIX1 header found"
      crpix1 = n.nan
   try:
      crpix2 = inhdr['crpix2']
   except:
      print "   No CRPIX2 header found"
      crpix2 = n.nan
   if n.isnan(crpix1)==False:
      inhdr['crpix1'] -= x1
      print "   Updating CRPIX1:  %8.2f --> %8.2f" % (crpix1,inhdr['crpix1'])
   if n.isnan(crpix2)==False:
      inhdr['crpix2'] -= y1
      print "   Updating CRPIX2:  %8.2f --> %8.2f" % (crpix2,inhdr['crpix2'])
      

   """ Write to output file """
   outhdu = pf.PrimaryHDU(data=outdat,header=inhdr)
   outhdu.verify('fix')
   print "imcopy: Writing to output file %s" % outfile
   outhdu.writeto(outfile,clobber=True)

   return

#---------------------------------------------------------------------------

def poststamp(infile,centx,centy,xsize,ysize,outfile):
   """ 
   Description: Given a fits file name, an image center, and an image size,
     creates a new fits file that is cutout of part of the original image.
   """

   im = Image(infile)
   im.poststamp_xy(centx,centy,(xsize,ysize),outfile)
   im.close()
   del im

   return

#-----------------------------------------------------------------------

def del_history_cards(hdu):
   """
   Deletes the history cards from the input HDU
   """

   print "del_history_cards is not yet implemented"

#-----------------------------------------------------------------------

def image_cutout_hdu(hdu,ra,dec,xsize,ysize,scale,hext=0,dext=0,verbose=True):
   """
   Given a central coordinate (RA,Dec), a size in pixels, and a pixel scale,
   creates an output cutout (as a HDU).

   The vast majority of the code is Matt Auger's (his image_cutout in 
   imagelib.py).
   Some modifications have been made by Chris Fassnacht.

   """

   """ Convert ra and dec to decimal degrees if necessary """
   if wcsmwa.is_degree(ra)==False:
      ra = wcsmwa.ra2deg(ra)
   if wcsmwa.is_degree(dec)==False:
      dec = wcsmwa.dec2deg(dec)

   """ Calculate the (x,y) that is associated with the requested center"""
   inhdr = hdu[hext].header.copy()
   x,y = wcsmwa.sky2pix(inhdr,ra,dec)

   """ 
   Get rough image size in pixels for the segment of input image, since the pixel
   scale for the output image does not necessarily match that of the input 
   image.
   Note that wcsinfo[2] is the CD matrix.  The rough scale is just the
   x-axis scale, assuming that the y-axis scale is the same.
   """
   wcsinfo = wcsmwa.parse_header(inhdr)
   inscale = sqrt(wcsinfo[2][0,0]**2 + wcsinfo[2][1,0]**2)*3600.
   wcsxsize = xsize * scale
   wcsysize = ysize * scale
   inpixxsize = int(wcsxsize / inscale)
   inpixysize = int(wcsysize / inscale)

   """ Summarize the request """
   if verbose:
      print " Requested center (RA,Dec): %11.7f %+10.6f" % (ra,dec)
      print " Requested center (x,y):    %8.2f %8.2f" % (x,y)
      print " Requested image size (arcsec): %6.2f %6.2f" % \
          (wcsxsize,wcsysize)
      print " Requested size in input pixels: %d %d" % (inpixxsize,inpixysize)

   """
   In order to account for rotations, etc., when cutting out the
   desired image section, start with a region that is larger
   (by a factor of 2, if the image is large enough).
   """
   x0 = max(0,int(x-inpixxsize))
   x1 = min(inhdr['naxis1'],int(x+inpixxsize))
   y0 = max(0,int(y-inpixysize))
   y1 = min(inhdr['naxis2'],int(y+inpixysize))
   if verbose:
      print " Cutting out image with x=%d--%d, y=%d--%d" % (x0,x1,y0,y1)

   """ Actually get the data in the large region """
   if  inhdr['naxis'] == 4:
      data = hdu[dext].data[0,0,y0:y1,x0:x1].copy()
   else:
      data = hdu[dext].data[y0:y1,x0:x1].copy()
   data[~n.isfinite(data)] = 0.

   """ Update the headers to reflect the cutout center"""
   inhdr.update('CRPIX1',inhdr['CRPIX1']-x0)
   inhdr.update('CRPIX2',inhdr['CRPIX2']-y0)

   """ Set up the output header and do the coordinate transform preparation """
   outheader = wcsmwa.make_header(ra,dec,xsize,ysize,scale)
   coords = n.indices((ysize,xsize)).astype(n.float32)
   skycoords = wcsmwa.pix2sky(outheader,coords[1],coords[0])
   ccdcoords = wcsmwa.sky2pix(inhdr,skycoords[0],skycoords[1])
   coords[1] = ccdcoords[0]
   coords[0] = ccdcoords[1]
   print coords.shape

   """ Create the output HDU, while at the same time transforming the coords """
   out = pf.PrimaryHDU(ndimage.map_coordinates(data,coords,
                                               output=n.float64,order=5))
   out.header = outheader.copy()

   return out

#-----------------------------------------------------------------------

def make_snr_image(infile, xmin=0, xmax=0, ymin=0, ymax=0, 
                   outfile=None, hext=0):
   """
   Reads in the data from the input image and estimates the RMS noise
   in the region defined by xmin, etc. (default is to use the entire image).
   Divides the image by that number and either writes the output SNR image
   to a file or returns it as a data array.
   """

   """ Read in the data """
   try:
      hdulist = open_fits(infile)
   except:
      print ""
      print "ERROR.  Could not open %s" % infile
      return

   data = hdulist[hext].data

   """ Set up the image section """
   x1,x2,y1,y2 = ccd.define_trimsec(hdulist[hext],xmin,xmax,ymin,ymax)

   """ Calculate the data in the requested region """
   datstat = data[y1:y2,x1:x2].copy()
   m,s = ccd.sigma_clip(datstat)
   print ""
   print "For image section [x1:x2,y1:y2] = [%d:%d,%d:%d], rms = %f" \
       % (x1,x2,y1,y2,s)

   """ Create SNR image and either write it out or return it """
   snrim = datstat / s
   del datstat
   del data
   del hdulist
   if outfile is None:
      print ""
      return snrim
   else:
      pf.PrimaryHDU(snrim).writeto(outfile)
      print "Wrote SNR image to output file %s" % outfile
      print ""
      del snrim
      return None

#-----------------------------------------------------------------------

def overlay_contours_hdu(hdu1, hdu2, ra, dec, imsize, pixscale, rms1=None,
                         rms2=None, fmax=10., title=None, showradec=True,
                         verbose=True):
   """
   The machinery used to do all of the work for overlay_contours.

   See overlay_contours for parameter explanations, etc.
   """

   """ 
   Create HDUs containing the postage stamp cutouts, sampled on the 
   output grid
   """
   outsize = int(imsize / pixscale)
   print "Image 1"
   print "------"
   ohdu1 = image_cutout_hdu(hdu1,ra,dec,outsize,outsize,pixscale,verbose=verbose)
   print "Image 2"
   print "------"
   ohdu2 = image_cutout_hdu(hdu2,ra,dec,outsize,outsize,pixscale,verbose=verbose)

   """ Get data from the HDUs """
   dat1 = ohdu1.data.copy()
   dat2 = ohdu2.data.copy()

   """ Set display limits """
   if rms1 is not None:
      vmin = 0.0
      vmax = fmax * rms1
   else:
      m1,s1 = ccd.sigma_clip(dat1)
      vmin = m1 - s1
      vmax = m1 + fmax * s1

   """ Invert greyscale for better display """
   dat1 *= -1.0
   tmp = -1.0 * vmin
   vmin = -1.0 * vmax
   vmax = tmp

   """ Set contour levels """
   if rms2 is None:
      m2,rms2 = ccd.sigma_clip(dat2)
   contbase = sqrt(3.)
   maxcont = int(log((dat2.max()/rms2),contbase))
   if maxcont < 3:
      clevs = n.array([-3.,3.,contbase**3])
   else:
      clevs = n.concatenate(([-contbase**2],
                             n.logspace(2.,maxcont,maxcont-1,base=contbase)))
   print "Contour levels: %f *" % rms2
   print clevs
   clevs *= rms2

   """ Actually plot the data and overlay """
   coords = n.indices((outsize,outsize)).astype(n.float32)
   maxi = outsize - 1
   #skyc   = wcsmwa.pix2sky(hdu1.header,coords[1],coords[0])
   pltc   = (coords - outsize/2.)*pixscale
   pltc[1] *= -1.
   plt.imshow(dat1,origin='bottom',vmin=vmin,vmax=vmax,cmap=plt.cm.gray,
              interpolation='none',aspect='equal',
              extent=(pltc[1][0,0],pltc[1][maxi,maxi],
                      pltc[0][0,0],pltc[0][maxi,maxi]))
   plt.xlabel(r"$\Delta \alpha$ (arcsec)")
   plt.ylabel(r"$\Delta \delta$ (arcsec)")
   
   plt.contour(pltc[1],pltc[0],dat2,clevs,colors='r')
   if title is not None:
      plt.title(title)
   if showradec:
      xloc = pltc[1][0,0] - 0.03*imsize
      yloc = pltc[0][0,0] + 0.04*imsize
      plt.text(xloc,yloc,"(%11.7f,%+10.6f)" % \
                  (ra,dec),fontsize=14,weight='bold')

   """ Clean up """
   del dat1,dat2,ohdu1,ohdu2

#---------------------------------------------------------------------------

def overlay_contours(infile1, infile2, ra, dec, imsize, pixscale=None, 
                     zeropos=None, fmax=10., hext1=0, 
                     hext2=0, rms2=None, ccolor2='r', 
                     infile3=None, hext3=0, rms3=None, ccolor3='b',
                     title=None, showradec=True,
                     verbose=True):
   """
   Creates a postage-stamp cutout (of size imgsize arcsec) of the data in the
    Image class and then overlays contours from the second image (infile2).

   Required inputs:
      infile1   - fits file containing the data for the first image
      infile2   - fits file containing the data for the second image
      ra        - single number containing RA for image center
                  (best if in decimal degrees)
      dec       - single number containing Dec for image center
                  (best if in decimal degrees)
      imsize    - length of one side of output image, in arcsec
   Optional inputs:
      pixscale  - pixel scale of output image, in arcsec/pix
                   If pixscale is None (the default) then just use the
                   native pixel scale of each of the input images.
      zeropos   - By default, which happens when zeropos=None, the (0,0)
                   point on the output image, as designated by the image
                   axis labels, will be at the center of the image.  However,
                   you can shift the (0,0) point to be somewhere else by
                   setting zeropos.  For example, zeropos=(0.5,0.3) will
                   shift the origin to the point that would have been
                   (0.5,0.3) if the origin were at the center of the image
      fmax      - upper range for plotting the greyscale in the first image,
                   expressed as the number of sigma above the clipped mean.
                   Default = 10.
      hext1     - HDU containing the actual image data for the first file
                   default=0
      hext2     - HDU containing the actual image data for the second file
                   default=0
      rms2      - user-requested rms for data in the second image. If set to 
                   None (the default) then calculate rms from the cutout data
                   themselves
      ccolor2   - color for the contours from infile2.  Default='r'
      infile3   - OPTIONAL name of a third image, to be used for a second set
                   of contours in a different line style.  Default=None
      hext3     - HDU containing the actual image data for the third file
                   default=0
      rms3      - user-requested rms for data in the optional third image. If 
                   set to  None (the default) then calculate rms from the
                   cutout data themselves
      ccolor3   - color for the contours from infile3.  Default='b' (black)
      title     - title for the figure.  The default value (None) will show
                   no title
      showradec - print the RA and Dec of the center of the image, in decimal
                   degrees, on the figure.  Default=True
      verbose   - print out useful information while running.  Default=True
   """

   """ Read the input images """
   try:
      im1 = Image(infile1)
   except:
      print ''
      print 'ERROR: Could not properly open %s' % infile1
      return
   print "   .... Done"
   try:
      im2 = Image(infile2)
   except:
      print ''
      print 'ERROR: Could not properly open %s' % infile2
      return
   print "   .... Done"

   """ 
   Make cutouts of the appropriate size for each of the input images
   For the first image this is done via a call to display
   """
   im1.display(hext=hext1,cmap='gray_inv',subimdef='radec',subimcent=(ra,dec),
               subimsize=(imsize,imsize),dispunits='radec',fmax=fmax,
               zeropos=zeropos)
   im2.def_subim_radec(ra,dec,imsize,outscale=pixscale)

   """ Set contour levels for the second image """
   im2.set_contours(rms2)

   """ 
   Change the axis labels to be offsets in arcseconds from a fiducial point
   in the image, which is set to be the origin.
   Default value for the origin is the center of the image.
   Override the default value by setting the zeropos parameter
   """
   im2.set_wcsextent(zeropos=zeropos)

   """ Plot the contours """
   plt.contour(im2.subim,im2.clevs,colors=ccolor2,extent=im2.extval)

   """ If there is a third image, plot contours from it """
   if infile3 is not None:
      try:
         im3 = Image(infile3)
      except:
         print ''
         print 'ERROR: Could not properly open %s' % infile3
         return
      im3.def_subim_radec(ra,dec,imsize,outscale=pixscale)
      im3.set_contours(rms3)
      im3.set_wcsextent(zeropos=zeropos)
      plt.contour(im3.subim,im3.clevs,colors=ccolor3,extent=im3.extval)

   """ Clean up """
   im1.close()
   im2.close()
   del im1,im2
   if infile3 is not None:
      im3.close()
      del im3

#---------------------------------------------------------------------------

def overlay_contours_old(infile1, infile2, ra, dec, imsize, pixscale, rms1=None,
                         rms2=None, fmax=10., title=None, showradec=True,
                         verbose=True):
   """
   Creates a postage-stamp cutout (of size imgsize arcsec) of the data in the
    Image class and then overlays contours from the second image (infile2).

   Inputs:
      infile1   - fits file containing the data for the first image
      infile2   - fits file containing the data for the second image
      ra        - single number containing RA for image center
                  (best if in decimal degrees)
      dec       - single number containing Dec for image center
                  (best if in decimal degrees)
      imsize    - length of one side of output image, in arcsec
      pixscale  - pixel scale of output image, in arcsec/pix
      rms1      - user-requested rms for data in the first image. If set to 
                   None (the default) then calculate rms from the cutout data
                   themselves
      rms2      - user-requested rms for data in the second image. If set to 
                   None (the default) then calculate rms from the cutout data
                   themselves
   """

   """ Read the input images """
   print ""
   print "Opening %s" % infile1
   try:
      hdu1 = open_fits(infile1)
   except:
      return
   print "   .... Done"
   print "Opening %s" % infile2
   try:
      hdu2 = open_fits(infile2)
   except:
      return
   print "   .... Done"

   """ 
   Produce the plot
   The overlay_contours_hdu code is separated from overlay_contours so that
   it can be called from other functions.
   """
   overlay_contours_hdu(hdu1,hdu2,ra,dec,imsize,pixscale,rms1,rms2,fmax,
                        title,showradec,verbose)

   """ Clean up """
   del hdu1,hdu2


#-----------------------------------------------------------------------

def calc_sky_from_seg(infile,segfile):
   """
   Description: Calculates the sky level in an image using only 
      those regions in which SExtractor's segmentation file has
      a value of 0.

   Inputs:
    infile:   input fits file
    segfile:  SExtractor segmentation file associated with infile.
   """

   """ Load data """
   indat = pf.getdata(infile)
   segdat = pf.getdata(segfile)

   """ Set mask region and select associated regions of indat """
   mask = segdat == 0  
   sky = indat[mask]  
   # NB: These preceding 2 lines could have been combined as
   #   sky = indat[segdat==0]

   """ Calculate statistics """
   print "Statistics of sky outside masked regions"
   print "----------------------------------------"
   print "  N_pix  = %d" % sky.size
   print "  Median = %f" % n.median(sky)
   print "  Mean   = %f" % n.mean(sky)

   return

#-----------------------------------------------------------------------

def display_image(infile, inhdu=0, cmap='gaia', fmin=-1.0, fmax=10.0,
                  funits='sigma'):
   """
   Displays the image data contained in an input fits file.
   Does this through a call to the Image class, which is returned.
   """

   # Read in the image
   image = Image(infile)

   # Display the image.  Note that for now this call does not include
   #  all of the possible parameters defined in the Image.display method
   #  (missing, e.g., wtfile, statsize, extent)
   image.display(hext=inhdu,cmap=cmap,fmin=fmin,fmax=fmax,funits=funits)

   return image

#-----------------------------------------------------------------------

def plot_cat(fitsfile, catfile, xcol=0, ycol=1, marksize=20., markcolor='g',
             inhdu=0, cmap='gray', fmin=-1., fmax=10., funits='sigma'):
   """
   Plots a catalog (e.g., one generated by SExtractor), on top of a fits
   image.

   Inputs:
      fitsfile  - input fits data file containing the image
      catfile   - input file containing the object catalog
      xcol      - column in the input file with the object x coordinates
                  (remember that the first column corresponds to xcol=0)
                  default value: 0
      ycol      - column in the input file with the object y coordinates
                  (remember that the second column corresponds to ycol=1)
                  default value: 1
      marksize  - size of circles marking the objects on the image, in points
                  default value: 20.0
      markcolor - color of circles marking the objects
                  default value: 'g'
      inhdu     - header-data unit containing the image data in the input fits
                  image.  The default value of 0 is appropriate for all simple
                  fits images (i.e., those without multiple extensions).
                  default value: 0
      cmap      - color map used to present the image data
                  default value: 'gray'
      fmin      - sets display range for input image
                  default value: -1.0 (1-sigma below the clipped mean)
      fmax      - sets display range for input image
                  default value: 10.0 (10-sigma below the clipped mean)
      funits    - units for display range.  Default value = 'sigma', to
                  calculate range in terms of sigma above and below the
                  clipped mean
   """

   """ Plot the image """
   try:
      display_image(fitsfile,inhdu=inhdu,cmap=cmap,fmin=fmin,fmax=fmax,
                    funits=funits)
   except:
      print ""
      print "Image display failed when called from plot_cat."
      print ""
      return
   nx = pf.getval(fitsfile,'naxis1')
   ny = pf.getval(fitsfile,'naxis2')

   """ Read in the catalog and extract the x and y coordinates """
   data = n.loadtxt(catfile)
   x = data[:,xcol]
   y = data[:,ycol]

   """ Mark the catalog objects """
   plt.plot(x,y,'o',ms=marksize,mec=markcolor,mfc="none")
   plt.xlim(0,nx-1)
   plt.ylim(0,ny-1)

#-----------------------------------------------------------------------

def read_wcsinfo(fitsfile, inhdu=0, verbose=True):
   """
   Reads the wcs information from a fits file and returns that information.

   NOT YET FUNCTIONAL
   """

   validwcs = False

   try:
      hdulist = open_fits(fitsfile)
   except:
      return
   hdr = hdulist[inhdu].header

   """ Get the CTYPE, CRPIX, CRVAL pairs first """
   ctype = []
   crpix = n.zeros(hdr['naxis'])
   crval = n.zeros(hdr['naxis'])
   for i in range(hdr['naxis']):
      typecard = 'ctype%d' % (i+1)
      pixcard = 'crpix%d' % (i+1)
      valcard = 'crval%d' % (i+1)
      ctype.append(hdr[typecard])
      crpix[i] = hdr[pixcard]
      crval[i] = hdr[valcard]

   """ Search for CD matrix """
   cdmatx = n.zeros((hdr['naxis'],hdr['naxis']))
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
      cdelt = n.zeros(hdr['naxis'])
      crota = n.zeros(hdr['naxis'])
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
         junk = 0.0

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
   

#-----------------------------------------------------------------------

def make_wcs_from_tel_pointing(infile, pixscale, rotatekey=None, 
                               rakey='ra', deckey='dec', hduext=0):
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
      pixscale  -  desired pixel scale for the output header (should be the
                   expected pixel scale for the camera)
      rotatekey -  [OPTIONAL] if the image was taken with a rotation, that can
                   be incorporated into the output header.  In that case,
                   rotatekey needs to be set to the keyword for the header
                   card that contains the rotation value (these card keywords
                   are very non-standardized).
      rakey     -  [OPTIONAL] header keyword for RA of pointing, if different
                   from the default value of 'RA'
      deckey    -  [OPTIONAL] header keyword for Dec of pointing, if different
                   from the default value of 'Dec'
      hduext    -  [OPTIONAL] HDU number to use.  Default=0 (i.e., the primary
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
      inhdr = hdu[hduext].header
   except:
      print ''
      print 'ERROR: Could not open header for HDU %d in %s' % (infile,hduext)
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

#-----------------------------------------------------------------------

def make_wcs_from_ref_tel(reffile, infile, pixscale, rotatekey=None,
                          hduext=0, rakey='ra', deckey='dec'):
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
      hduext    -  HDU number to use.  Default=0 (i.e., the primary HDU)
      rakey     -  [OPTIONAL] header keyword for RA of pointing, if different
                   from the default value of 'RA'
      deckey    -  [OPTIONAL] header keyword for Dec of pointing, if different
                   from the default value of 'Dec'
   """

   print ""

   """ Open the files """
   try:
      refhdu = open_fits(reffile)
   except:
      print "ERROR.  Could not read reference file %s" % reffile
      return
   refhdr = refhdu[hduext].header

   try:
      inhdu = open_fits(infile,mode='update')
   except:
      return
   inhdr = inhdu[hduext].header

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

