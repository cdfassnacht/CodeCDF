"""
imfuncs.py - A library of functions to do various basic image processing
             operations

NB: Some of these functions are (slowly) being incorporated into the
Image class, the code for which is at the beginning of this file.

Functions:
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

try:
   from astropy.io import fits as pf
except:
   import pyfits as pf
import numpy as n
from scipy import ndimage
import matplotlib.pyplot as plt
from math import log,sqrt
import wcs
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

      """ Initialize display parameters """
      self.found_rms = False
      self.mean_clip = 0.0
      self.rms_clip = 0.0
      self.statsize = 2048

      """ Initialize other parameters """
      self.overlay_im = None

   #-----------------------------------------------------------------------

   def close(self):
      """
      Closes the image
      """

      self.hdu.close()

   #-----------------------------------------------------------------------

   def read_overlay_image(file2name):
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
         newwcs = wcs.ma

   #-----------------------------------------------------------------------

   def get_subim_bounds(self, hdu, subimsize, subimcent):
      """
      Defines a subimage based on a subimage size and center
      """

      hdr = self.hdu[hdu].header
      nx = hdr['naxis1']
      ny = hdr['naxis1']

      """ Define subimage center """
      if subimcent is None:
         self.subcentx = int((nx+1.)/2.)
         self.subcenty = int((ny+1.)/2.)
      else:
         self.subcentx = int(subimcent[0])
         self.subcenty = int(subimcent[1])

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

   def def_subim_xy(self, hext=0):
      """

      Selects the data in the subimage defined by the bounds x1, x2, y1, y2.
      These bounds are all contained within the Image class itself, and
      were either set directly (e.g., by a call to imcopy) or by the
      get_subim_bounds function (which takes a subimage center and size)

      Inputs:
         hext    - Image HDU number that contains the full image

      """

      """ Cut out the subimage based on the bounds """
      self.subim = self.hdu[hext].data[self.suby1:self.suby2,
                                       self.subx1:self.subx2].copy()
      self.subimhdr = self.hdu[hext].header.copy()

      """ 
      Update the header info, including updating the CRPIXn values if they
      are present.
      """
      self.subimhdr.update('ORIG_IM','Copied from %s with region[%d:%d,%d:%d]'\
                              % (self.infile,self.subx1,self.subx2,self.suby1,
                                 self.suby2))

      """ Update the headers to reflect the cutout center"""
      try:
         self.subimhdr.update('CRPIX1',hdr['CRPIX1']-self.subx1)
      except:
         pass
      try:
         self.subimhdr.update('CRPIX2',hdr['CRPIX2']-self.suby1)
      except:
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
         self.subim     = self.hdu[hext].data.copy()
         self.subimhdr  = self.hdu[hext].header.copy()
         self.subsizex = self.hdu[hext].data.shape[1]
         self.subsizey = self.hdu[hext].data.shape[0]
         return

      """ Convert ra and dec to decimal degrees if necessary """
      if wcs.is_degree(ra)==False:
         ra = wcs.ra2deg(ra)
      if wcs.is_degree(dec)==False:
         dec = wcs.dec2deg(dec)

      """ Calculate the (x,y) that is associated with the requested center"""
      inhdr = self.hdu[hext].header.copy()
      x,y = wcs.sky2pix(inhdr,ra,dec)

      """ 
      Get rough image size in pixels for the segment of input image, since the 
      pixel scale for the output image does not necessarily match that of the 
      input image.
      Note that wcsinfo[2] is the CD matrix.  The rough scale is just the
      x-axis scale, assuming that the y-axis scale is the same.
      """
      wcsinfo = wcs.parse_header(inhdr)
      inscale = sqrt(wcsinfo[2][0,0]**2 + wcsinfo[2][1,0]**2)*3600.
      if ysize is None:
         ysize = xsize
      inpixxsize = int(xsize / inscale)
      inpixysize = int(ysize / inscale)
      if outscale is None:
         outscale = inscale
      self.subsizex = int(xsize / outscale)
      self.subsizey = int(ysize / outscale)

      """ Summarize the request """
      if verbose:
         print " Requested center (RA,Dec): %11.7f %+10.6f" % (ra,dec)
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
      x1 = min(inhdr['naxis1'],int(x+inpixxsize))
      y0 = max(0,int(y-inpixysize))
      y1 = min(inhdr['naxis2'],int(y+inpixysize))
      if verbose:
         print " Cutting out image with x=%d--%d, y=%d--%d" % (x0,x1,y0,y1)

      """ Actually get the data in the large region """
      if  inhdr['naxis'] == 4:
         data = self.hdu[dext].data[0,0,y0:y1,x0:x1].copy()
      else:
         data = self.hdu[dext].data[y0:y1,x0:x1].copy()
      data[~n.isfinite(data)] = 0.

      """ Update the headers to reflect the cutout center"""
      inhdr.update('CRPIX1',inhdr['CRPIX1']-x0)
      inhdr.update('CRPIX2',inhdr['CRPIX2']-y0)

      """ 
      Set up the output header and do the coordinate transform preparation 
      """
      outheader = wcs.make_header(ra,dec,self.subsizex,self.subsizey,outscale,
                                  docdmatx=docdmatx)
      coords = n.indices((self.subsizey,self.subsizex)).astype(n.float32)
      skycoords = wcs.pix2sky(outheader,coords[1],coords[0])
      ccdcoords = wcs.sky2pix(inhdr,skycoords[0],skycoords[1])
      coords[1] = ccdcoords[0]
      coords[0] = ccdcoords[1]
      print coords.shape

      """ Transform the coordinates """
      self.subim = ndimage.map_coordinates(data,coords,output=n.float64,order=5)
      self.subimhdr = outheader.copy()

      """ Clean up """
      del data,outheader,coords,skycoords,ccdcoords

   #-----------------------------------------------------------------------

   def poststamp_xy(self, centx, centy, xsize, ysize, outfile, hext=0):
      """
      Creates a new fits file that is a cutout of the original image.  For
      this method, the image center is defined by its (x,y) coordinate
      rather than (ra,dec).

      Inputs:
         centx   - x coordinate of cutout center
         centy   - y coordinate of cutout center
         xsize   - size of cutout in x direction
         ysize   - size of cutout in x direction
         outfile - name of output file
         hext    - HDU containing the image data in the input image (default=0)
      """

      print ""
      print "Input file:  %s" % self.infile
      print "Output file: %s" % outfile

      """ Read in relevant data """
      subimsize = (xsize,ysize)
      subimcent = (centx,centy)
      self.get_subim_bounds(hext,subimsize,subimcent)
      data = self.hdu[hext].data[self.suby1:self.suby2,
                                    self.subx1:self.subx2].copy()
      hdr = self.hdu[hext].header.copy()
      print "Cutout image center (x,y): (%d, %d)" % \
          (self.subcentx,self.subcenty)
      print "Cutout image size (x y): %dx%d" % \
          (self.subsizex,self.subsizey)
      print ''

      """ 
      Update the header info, including updating the CRPIXn values if they
      are present.
      """
      hdr.update('ORIG_IM','Copied from %s' % self.infile)
      hdr.update('ORIG_REG','Region in original image: [%d:%d,%d:%d]' % \
                    (self.subx1,self.subx2,self.suby1,self.suby2))
      """ Update the headers to reflect the cutout center"""
      try:
         hdr.update('CRPIX1',hdr['CRPIX1']-self.subx1)
      except:
         pass
      try:
         hdr.update('CRPIX2',hdr['CRPIX2']-self.suby1)
      except:
         pass

      """ Write to the output file and clean up"""
      pf.PrimaryHDU(data,hdr).writeto(outfile)
      print "Wrote postage stamp cutout to %s" % outfile
      del data,hdr

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
      #newhdr.update('RA',self.subimhdr['ra'])
      #newhdr.update('DEC',self.subimhdr['dec'])
      #newhdr.update('CTYPE1',self.subimhdr['ctype1'])
      #newhdr.update('CTYPE2',self.subimhdr['ctype2'])
      #newhdr.update('CRVAL1',self.subimhdr['crval1'])
      #newhdr.update('CRPIX1',self.subimhdr['crpix1'])
      #newhdr.update('CRVAL2',self.subimhdr['crval2'])
      #newhdr.update('CRPIX2',self.subimhdr['crpix2'])
      #newhdr.update('CDELT1',self.subimhdr['cdelt1'])
      #newhdr.update('CDELT2',self.subimhdr['cdelt2'])
      #newhdr.update('PC1_1',self.subimhdr['pc1_1'])
      #newhdr.update('PC1_2',self.subimhdr['pc1_2'])
      #newhdr.update('PC2_1',self.subimhdr['pc2_1'])
      #newhdr.update('PC2_2',self.subimhdr['pc2_2'])
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

   def display(self, hext=0, wtfile=None, cmap='gaia', absrange=None, siglow=1.0,
               sighigh=10.0, statsize=2048, title=None, subimdef='xy', 
               subimcent=None, subimsize=None, subimunits='pixels', 
               dispunits='pixels'):
      # NB: Need to add extent parameter for call to imshow
      """
      
      """
      print ""
      print "Input file:  %s" % self.infile

      """ Read in relevant data """
      if subimdef == 'radec':
         if subimcent == None:
            ra = None
            dec = None
         else:
            ra = subimcent[0]
            dec = subimcent[1]
         if subimsize == None:
            xsize = None
            ysize = None
         else:
            xsize = subimsize[0]
            ysize = subimsize[1]
         self.def_subim_radec(ra,dec,xsize,ysize,hext=hext)
      else:
         self.get_subim_bounds(hext,subimsize,subimcent)
         self.def_subim_xy(hext)
         #data = self.hdu[hext].data[self.suby1:self.suby2,
         #                              self.subx1:self.subx2].copy()
         print "Display image center (x,y): (%d, %d)" % \
             (self.subcentx,self.subcenty)
      print "Displayed image size (x y): %dx%d" % \
          (self.subsizex,self.subsizey)
      print ''

      """ 
      Take the absolute range to display, if requested, otherwise clip the 
      data and set display limits from the clipped values 
      """
      if absrange is not None:
         print "Taking requested absolute display limits"
         print "----------------------------------------"
         vmin = absrange[0]
         vmax = absrange[1]
         print "  vmin (absolute): %f" % vmin
         print "  vmax (absolute): %f" % vmax
      else:
         if self.found_rms == False:
            print "Calculating display limits"
            print "--------------------------"
            self.mean_clip,self.rms_clip = ccd.sigma_clip(self.subim,
                                                          verbose=True)
            self.found_rms = True
         vmin = self.mean_clip - siglow*self.rms_clip
         vmax = self.mean_clip + sighigh*self.rms_clip
         print " Clipped mean: %f" % self.mean_clip
         print " Clipped rms:  %f" % self.rms_clip
         print " vmin (mean - %2d sigma):  %f" % (siglow,vmin)
         print " vmax (mean + %2d sigma):  %f" % (sighigh,vmax)

      """ Set the color map """
      if cmap == 'gray':
         cmap = plt.cm.gray
      elif cmap == 'gray_inv':
         cmap = plt.cm.gray_r
      elif cmap == 'heat' or cmap == 'hot':
         cmap = plt.cm.hot
      elif cmap == 'Yl_Or_Br' or cmap == 'gaia':
         cmap = plt.cm.YlOrBr_r
      elif cmap == 'jet':
         cmap = plt.cm.jet
      else:
         print ' WARNING - Requested unknown color map.  Using gaia colors'
         cmap = plt.cm.YlOrBr_r

      """ Set the displayed axes to be in WCS offsets, if requested """
      coords = n.indices(self.subim.shape).astype(n.float32)
      if dispunits == 'radec':
         inhdr = self.hdu[hext].header.copy()
         # INSERT ERROR CHECKING HERE
         wcsinfo = wcs.parse_header(inhdr)
         inscale = sqrt(wcsinfo[2][0,0]**2 + wcsinfo[2][1,0]**2)*3600.
         xxsize = self.subimsize[0] / inscale
         yysize = self.subimsize[1] / inscale
         pltc   = (coords - self.subimsize/2.)*inscale
         pltc[1] *= -1.
      else:
         pltc = coords

      """ Display the image """
      maxi = self.subim.shape
      plt.imshow(self.subim,origin='bottom',cmap=cmap,vmin=vmin,vmax=vmax,
                 interpolation='none',
                 extent=(pltc[1][0,0],pltc[1][maxi[1]-1,maxi[1]-1],
                         pltc[0][0,0],pltc[0][maxi[0]-1,maxi[0]-1]))
   #plt.xlabel(r"$\Delta \alpha$ (arcsec)")
   #plt.ylabel(r"$\Delta \delta$ (arcsec)")
      if title is not None:
         plt.title(title)
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

   """ Make output RMS file, if requesed """
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

   try:
      hdu_list = open_fits(infile)
   except:
      return

   """ Get info about input image """
   inhdr = hdu_list[0].header.copy()
   xmax = inhdr["NAXIS1"]
   ymax = inhdr["NAXIS2"]
   print "poststamp: Image %s has dimensions %d x %d" % (infile,xmax,ymax)

   """ Check to make sure requested size isn't bigger than the input image """
   """ (to be done) """

   """ Make sure that everything is in integer format """
   cx = int(centx)
   cy = int(centy)
   sz = [int(xsize),int(ysize)]
   print sz
   halfx = int(xsize/2.0)
   halfy = int(ysize/2.0)

   """ Calculate output region 
   For now does not deal with regions partially outside the input file
   """
   x1 = cx - halfx
   x2 = cx + halfx
   y1 = cy - halfy
   y2 = cy + halfy

   """ Update the header with cutout information """
   inhdr.update('ORIG_IM','Copied from %s with region[%d:%d,%d:%d]' % \
                   (infile,x1,x2,y1,y2))

   """ """
   outdat = n.zeros(sz)
   outdat = hdu_list[0].data[y1:y2,x1:x2]
   outhdu = pf.PrimaryHDU(data=outdat,header=inhdr)
   outhdu.writeto(outfile)

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
   if wcs.is_degree(ra)==False:
      ra = wcs.ra2deg(ra)
   if wcs.is_degree(dec)==False:
      dec = wcs.dec2deg(dec)

   """ Calculate the (x,y) that is associated with the requested center"""
   inhdr = hdu[hext].header.copy()
   x,y = wcs.sky2pix(inhdr,ra,dec)

   """ 
   Get rough image size in pixels for the segment of input image, since the pixel
   scale for the output image does not necessarily match that of the input 
   image.
   Note that wcsinfo[2] is the CD matrix.  The rough scale is just the
   x-axis scale, assuming that the y-axis scale is the same.
   """
   wcsinfo = wcs.parse_header(inhdr)
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
   outheader = wcs.make_header(ra,dec,xsize,ysize,scale)
   coords = n.indices((ysize,xsize)).astype(n.float32)
   skycoords = wcs.pix2sky(outheader,coords[1],coords[0])
   ccdcoords = wcs.sky2pix(inhdr,skycoords[0],skycoords[1])
   coords[1] = ccdcoords[0]
   coords[0] = ccdcoords[1]
   print coords.shape

   """ Create the output HDU, while at the same time transforming the coords """
   out = pf.PrimaryHDU(ndimage.map_coordinates(data,coords,
                                               output=n.float64,order=5))
   out.header = outheader.copy()

   return out

#-----------------------------------------------------------------------

def image_cutout(infile, ra, dec, imsize, pixscale, nohist=False,
                 verbose=True):
   """
   This function takes as an input a fits file, a central position (RA,Dec),
   an image size in pixels, and a pixel scale and produces an output postage
   stamp fits file.  All of the heavy lifting is actually done by
   image_cutout_hdu, which is a very slightly modified version of Matt
   Auger's image_cutout from imagelib.py.

   Inputs:
      infile   - input fits file
      ra       - RA of postage stamp center
      dec      - RA of postage stamp center
      imsize   - image size in pixels
      pixscale - pixel scale in arcsec/pix
      nohist   - [OPTIONAL] Boolean flag set to True to delete the HISTORY cards
                 in the fits header.  Default=False.
      verbose  - [OPTIONAL] Boolean flag set to True for verbose output.
                 Default=True
   """

   if verbose:
      print ""
      print "     Input image              RA         Dec    "
      print "------------------------- ----------- ----------"
      print "%-25s %11.7f %+10.6f" % (infile,ra,dec)

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
                         rms2=None, sighigh=10., title=None, showradec=True,
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
      vmax = sighigh * rms1
   else:
      m1,s1 = ccd.sigma_clip(dat1)
      vmin = m1 - s1
      vmax = m1 + sighigh*s1

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
   #skyc   = wcs.pix2sky(hdu1.header,coords[1],coords[0])
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

def overlay_contours(infile1, infile2, ra, dec, imsize, pixscale=None, rms1=None,
                     rms2=None, sighigh=10., title=None, showradec=True,
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
                  If pixscale is None (the default) then just use the
                  native pixel scale of each of the input images.
      rms1      - user-requested rms for data in the first image. If set to 
                   None (the default) then calculate rms from the cutout data
                   themselves
      rms2      - user-requested rms for data in the second image. If set to 
                   None (the default) then calculate rms from the cutout data
                   themselves
   """

   """ Read the input images """
   try:
      im1 = Image(infile1)
   except:
      return
   print "   .... Done"
   try:
      im2 = Image(infile2)
   except:
      return
   print "   .... Done"

   """ 
   Make cutouts of the appropriate size for each of the input images
   For the first image this is done via a call to display
   """
   #im1.display(cmap='gray_inv',subimdef='radec',subimcent=(ra,dec),
   #            subimsize=(imsize,imsize))
   im1.def_subim_radec(ra,dec,imsize,outscale=pixscale)
   im2.def_subim_radec(ra,dec,imsize,outscale=pixscale)

   """ Display the first image, with axes in dRA,dDec """
   """ Set contour levels for the second image """
   if rms2 is None:
      m2,rms2 = ccd.sigma_clip(im2.subim)
   contbase = sqrt(3.)
   maxcont = int(log((im2.subim.max()/rms2),contbase))
   if maxcont < 3:
      clevs = n.array([-3.,3.,contbase**3])
   else:
      clevs = n.concatenate(([-contbase**2],
                             n.logspace(2.,maxcont,maxcont-1,base=contbase)))
   print "Contour levels: %f *" % rms2
   print clevs
   clevs *= rms2

   """ Clean up """
   im1.close()
   im2.close()
   del im1,im2


#---------------------------------------------------------------------------

def overlay_contours_old(infile1, infile2, ra, dec, imsize, pixscale, rms1=None,
                         rms2=None, sighigh=10., title=None, showradec=True,
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
   overlay_contours_hdu(hdu1,hdu2,ra,dec,imsize,pixscale,rms1,rms2,sighigh,
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

def display_image(infile,inhdu=0,cmap='gray',siglow=1.0,sighigh=10.0):
   """
   Displays the image data contained in an input fits file.
   Does this through a call to the Image class, which is returned.
   """

   # Read in the image
   image = Image(infile)

   # Display the image.  Note that for now this call does not include
   #  all of the possible parameters defined in the Image.display method
   #  (missing, e.g., wtfile, statsize, extent)
   image.display(hdu=inhdu, cmap=cmap, siglow=siglow, sighigh=sighigh)

   return image

#-----------------------------------------------------------------------

def plot_cat(fitsfile, catfile, xcol=0, ycol=1, marksize=20., markcolor='g',
             inhdu=0, cmap='gray', siglow=1.0, sighigh=10.0):
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
      siglow    - sets display range for input image
                  default value: 1.0 (1-sigma below the clipped mean)
      siglow    - sets display range for input image
                  default value: 10.0 (10-sigma below the clipped mean)
   """

   """ Plot the image """
   try:
      display_image(fitsfile,inhdu=inhdu,cmap=cmap,siglow=siglow,sighigh=sighigh)
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

def make_wcs_from_tel_pointing(inhdr, pixscale, rotatekey=None, 
                               rakey='ra', deckey='dec'):
   """
   Many fits files produced by cameras on ground-based telescopes
   do not come with full WCS header cards.  However, most do store the
   (rough) telescope pointing information in the RA and DEC header cards.
   This function uses those and the input pixel scale to create a fake
   output header.
   This function returns the temporary header, the WCS information in which 
   can either be copied into the input file or into another file that is
   using the input file as a reference.

   Inputs:
      inhdr     -  header of the input fits file
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
   """

   print ""
   """ Read in the RA and Dec pointing values """
   try:
      ra_tel = inhdr[rakey]
   except:
      print "RA designated by keyword %s was not found in header" % rakey
      print ""
      return
   try:
      dec_tel = inhdr[deckey]
   except:
      print "Dec designated by keyword %s was not found in header" % deckey
      print ""
      return

   """ Get size of image """
   try:
      xsize = inhdr['naxis1']
   except:
      print "Problem reading x-axis size"
      print ""
      return
   try:
      ysize = inhdr['naxis2']
   except:
      print "Problem reading x-axis size"
      print ""
      return

   """ Create the temporary output header """
   try:
      outhdr = wcs.make_header(ra_tel,dec_tel,xsize,ysize,pixscale)
   except:
      print "ERROR. Could not create output header"
      print ""
      return

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
      rothdr = wcs.rotate_header(outhdr,rot)
      outhdr = rothdr.copy()

   print "Created header from telescope pointing info."
   print outhdr.ascardlist()
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
   refinfo = wcs.parse_header(refwcs)

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

