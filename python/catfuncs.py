"""
catfuncs.py

A library containing functions that are useful for working with catalogs
of objects.  These will primarily be produced by SExtractor, but do not
necessarily have to.

"""

import numpy as n
import pyfits as pf
import imfuncs as imf
from matplotlib import pyplot as plt
from math import pi
from astrom_simple import select_good_ast

#------------------------------------------------------------------------------

class Secat:

   """

   The __init__ method has been changed to return something like a record
   array, which has the same number of rows as the old 2D float array, but
   which stores each row as a single tuple.  It is thus a 1D array, sort of
   like a structure array in C.  The columns can be accessed by field name,
   which for now is just 'f0', 'f1', etc., unless the input catalog is in
   SExtractor's FITS LDAC format, in which case the field names actually
   correspond to the SExtractor variable names.

   The code used to expect the old 2D float array format.  It should have
   all been updated, but there may still be some issues.

   """

   def __init__(self, infile, catformat='ldac', verbose=True, namecol=None,
                racol=None, deccol=None):
      """
      This method gets called when the user types something like
         secat = Secat(infile)

      Inputs:
         infile - input file containing the catalog
      """

      """
      Define a flag for successful reading of input catalog
      """

      read_success = True

      """ Set a flag showing whether the file has been modified """
      self.modified = False

      """
      Start by loading the catalog information
      """
      if verbose:
         print ""
         print "Loading data from catalog file %s" % infile
         print "-----------------------------------------------"
         print "Expected catalog format: %s" % catformat
         print ""

      """ ASCII format """
      if catformat=='ascii':
         try:
            """ Set up the data format in the catalog """
            foo = n.loadtxt(infile,dtype='S30')
            ncols = foo.shape[1]
            del foo
            coltypes = n.ones(ncols,dtype='S3')
            coltypes[:] = 'f8'
            if namecol is not None:
               print "Object name column: %d" % namecol
               coltypes[namecol] = 'S30'
            colstr = ''
            for i in range(ncols):
               colstr = '%s,%s' % (colstr,coltypes[i])
            colstr = colstr[1:]
            dt = n.dtype(colstr)

            """ Actually read in the data """
            self.informat = 'ascii'
            self.data = n.loadtxt(infile,dtype=dt)

            """ Set the field names """
            if racol is not None:
               self.rafield = 'f%d' % racol
            else:
               self.rafield = None
            if deccol is not None:
               self.decfield = 'f%d' % deccol
            else:
               self.decfield = None
            if namecol is not None:
               self.namefield = 'f%d' % namecol
            else:
               self.namefield = None

         except:
            print "  ERROR. Problem in loading file %s" % infile
            print "  Check to make sure filename matches an existing file."
            print "  "
            print "  This may have failed if there is a string column in"
            print "   the input catalog (e.g., for an object name).  "
            print "  If this is the case, use the namecol to indicate which "
            print "   column contains the string values (column numbers are "
            print "   zero-indexed)"
            print ""
            print "  This also may have failed if the input file is in the"
            print "   SExtractor FITS LDAC format.  Checking that..."
            print ""
            read_success = False

      """ FITS LDAC format """
      if catformat.lower()=='ldac' or read_success==False:
         try:
            self.hdu = pf.open(infile,mode='update')
         except:
            print "  ERROR. Problem in loading file %s" % infile
            print "  Check to make sure filename matches an existing file."
            print "  "
            return

         self.informat = 'ldac'
         self.data = self.hdu[2].data
         ncols = self.hdu[2].header['tfields']

         """ Set the field names """
         self.rafield = 'alpha_j2000'
         self.decfield = 'delta_j2000'

      if verbose:
         print "Number of rows:    %d" % self.data.shape[0]
         print "Number of columns: %d" % ncols
      self.infile = infile

   #-----------------------------------------------------------------------

   def close_ldac(self):
      """
      Closes the catalog.  If the catalog is in fits format and it has
      been modified (as shown by the modified parameter in this Secat class)
      then use flush rather than close.
      """

      """ Close the file if it is in the expected format """
      if self.informat == 'ldac':
         if self.modified:
            self.hdu.flush()
            print 'Updating input fits LDAC file: %s' % self.infile
         else:
            self.hdu.close()
      else:
         print ''
         print 'WARNING. Calling close_ldac but file is not in ldac format'
         print ''

   #-----------------------------------------------------------------------

   def make_reg_file(self, outfile, fluxcol=None, fluxerrcol=None, labcol=None, 
                     plot_high_snr=False):
      """
      Uses the RA and Dec info in the catalog to make a region file that
      can be used with ds9.
      """

      """ 
      Start by putting the RA and Dec info into a somewhat more convenient
      format
      """

      self.get_radec()
      if self.ra is None:
         print ""
         print "ERROR: No columns for RA or Dec given."
         print "Cannot make region file unless those columns are specified."
         return

      """ 
      If the flux information is given, then report on high SNR detections
      """

      ngood = 0
      if fluxcol is not None and fluxerrcol is not None:
         if self.informat == 'ldac':
            if type(fluxcol) is int:
               flux = self.data.field(fluxcol)
               fluxerr = self.data.field(fluxerrcol)
            else:
               flux = self.data[fluxcol]
               fluxerr = self.data[fluxerrcol]
         else:
            flux = self.data['f%d' %fluxcol]
            fluxerr = self.data['f%d' % fluxerrcol]
         snr = flux / fluxerr
         ragood  = self.ra[snr>10.]
         decgood = self.dec[snr>10.]
         ntot = self.ra.size
         ngood = ragood.size
         print "Out of %d total objects, %d have SNR>10" %(ntot,ngood)

      """ Write the output region file """
      f = open(outfile,'w')
      f.write('global color=green\n')
      for i in range(self.ra.size):
         f.write('fk5;circle(%10.6f,%+10.6f,0.0007)\n'% (self.ra[i],self.dec[i]))
      if plot_high_snr and ngood>0:
         f.write('global color=red\n')
         for i in range(ragood.size):
            f.write('fk5;circle(%10.6f,%+10.6f,0.0011)\n' \
                       %(ragood[i],decgood[i]))

      """ Add labels if requested """
      if labcol is not None:
         if self.informat == 'ldac':
            if type(labcol) is int:
               lab = self.data.field(labcol)
            else:
               lab = self.data[labcol]
         else:
            lab = self.data['f%d' % labcol]
         cosdec = n.cos(pi * self.dec / 180.)
         xx = self.ra + 0.0012 * cosdec
         yy = self.dec + 0.0012
         f.write('global color=green\n')
         for i in range(self.ra.size):
            f.write('fk5;text(%10.6f,%+10.6f) # text={%s}\n'% \
                       (xx[i],yy[i],str(lab[i])))

      """ Wrap up """
      print "Wrote region file %s" % outfile
      f.close()


   #-----------------------------------------------------------------------

   def get_radec(self):
      """
      Extracts the RA and Dec information from the data container.  This
      is not necessary, and some of the data containers may not even have
      WCS info, but extracting the coordinates if it does simplifies some
      later tasks.
      """

      """ Extract the information into the new containers """
      self.ra = None
      self.dec = None
      if self.rafield is not None:
         try:
            self.ra  = self.data[self.rafield].copy()
         except:
            try:
               self.ra = self.data['x_world'].copy()
            except:
               self.ra = None
      if self.decfield is not None:
         try:
            self.dec = self.data[self.decfield].copy()
         except:
            try:
               self.dec = self.data['y_world'].copy()
            except:
               self.dec = None


   #-----------------------------------------------------------------------

   #def plot_radec(self, symb='bo'):
   #-----------------------------------------------------------------------

   def plot_fwhm(self, fwhmcol='fwhm_image', magcol='mag_auto', 
                 xlim=(0,15), ylim=(28,16)):
      """
      Plots FWHM vs. magnitude.  This can be used to find the stellar locus
      and, thus, determine the seeing.

      Inputs:
         fwhmcol - column name for the FWHM data.  Default = 'fwhm_image'
         magcol  - column name for the magnitude data.  Default = 'mag_auto'
         xlim    - initial limits for FWHM axis on plot.  Default = (0,15)
         ylim    - initial limits for mag axis on plot.  Default = (28,16)
      """

      plt.plot(self.data[fwhmcol],self.data[magcol],'bo')
      plt.xlim(xlim)
      plt.ylim(ylim)
      plt.xlabel('FWHM (pixels)')
      plt.ylabel('Magnitude')
      plt.show()

   #-----------------------------------------------------------------------

   def plot_nhist(self, magcol='mag_auto', fwhmmin=0., fwhmcol='fwhm_image', 
                  magmin=15, magmax=28):
      """
      Plots a histogram of galaxy magnitudes (similar to a log N-log S plot)
      that can be used to determine the magnitude to which the catalog is 
      complete.  A minimum FWHM can be set in order to select objects that
      are likely to be galaxies, but this is not required.

      Inputs:
         magcol  - column containing the magnitudes. Default = 'mag_auto'
         fwhmmin - minimum FWHM to be used as a proxy for a star-galaxy
                   separation.  The default (0.) selects all object in the
                   catalog.
         fwhmcol - column containing the fwhm info. Default = 'fwhm_image'
         magmin  - minimum magnitude to use for the plot. Default=15
         magmax  - maximum magnitude to use for the plot. Default=28
      """

      """ Get the magnitudes to be plotted """
      if fwhmmin>0.:
         mag = self.data[magcol][self.data[fwhmcol]>fwhmmin]
      else:
         mag = self.data[magcol]

      """ Plot the histogram """
      nbins = int(2 * (magmax - magmin))
      plt.hist(mag,range=(magmin,magmax),bins=nbins)

   #-----------------------------------------------------------------------

   def match_radec(self, ra2, dec2, rmatch, dra2=0., ddec2=0., doplot=True):
      """

      *** UNDER CONSTRUCTION!  DO NOT USE YET. ***

      Given a list of ra,dec coordinates (ra2, dec2), possibly from a second
      catalog, and a match tolerance, find the matches to the catalog
      contained in self.data

      Inputs:
        ra2      - RA (decimal degrees) for catalog 
        dec2     - Dec (decimal degrees) for second catalog
        rmatch   - max distance for a valid match (arcsec)
        dra2     - optional offset in ARCSEC to apply to ra2, if there is a known
                   offset between the catalogs (default=0.0)
        ddec2    - optional offset in ARCSEC to apply to dec2, if there is a
                   known offset between the catalogs (default=0.0)
      """

      print ""
      print "Matching catalogs: basic info"
      print "--------------------------------------------"
      print " Catalog 1: %d coordinates" % self.ra.size
      print " Catalog 2: %d coordinates" % ra2.size

      """ Initialize containers for output information """
      ramatch  = n.zeros(self.ra.size)
      decmatch = n.zeros(self.ra.size)
      self.nmatch   = n.zeros(self.ra.size,dtype=int)
      self.matchdx  = n.zeros(self.ra.size)
      self.matchdy  = n.zeros(self.ra.size)
      self.indmatch = n.ones(self.ra.size,dtype=int) * -1

      """ Correct for known shifts """
      ra2 = ra2.copy() + dra2/(3600.*n.cos(dec2))
      dec2 = dec2.copy() + ddec2/3600.

      """ Loop over catalog """
      print ""
      print "Searching for matches..."
      print "------------------------------"
      for i in range(self.ra.size):
         dx,dy = coords.sky_to_darcsec(self.ra[i],self.dec[i],ra2,dec2)
         dpos = n.sqrt(dx**2 + dy**2)
         isort = n.argsort(dpos)
         if dpos[isort[0]]<=rmatch:
            ramatch[i]  = self.ra[i]
            decmatch[i] = self.dec[i]
            self.matchdx[i]  = dx[isort[0]]
            self.matchdy[i]  = dy[isort[0]]
            self.nmatch[i]   = dpos[dpos<=rmatch].size
            self.indmatch[i] = isort[0]
         del dx,dy,dpos
      print " Number of matches between the catalogs:  %d" % \
          (self.nmatch>0).sum()
      mra  = ramatch[self.nmatch>0]
      mdec = decmatch[self.nmatch>0]
      mdx  = self.matchdx[self.nmatch>0]
      mdy  = self.matchdy[self.nmatch>0]
      mdx0 = n.median(mdx)
      mdy0 = n.median(mdy)
      print " Median offset for matches (RA):  %+6.2f arcsec" % mdx0
      print " Median offset for matches (Dec): %+6.2f arcsec" % mdy0

      """ Plot up some offsets, if desired """
      if doplot:
         plt.figure(1)
         plt.scatter(mdx,mdy)
         plt.axis('scaled')
         plt.xlabel(r'$\Delta \alpha$ (arcsec)')
         plt.ylabel(r'$\Delta \delta$ (arcsec)')
         plt.title('Offsets between matched sources (rmatch = %5.2f)' % rmatch)
         plt.axvline(0.0,color='r')
         plt.axhline(0.0,color='r')
         plt.plot(n.array([mdx0]),n.array([mdy0]),'r*',ms=20)
         plt.xlim(-1.1*rmatch,1.1*rmatch)
         plt.ylim(-1.1*rmatch,1.1*rmatch)

         plt.figure(2)
         #
         ax1 = plt.subplot(221)
         plt.scatter(mra,mdy)
         plt.setp(ax1.get_xticklabels(), visible=False)
         plt.ylabel(r'$\Delta \delta$ (arcsec)')
         plt.axhline(0.0,color='r')
         #
         ax2 = plt.subplot(223, sharex=ax1)
         plt.scatter(mra,mdx)
         plt.xlabel(r'$\alpha$')
         plt.ylabel(r'$\Delta \alpha$ (arcsec)')
         plt.axhline(0.0,color='r')
         #
         ax3 = plt.subplot(222, sharey=ax1)
         plt.scatter(mdec,mdy)
         plt.axhline(0.0,color='r')
         plt.setp(ax3.get_xticklabels(), visible=False)
         plt.setp(ax3.get_yticklabels(), visible=False)
         #
         ax4 = plt.subplot(224)
         plt.scatter(mdec,mdx)
         plt.xlabel(r'$\delta$')
         plt.axhline(0.0,color='r')
         plt.setp(ax4.get_yticklabels(), visible=False)

         plt.show()

      """ Clean up """
      del ramatch,decmatch
      del mdx,mdy,mra,mdec

   #-----------------------------------------------------------------------

   def print_ccmap(self, outfile, verbose=True):
      """
      Prints out a file that can be used as the input for the pyraf ccmap
      task.  This file has 4 columns:  x  y  RA  Dec

      Inputs:
         outfile   -  output file to be used as input for ccmap
         verbose   -  print task info
      """
      if verbose:
         print ""
         print "Printing to file for use in ccmap:  %s" % outfile
         print ""
      f = open(outfile,'w')
      f.write('# (x,y) catalog: %s\n' % self.infile)
      f.write('# Astrometric catalog: %s\n' % self.matchcat)
      f.write('# Columns are x y RA Dec\n')
      for i in range(self.nmatch):
         f.write('%8.2f %8.2f  %11.7f %+11.7f\n' % \
                    (self.matchx[i],self.matchy[i],self.matchra[i],
                     self.matchdec[i]))
      f.close()

   #-----------------------------------------------------------------------

   def find_closest_xy(self, xast, yast, xcol, ycol):

      """
      Finds the closest match, in (x,y) space to each member of the astrometric
      catalog (represented by xast,yast).
      """

      self.matchind = n.zeros(xast.size, dtype=int)

      xfield = 'f%d' % xcol
      yfield = 'f%d' % ycol

      for i in range(xast.size):
         dx = xast[i] - self.data[xfield]
         dy = yast[i] - self.data[yfield]
         dpos = dx**2 + dy**2
         sindex = n.argsort(dpos)
         self.matchind[i] = sindex[0]

      self.matchdx = xast - self.data[xfield][self.matchind]
      self.matchdy = yast - self.data[yfield][self.matchind]

   #-----------------------------------------------------------------------

   def match_xy(self, xa, ya, max_offset=None, xcol=8, ycol=9, verbose=True):

      """
      Find the closest match to each astrometric catalog object and calculate the
      offsets.
      Do two loops, to deal with possible confusion of sources on first pass
      through
      """
      dxmed = 0
      dymed = 0
      for i in range(2):
         if verbose:
            print ''
            print 'Pass %d' % (i+1)
            print '------------------------'
         xa0 = xa - dxmed
         ya0 = ya - dymed
         self.find_closest_xy(xa0,ya0,xcol,ycol)
         dxmed = n.median(self.matchdx)
         dymed = n.median(self.matchdy)
         if max_offset is not None:
            dpos = n.sqrt(self.matchdx**2 + self.matchdy**2)
            goodmask = dpos<max_offset
            if verbose:
               print "Applying a maximum offset cut of %7.1f pixels" % max_offset
               print "Median shifts before clipping: %7.2f %7.2f" % (dxmed,dymed)
         else:
            goodmask = n.ones(xa.size,dtype=bool)
         dxm  = self.matchdx[goodmask]
         dym  = self.matchdy[goodmask]
         dxmed = n.median(dxm)
         dymed = n.median(dym)
         if verbose:
            print "Median shifts after pass:   %7.2f %7.2f" % (dxmed,dymed)

      """
      Transfer information into object and clean up
      """
      if verbose:
         print ''
         print 'Found %d astrometric objects within FOV of image' % xa.size
         print 'Matched %d objects to astrometric catalog.' % dxm.size
      self.nmatch = dxm.size
      self.goodmask = goodmask.copy()
      self.matchind = self.matchind[goodmask]
      del xa0,ya0,goodmask

   #-----------------------------------------------------------------------

   def match_fits_to_ast(self, fitsfile, astcat, outfile=None, max_offset=None, 
                         racol=1, deccol=2, xcol=8, ycol=9, 
                         doplot=True, edgedist=50., imhdu=0, verbose=True):

      """
      Given the fits file from which this object (self) was defined 
      and an astrometric catalog, find the closest matches of the 
      astrometric objects to those contained in this object, using the WCS 
      information in the fits header.
      """

      if(verbose):
         print "Running match_fits_to_ast with:"
         print "   fitsfile = %s" % fitsfile
         print "   astcat   = %s" % astcat
      self.infits = fitsfile
      self.matchcat = astcat

      """
      Start by opening the fits file and reading the appropriate columns from
      the catalogs
      """
      hdulist = imf.open_fits(fitsfile)
      hdr = hdulist[imhdu].header
      if verbose:
         print ""
         hdulist.info()

      """
      Select the astrometric catalog objects that fall within the fits file FOV
      (at least with its current WCS)
      """
      raa,deca,xa,ya,astmask = select_good_ast(astcat,hdr,racol,deccol,edgedist)
      if verbose:
         print 'Found %d astrometric objects within FOV of image' % raa.size
         
      """
      Find the closest match to each astrometric catalog object
      """
      self.match_xy(xa,ya,max_offset,xcol,ycol,verbose)

      """ Transfer info about matches into the object """
      xfield = 'f%d' % xcol
      yfield = 'f%d' % ycol
      self.astmask  = astmask.copy()
      self.matchx   = self.data[xfield][self.matchind].copy()
      self.matchy   = self.data[yfield][self.matchind].copy()
      self.matchra  = raa[self.goodmask].copy()
      self.matchdec = deca[self.goodmask].copy()
      self.matchdx  = self.matchdx[self.goodmask]
      self.matchdy  = self.matchdy[self.goodmask]

      """ Plot the offsets if desired """
      if doplot:
         dxmed = n.median(self.matchdx)
         dymed = n.median(self.matchdx)
         plt.figure()
         plt.scatter(self.matchdx,self.matchdy)
         plt.xlabel('x offset (pix)')
         plt.ylabel('y offset (pix)')
         plt.axhline(color='k')
         plt.axvline(color='k')
         plt.axvline(dxmed,color='r')
         plt.axhline(dymed,color='r')
         print ""
         print "Black lines represent x=0 and y=0 axes"
         print "Red lines show median offsets of dx_med=%7.2f and dy_med=%7.2f" \
             % (dxmed,dymed)
         #plt.show()

      """ Write the output file, in a format appropriate for input to ccmap """
      if outfile is not None:
         self.print_ccmap(outfile,verbose)

      """ Clean up """
      hdulist.close()
      del hdr,raa,deca,xa,ya,astmask

#------------------------------------------------------------------------------

