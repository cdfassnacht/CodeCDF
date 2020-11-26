"""
catfuncs.py

A library containing functions that are useful for working with catalogs
of objects.  These will primarily be produced by SExtractor, but do not
necessarily have to be.

"""

import os, sys
import numpy as np
from math import pi, fabs
from matplotlib import pyplot as plt
import astropy
from astropy import units as u
from astropy.io import ascii
from astropy.table import Table, vstack
try:
   from astropy.io import fits as pf
except:
   import pyfits as pf
if astropy.__version__[:3] == '0.3':
   from astropy.coordinates import ICRS as SkyCoord
else:
   from astropy.coordinates import SkyCoord
from specim.imfuncs import image as imf
from astrom_simple import select_good_ast

# ===========================================================================

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

   def __init__(self, incat, catformat='ldac', verbose=True, namecol=None,
                racol=None, deccol=None, rafield=None, decfield=None,
                usecols=False):
      """
      This method gets called when the user types something like
         secat = Secat(infile)

      Inputs:
         incat     - input table data.  This can be in one of two forms:
                       1. a file containing the catalog (most common)
                       2. a Table instance containing the catalog data
         catformat - format of the input file, if incat is a filename.  
                     The options are:
                      ascii    - 
                      asciitab -
                      ldac     - 
                      csv      -
                      sdssfits - 
                      secat    -
                     NOTE: this parameter is not used if incat is a Table
                      rather than the name of an input file
      """

      """ Set a flag showing whether the file has been modified """
      self.modified = False

      """ Set other default values """
      self.radec = None
      self.rafield = None
      self.decfield = None
      self.centpos = None
      self.galmask = None
      self.starmask = None

      """
      Start by loading the catalog information
      """
      incattype = (str(type(incat))).split('.')[-1]
      if incattype[0:3] == 'Tab' or incattype[0:4] == 'FITS':
         self.data = incat.copy()
         if rafield:
            self.rafield = rafield
         if decfield:
            self.decfield = decfield
         self.nrows = len(incat)
         self.ncols = len(incat.columns)
         self.catformat = 'Table'

      elif isinstance(incat, str):
         self.load_from_file(incat, catformat=catformat, verbose=verbose,
                             namecol=namecol, racol=racol, deccol=deccol,
                             rafield=rafield, decfield=decfield,
                             usecols=usecols)
      else:
         print('')
         print('ERROR: input catalog must either be a filename or a Table')
         print(' Input catalog type is' + str(type(incat)))
         print('')
         return None

   #-----------------------------------------------------------------------

   def load_from_file(self, incat, catformat='ldac', verbose=True, namecol=None,
                      racol=None, deccol=None, rafield=None, decfield=None,
                      usecols=False):

      if verbose:
         print('')
         print("Loading data from catalog file %s" % incat)
         print("-----------------------------------------------")
         print("Expected catalog format: %s" % catformat)
         print('')

      """
      Define a flag for successful reading of input catalog
      """
      read_success = True

      """ Read in catalog in a manner appropriate to the catformat """
      if catformat == 'secat':
         try:
            self.data = ascii.read(incat)
            ncols = len(self.data.colnames)
            nrows = len(self.data)
            """ Set the field names """
            self.rafield  = 'ALPHA_J2000'
            self.decfield = 'DELTA_J2000'
         except:
            print("  ERROR. Problem in loading file %s" % incat)
            print("  Check to make sure filename matches an existing file.")
            print('')
            print("  This also may have failed if the input file is in the")
            print("   SExtractor FITS LDAC format.  Checking that...")
            print('')
            read_success = False

      elif catformat == 'asciitab':
         f = open(incat)
         foo = f.readline()
         f.close()
         if foo[0] == '#':
            try:
               self.data = ascii.read(incat, guess=False,
                                      format='commented_header')
            except:
               print('')
               print('ERROR: Could not read data from %s' % incat)
               print(' Tried using "commented_header" format but failed') 
               print(' Please check input file.')
               print('')
               raise IOError
         else:
            try:
               self.data = ascii.read(incat)
            except:
               print('')
               print('ERROR: Could not properly read data from %s' % incat)
               print('Tried using the automatic formatting but failed')
               print(' Please check input file.')
               print('')
               raise IOError
         ncols = len(self.data.colnames)
         nrows = len(self.data)
         """ Set the field names """
         if rafield:
            self.rafield = rafield
         if decfield:
            self.decfield = decfield

      elif catformat=='ascii':
         """ ASCII format """
         try:
            """ Set up the data format in the catalog """
            foo = np.loadtxt(incat,dtype='S30')
            ncols = foo.shape[1]
            del foo
            coltypes = np.ones(ncols,dtype='S3')
            coltypes[:] = 'f8'
            if namecol is not None:
               print("Object name column: %d" % namecol)
               coltypes[namecol] = 'S30'
            colstr = ''
            for i in range(ncols):
               colstr = '%s,%s' % (colstr,coltypes[i])
            colstr = colstr[1:]
            dt = np.dtype(colstr)

            """ Actually read in the data """
            self.informat = 'ascii'
            self.data = np.loadtxt(incat,dtype=dt)
            nrows = self.data.shape[0]

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
            print("  ERROR. Problem in loading file %s" % incat)
            print("  Check to make sure filename matches an existing file.")
            print("  ")
            print("  This may have failed if there is a string column in")
            print("   the input catalog (e.g., for an object name).  ")
            print("  If this is the case, use the namecol to indicate which ")
            print("   column contains the string values (column numbers are ")
            print("   zero-indexed)")
            print('')
            print("  This also may have failed if the input file is in the")
            print("   SExtractor FITS LDAC format.  Checking that...")
            print('')
            read_success = False

      elif catformat.lower()=='ldac' or read_success==False:
         try:
            self.data = Table.read(incat, format='fits', hdu=2)
         except:
            print("  ERROR. Problem in loading file %s" % incat)
            print("  Check to make sure filename matches an existing file.")
            print('')
            return
         self.informat = 'ldac'
         nrows = len(self.data)
         ncols = len(self.data.columns)

         """ Set the field names """
         self.rafield = 'ALPHA_J2000'
         self.decfield = 'DELTA_J2000'

      elif catformat.lower()=='csv' or read_success==False:
         try:
            self.data = Table.read(incat)
         except:
            print("  ERROR. Problem in loading file %s" % incat)
            print("  Check to make sure filename matches an existing file.")
            print('')
            return
         self.informat = 'csv'
         nrows = len(self.data)
         ncols = len(self.data.columns)

         """ Set the field names """
         self.rafield = 'raStack'
         self.decfield = 'decStack'

      elif catformat.lower()=='sdssfits' or read_success==False:
         try:
            self.data = Table.read(incat, format='fits', hdu=1)
         except:
            print("  ERROR. Problem in loading file %s" % incat)
            print("  Check to make sure filename matches an existing file.")
            print('')
            return
         self.informat = 'sdss'
         nrows = len(self.data)
         ncols = len(self.data.columns)

         """ Set the field names """
         self.rafield = 'ra'
         self.decfield = 'dec'

         """
         Split the stars from the galaxies, according to the SDSS 
         classification.
         In the SDSS scheme, type=3 is a galaxy and type=6 is a star
         """
         if 'type' in self.data.colnames:
            objclass = self.data['type']
         elif 'type_r' in self.data.colnames:
            objclass = self.data['type_r']
         else:
            objclass = None
   
         if objclass is not None:
            self.galmask = objclass == 3
            self.starmask = objclass == 6
         else:
            self.galmask = np.ones(dtype=bool)
            self.starmask = None
         print('Read SDSS catalog from %s' % incat)
         print('Number of galaxies: %5d' % self.galmask.sum())
         if self.starmask is not None:
            print('Number of stars:    %5d' % self.starmask.sum())

      else:
         print('')
         print('Unrecognized format.  Must be one of:')
         print('  ascii, secat, ldac, sdssfits')
         sys.exit()

      if verbose:
         print("Number of rows:    %d" % nrows)
         print("Number of columns: %d" % ncols)
         if self.rafield is not None:
            print('RA field name:  %s' % self.rafield)
         if self.decfield is not None:
            print('Dec field name: %s' % self.decfield)

      self.infile = incat
      self.catformat = catformat
      self.nrows = nrows
      self.ncols = ncols
      
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
            print('Updating input fits LDAC file: %s' % self.infile)
         else:
            self.hdu.close()
      else:
         print('')
         print('WARNING. Calling close_ldac but file is not in ldac format')
         print('')

   #----------------------------------------------------------------------

   def make_magmask(self, magname, mfaint=None, mbright=None):
      """
      Makes a mask that is True for magnitudes brighter than mfaint and
       fainter than mbright.
      Note that one of these two could have the value None, in which case
       it would be ignored.  For example, if mbright is None, then the mask
       will be True for all galaxies brighter than mfaint

      Inputs:
       magname - the name in the catalog for the column that represents the
                  object magnitudes.  This could be something like, e.g., 'r' 
                  or 'MAG_AUTO'
      """

      mag = self.data[magname].astype(float)
      if mfaint is None and mbright is None:
         self.magmask = np.ones(self.nrows,dtype=bool)
      elif mfaint is None:
         self.magmask = mag >= mbright
      elif mbright is None:
         self.magmask = mag <= mfaint
      else:
         self.magmask = (mag >= mbright) & (mag <= mfaint)

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

      """ Make sure that values are in the correct range """
      self.ra = self.ra.astype(float)
      self.dec = self.dec.astype(float)
      print(self.ra)
      print(type(self.ra))
      ramask = (self.ra >= 0.) & (self.ra < 360.)
      decmask = (self.dec >= -90.) & (self.dec <= 90.)
      mask = (ramask) & (decmask)
      self.ra = self.ra[mask]
      self.dec = self.dec[mask]
      self.data = self.data[mask]
      
      """ 
      Put data into a SkyCoords container for easy coordinate-based calculations
      """
      if (self.ra is not None) and (self.dec is not None):
         """ 
         Determine whether the RA coordinate is in hours or in degrees
         For now this test is made simple: if the format of the RA data is
          a string, then the data is expected to be in hms format, otherwise
          expect decimal degrees.
         """
         if type(self.ra[0]) is str or type(self.ra[0]) is np.string_:
            raunit = u.hourangle
         else:
            raunit = u.deg
         self.radec = SkyCoord(self.ra,self.dec,unit=(raunit,u.deg))

      """ 
      If radec has not been set, report this information
      (consider raising an exception in future versions of the code)
      """
      if self.radec is None:
         print('')
         print('WARNING: get_radec was called but RA and Dec information was')
         print('         not found.  Please check the format of your catalog.')
         print('')

   #----------------------------------------------------------------------

   def read_centpos(self, posfile, verbose=False):
      """
      Reads a position in from a file and converts it to an astropy SkyCoord
      format.

      NOTE: right now this is hard-wired to just read in a file in the 
      standard Keck starlist format:

          label rahr ramin rasec decdeg decamin decasec equinox

      which matches the *.pos in CDF's Lenses/Data/* directories

      Inputs:
       posfile - file containing the RA and dec
      """

      """ Read the data from the file """
      posinfo = ascii.read(posfile, 
                           names=['object', 'rahr', 'ramin', 'rasec',
                                  'decdeg', 'decamin', 'decasec', 'equinox'])

      """ Convert to SkyCoord format """
      i = posinfo[0]
      ra  = '%d:%d:%f' % (i['rahr'], i['ramin'], i['rasec'])
      dec = '%d:%d:%f' % (i['decdeg'], i['decamin'], i['decasec'])
      radec = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))

      """ Print out information if requested """
      if verbose:
         print('')
         print('Object: %s' % i['object'])
         print('Coordinates: %s %s' % \
                  (radec.ra.to_string(unit=u.hourangle, decimal=False,
                                      sep=':', precision=3, pad=True),
                   radec.dec.to_string(decimal=False, sep=':', precision=3,
                                       alwayssign=True, pad=True)))

      """ Save the output within the class """
      self.centpos = radec

   #----------------------------------------------------------------------

   def sort_by_pos(self, centpos):
      """
      Sorts the catalog in terms of distance from some central position, which
      is passed as a SkyCoord variable (from astropy.coordinates)
      """

      """ First check to see that we actual have WCS information  """
      if self.radec is None:
         print('')
         print('ERROR: sort_by_pos.  No WCS information in catalog')
         print('')
         return

      """ Otherwise, use the SkyCoords functionality to easily sort """
      sep   = self.radec.separation(centpos)
      try:
         offsets = centpos.spherical_offsets_to(self.radec)
      except:
         offsets = None
      ind   = np.argsort(sep.arcsec)
      if self.catformat == 'ascii':
         self.data = self.data[ind,:]
      else:
         print(self.catformat)
         self.data = self.data[ind]

      """ Also sort things that are outside the data table """
      self.radec   = self.radec[ind]
      self.ra      = self.ra[ind]
      self.dec     = self.dec[ind]
      self.sep     = sep[ind]
      if offsets is None:
         self.dx = None
         self.dy = None
      else:
         self.dx = (offsets[0].arcsecond)[ind]
         self.dy = (offsets[1].arcsecond)[ind]
      self.sortind = ind


   #-----------------------------------------------------------------------

   def make_reg_file(self, outfile, rcirc, color='green', fluxcol=None, 
                     fluxerrcol=None, labcol=None, labdx=0.0012, labdy=0.0012,
                     plot_high_snr=False, mask=None, snrgood=10.):
      """
      Uses the RA and Dec info in the catalog to make a region file that
      can be used with ds9.
      """

      """ 
      Start by putting the RA and Dec info into a somewhat more convenient
      format
      """

      self.get_radec()
      if self.radec is None:
         print('')
         print("ERROR: Could not read RA and Dec information from input file")
         return

      """ 
      If the flux information is given, then report on high SNR detections
      """

      ngood = 0
      snrmask = None
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
         snrmask = snr>snrgood

      """ Mask the input data if requested """
      print('Total objects in catalog:        %d' % len(self.radec))
      if mask is not None:
         radec  = self.radec[mask]
         selmask = mask
         print('Objects selected by mask:        %d' % len(radec))
      else:
         radec = self.radec
         selmask = np.ones(len(radec), dtype=bool)
      ntot   = len(radec)
      radeg  = radec.ra.degree
      decdeg = radec.dec.degree

      """ Select the high SNR objects, if requested """
      if snrmask is not None:
         selsnrmask = (selmask) & (snrmask)
         radecgood = self.radec[selsnrmask]
         ngood = len(radecgood)
         print('Of those, objects with SNR>%.1f: %d' % (snrgood,ngood))
         gradeg  = radecgood.ra.degree
         gdecdeg = radecgood.dec.degree
         
      """ Write the output region file """
      f = open(outfile,'w')
      f.write('global color=%s\n' % color)
      for i in range(ntot):
         f.write('fk5;circle(%10.6f,%+10.6f,%.1f")\n' % 
                 (radeg[i],decdeg[i],rcirc))
      if plot_high_snr and ngood>0:
         f.write('global color=red\n')
         for i in range(ngood):
            f.write('fk5;circle(%10.6f,%+10.6f,0.0011)\n' \
                       %(gradeg[i],gdecdeg[i]))

      """ Add labels if requested """
      if labcol is not None:
         lab = self.data[labcol][selmask]
         cosdec = np.cos(pi * self.dec[selmask] / 180.)
         xx = self.ra[selmask] + labdx * cosdec
         yy = self.dec[selmask] + labdy
         f.write('global color=%s\n' % color)
         for i in range(ntot):
            f.write('fk5;text(%10.6f,%+10.6f) # text={%s}\n'% \
                       (xx[i],yy[i],str(lab[i])))

      """ Wrap up """
      print("Wrote region file %s" % outfile)
      f.close()

   #-----------------------------------------------------------------------

   def _print_autoslit_infile(self, filename, galmask, starmask,
                              maskcent, magcol, smagcol, objroot=None,
                              add_lens=False, ndigits=4):
      """

      This method is called by lrismask_prep.  It prints out a file that
      is used as input for the autoslit3 code.

      """

      """ Set up the object name column """
      ind = np.arange(len(self.data)) + 1
      objname = np.zeros(len(self.data), dtype='U16')
      if objroot is not None:
         for i, obj in enumerate(ind):
            objname[i] = '%s_%04d' % (objroot, obj)
      else:
         for i, obj in enumerate(ind):
            objname[i] = '%04d' % obj

      """ Set up the output table """
      outnames = ['name', 'priority', 'mag', 'rahr', 'ramin', 'rasec', 'decdeg',
                  'decamin', 'decasec', 'epoch', 'equinox', 'pm_ra', 'pm_dec']
      tmpgal = self.data[galmask]
      radec = self.radec[galmask]
      pri = np.ones(len(tmpgal)) * 100
      tmpgal[magcol][tmpgal[magcol]<0.] = 99.
      epoch = np.ones(len(tmpgal)) * 2000.
      pm = np.zeros(len(tmpgal))
      galtab = Table([objname[galmask], pri, tmpgal[magcol], radec.ra.hms.h,
                      radec.ra.hms.m, radec.ra.hms.s, radec.dec.dms.d,
                      np.fabs(radec.dec.dms.m), np.fabs(radec.dec.dms.s),
                      epoch, epoch, pm, pm], names=outnames)

      """ Set up the guidestar table """
      tmpstar = self.data[starmask]
      radec = self.radec[starmask]
      pri = np.ones(len(tmpstar)) * -1
      epoch = np.ones(len(tmpstar)) * 2000.
      pm = np.zeros(len(tmpstar))
      starname = np.zeros(len(tmpstar), dtype='U16')
      starind = ind[starmask]
      for i, si in enumerate(starind):
         starname[i] = 'S%04d' % si
      startab = Table([starname, pri, tmpstar[smagcol], radec.ra.hms.h,
                       radec.ra.hms.m, radec.ra.hms.s, radec.dec.dms.d,
                       np.fabs(radec.dec.dms.m), np.fabs(radec.dec.dms.s),
                       epoch, epoch, pm, pm], names=outnames)

      """
      Make a mini-table containing the mask center and, optionally, the lens.
      These one or two lines will be at the top of the output file
      """
      if add_lens:
         hdrrows = 2
      else:
         hdrrows = 1
      hdrtab = Table(np.zeros((hdrrows, (len(outnames)-1))),
                     names=outnames[1:])
      hdrtab.add_column(galtab['name'][:hdrrows], index=0)
      ra = maskcent.ra
      dec = maskcent.dec
      hdrtab[0] = ['CENTER', 9999, 0., ra.hms.h, ra.hms.m, ra.hms.s,
                   dec.dms.d, fabs(dec.dms.m), fabs(dec.dms.s),
                   2000., 2000., 0., 0.]
      if add_lens:
         cra = self.centpos.ra.hms
         cdec = self.centpos.dec.dms
         hdrtab[1] = ['Lens', 900, 20., cra.h, cra.m, cra.s, cdec.d,
                      fabs(cdec.m), fabs(cdec.s), 2000., 2000., 0., 0.]

      """ Combine the two tables """
      outtab = vstack([hdrtab, galtab, startab])
      
      """ Set up the formatting """
      outtab['priority'].format = '%4d'
      outtab['mag'].format = '%5.2f'
      outtab['rahr'].format = '%02d'
      outtab['ramin'].format = '%02d'
      outtab['rasec'].format = '%06.3f'
      outtab['decdeg'].format = '%+03d'
      outtab['decamin'].format = '%02d'
      outtab['decasec'].format = '%05.2f'
      outtab['epoch'].format = '%6.1f'
      outtab['equinox'].format = '%6.1f'
      outtab['pm_ra'].format = '%3.1f'
      outtab['pm_dec'].format = '%3.1f'
      print('')
      print(outtab)

      """ Save the table to the output file """
      print('')
      print('Writing selected objects to %s' % filename)
      outtab.write(filename, format='ascii.no_header', overwrite=True)
      
   #-----------------------------------------------------------------------

   def lrismask_prep(self, maskcent, PA, mask_w=3., mask_h=7., galmagcol='r',
                     galmagrange=None, galmask='default', starmask='default',
                     smagcol='g', smaglim=[17., 19.], galreg='maskobj.reg',
                     objcolor='green', rcirc=1., starreg='maskstar.reg',
                     add_lens=True, outfile=None, objroot=None):
      """

      Code to select objects for a LRIS slitmask

      """

      """
      Make sure that the catalog has been sorted, since this code needs to
      use the self.dx and self.dy arrays
      """
      if self.dx is None:
         print('')
         print('Position offsets are missing.  You must run sort_by_pos first')
         print('')
         return

      """ Get the offset of the mask center from the origin of (dx,dy) """
      maskcoord = SkyCoord(maskcent[0], maskcent[1], unit=(u.deg, u.deg))
      maskxy = self.centpos.spherical_offsets_to(maskcoord)
      print('Requested mask center has (dx,dy) = (%+7.2f, %+7.2f)' %
            (maskxy[0].arcsec, maskxy[1].arcsec))

      """ 
      Translate and rotate the object (dx,dy) to a masked-centered frame in
      which the mask is rectilinear.  That is, first make the mask center
      the origin and then rotate by the NEGATIVE of the mask PA.
      """
      rot_ang = -1. * PA * pi / 180.
      xtmp = (self.dx - maskxy[0].arcsec)
      ytmp = (self.dy - maskxy[1].arcsec)
      xx = xtmp * np.cos(rot_ang) + ytmp * np.sin(rot_ang)
      yy = -xtmp * np.sin(rot_ang) + ytmp * np.cos(rot_ang)

      """ Select only the objects that lie within the mask """
      star_extra = 75.
      xbound = 0.5 * mask_w * 60.
      ybound = 0.5 * mask_h * 60.
      objmask = (xx >= -xbound) & (xx <= xbound) & (yy >= -ybound) \
         & (yy <= ybound)
      starsonmask = (xx >= -xbound-star_extra) & (xx <= xbound+star_extra) & \
         (yy >= -ybound) & (yy <= ybound)

      """ Choose the galaxies and stars that fall within the mask """
      if galmask == 'default':
         if self.galmask is not None:
            totgalmask = np.logical_and(self.galmask, objmask)
         else:
            totgalmask = objmask
      elif galmask is not None:
         totgalmask = np.logical_and(galmask, objmask)
      else:
         totgalmask = objmask
         
      if starmask == 'default':
         if self.starmask is not None:
            totstarmask = np.logical_and(self.starmask, starsonmask)
         else:
            totstarmask = starsonmask
      elif starmask is not None:
         totstarmask = np.logical_and(starmask, starsonmask)
      else:
         totstarmask = starsonmask
      stardata = self.data[smagcol].astype(float)
      guidestarmask = totstarmask & (stardata >= smaglim[0]) & \
         (stardata <= smaglim[1])

      """ Make ds9 region files for diagnostic checks """
      print('')
      print('Making ds9 region file for possible targets for slitmask')
      self.make_reg_file(galreg, rcirc, color=objcolor, mask=totgalmask)
      print('')
      print('Making ds9 region file for possible guide stars for slitmask')
      self.make_reg_file(starreg, rcirc*1.5, color='red', mask=guidestarmask)

      """ Write out the autoslit3 input file """
      if outfile is not None:
         self._print_autoslit_infile(outfile, totgalmask, guidestarmask,
                                     maskcoord, galmagcol, smagcol,
                                     objroot=objroot, add_lens=add_lens)
                     
   #-----------------------------------------------------------------------

   #def plot_radec(self, symb='bo'):
   #-----------------------------------------------------------------------

   def plot_fwhm(self, fwhmcol='FWHM_IMAGE', magcol='MAG_AUTO', 
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

      try:
         fwhm = self.data[fwhmcol]
      except KeyError:
         print('')
         print('Catalog does not contain a %d column' % fwhmcol)
         print('')
         return
      try:
         mag = self.data[magcol]
      except KeyError:
         print('')
         print('Catalog does not contain a %d column' % magcol)
         print('')
         return
      plt.plot(fwhm,mag,'bo')
      plt.xlim(xlim)
      plt.ylim(ylim)
      plt.xlabel('FWHM (pixels)')
      plt.ylabel('Magnitude')
      plt.show()

   #-----------------------------------------------------------------------

   def plot_nhist(self, magcol='MAG_AUTO', usestarmask=False,
                  magmin=15, magmax=28, color='b', alpha=1.):
      """
      Plots a histogram of galaxy magnitudes (similar to a log N-log S plot)
      that can be used to determine the magnitude to which the catalog is 
      complete.  A minimum FWHM can be set in order to select objects that
      are likely to be galaxies, but this is not required.

      Inputs:
         magcol      - column containing the magnitudes. Default = 'mag_auto'
         usestarmask - when this parameter is set to True, then use the
                       starmask mask to select the galaxies.
                       NOTE: this means that the set_starmask method has to
                       have been run for each of the input catalogs, or
                       else all objects will be plotted
                       Default=False
         magmin      - minimum magnitude to use for the plot. Default=15
         magmax      - maximum magnitude to use for the plot. Default=28
      """

      """ Get the magnitudes to be plotted """
      if usestarmask:
         if self.starmask is None:
            print('')
            print('WARNING: you have set usestarmask=True but there are')
            print(' no stars selected by the starmask')
            print('Please make sure that set_starmask has been run BEFORE'
                  'runing plot_nhist')
            print(' if you want to use this mask')
            print('')
            return
         else:
            mag = (self.data[self.starmask==False][magcol]).copy()
      else:
         mag = self.data[magcol].copy()

      mag = mag[np.isfinite(mag)]

      """ Plot the histogram """
      nbins = int(2 * (magmax - magmin))
      plt.hist(mag,range=(magmin, magmax), bins=nbins, color=color, 
               alpha=alpha)

      del(mag)

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

      print('')
      print("Matching catalogs: basic info")
      print("--------------------------------------------")
      print(" Catalog 1: %d coordinates" % self.ra.size)
      print(" Catalog 2: %d coordinates" % ra2.size)

      """ Initialize containers for output information """
      ramatch  = np.zeros(self.ra.size)
      decmatch = np.zeros(self.ra.size)
      self.nmatch   = np.zeros(self.ra.size,dtype=int)
      self.matchdx  = np.zeros(self.ra.size)
      self.matchdy  = np.zeros(self.ra.size)
      self.indmatch = np.ones(self.ra.size,dtype=int) * -1

      """ Correct for known shifts """
      ra2 = ra2.copy() + dra2/(3600.*np.cos(dec2))
      dec2 = dec2.copy() + ddec2/3600.

      """ Loop over catalog """
      print('')
      print("Searching for matches...")
      print("------------------------------")
      for i in range(self.ra.size):
         dx,dy = coords.sky_to_darcsec(self.ra[i],self.dec[i],ra2,dec2)
         dpos = np.sqrt(dx**2 + dy**2)
         isort = np.argsort(dpos)
         if dpos[isort[0]]<=rmatch:
            ramatch[i]  = self.ra[i]
            decmatch[i] = self.dec[i]
            self.matchdx[i]  = dx[isort[0]]
            self.matchdy[i]  = dy[isort[0]]
            self.nmatch[i]   = dpos[dpos<=rmatch].size
            self.indmatch[i] = isort[0]
         del dx,dy,dpos
      print(" Number of matches between the catalogs:  %d" % 
          (self.nmatch>0).sum())
      mra  = ramatch[self.nmatch>0]
      mdec = decmatch[self.nmatch>0]
      mdx  = self.matchdx[self.nmatch>0]
      mdy  = self.matchdy[self.nmatch>0]
      mdx0 = np.median(mdx)
      mdy0 = np.median(mdy)
      print(" Median offset for matches (RA):  %+6.2f arcsec" % mdx0)
      print(" Median offset for matches (Dec): %+6.2f arcsec" % mdy0)

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
         plt.plot(np.array([mdx0]),np.array([mdy0]),'r*',ms=20)
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
         print('')
         print("Printing to file for use in ccmap:  %s" % outfile)
         print('')
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

      self.matchind = np.zeros(xast.size, dtype=int)

      xfield = 'f%d' % xcol
      yfield = 'f%d' % ycol

      for i in range(xast.size):
         dx = xast[i] - self.data[xfield]
         dy = yast[i] - self.data[yfield]
         dpos = dx**2 + dy**2
         sindex = np.argsort(dpos)
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
            print('')
            print('Pass %d' % (i+1))
            print('------------------------')
         xa0 = xa - dxmed
         ya0 = ya - dymed
         self.find_closest_xy(xa0,ya0,xcol,ycol)
         dxmed = np.median(self.matchdx)
         dymed = np.median(self.matchdy)
         if max_offset is not None:
            dpos = np.sqrt(self.matchdx**2 + self.matchdy**2)
            goodmask = dpos<max_offset
            if verbose:
               print("Applying a maximum offset cut of %7.1f pixels"
                     % max_offset)
               print("Median shifts before clipping: %7.2f %7.2f"
                     % (dxmed,dymed))
         else:
            goodmask = np.ones(xa.size,dtype=bool)
         dxm  = self.matchdx[goodmask]
         dym  = self.matchdy[goodmask]
         dxmed = np.median(dxm)
         dymed = np.median(dym)
         if verbose:
            print("Median shifts after pass:   %7.2f %7.2f" % (dxmed,dymed))

      """
      Transfer information into object and clean up
      """
      if verbose:
         print('')
         print('Found %d astrometric objects within FOV of image' % xa.size)
         print('Matched %d objects to astrometric catalog.' % dxm.size)
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
         print("Running match_fits_to_ast with:")
         print("   fitsfile = %s" % fitsfile)
         print("   astcat   = %s" % astcat)
      self.infits = fitsfile
      self.matchcat = astcat

      """
      Start by opening the fits file and reading the appropriate columns from
      the catalogs
      """
      hdulist = imf.open_fits(fitsfile)
      hdr = hdulist[imhdu].header
      if verbose:
         print('')
         hdulist.info()

      """
      Select the astrometric catalog objects that fall within the fits file FOV
      (at least with its current WCS)
      """
      raa,deca,xa,ya,astmask = select_good_ast(astcat,hdr,racol,deccol,edgedist)
      if verbose:
         print('Found %d astrometric objects within FOV of image' % raa.size)
         
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
         dxmed = np.median(self.matchdx)
         dymed = np.median(self.matchdx)
         plt.figure()
         plt.scatter(self.matchdx,self.matchdy)
         plt.xlabel('x offset (pix)')
         plt.ylabel('y offset (pix)')
         plt.axhline(color='k')
         plt.axvline(color='k')
         plt.axvline(dxmed,color='r')
         plt.axhline(dymed,color='r')
         print('')
         print("Black lines represent x=0 and y=0 axes")
         print("Red lines show median offsets of dx_med=%7.2f and "
               "dy_med=%7.2f" % (dxmed,dymed))
         #plt.show()

      """ Write the output file, in a format appropriate for input to ccmap """
      if outfile is not None:
         self.print_ccmap(outfile,verbose)

      """ Clean up """
      hdulist.close()
      del hdr,raa,deca,xa,ya,astmask

   # -----------------------------------------------------------------------

   def set_starmask(self, mask):
      """
      Takes the input mask and assigns it to an internal mask associated
      with this instance of the Secat class.
      The internal mask is called starmask and has values of True for objects
      that have been identified as stars by the provided external mask.

      The external mask will be based on some characteristics in the catalog.
      Examples could be objects with a SExtractor CLASS_STAR value greater
       than 0.7, or a FWHM_IMAGE less than a certain value, or ...
      """

      self.starmask = mask

   # -----------------------------------------------------------------------

   def set_galmask(self, mask):
      self.galmask = mask

#------------------------------------------------------------------------------

