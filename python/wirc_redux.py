"""
wirc_redux.py

Author: Chris Fassnacht

Functions to aid in the processing of data obtained with the WIRC instrument 
on the Palomar 200-Inch Telescope. 

"""

import numpy as n
import pyfits as pf
import imfuncs as im
import astrom_scamp as ast
import coords, os
import wcs as wcsmwa

def wcs_clean(infiles):
   """
   Function to do the cleaning.

   This function does two things to the input files:
      1. Removes extraneous header cards.  In particular, removes:
         wcsname, wcsaxes, cdelt1, cdelt2, crota1, crota2, pc1_1, pc1_2, etc.
      2. Creates a new CD matrix, with the correct pixel scale and rotation,
         and writes it into the appropriate header cards.

   These steps are necessary because the WIRC raw files (at least those
   obtained in March 2012) have too many types of WCS header cards, and they
   are not all consistent.  This is confusing ds9, swarp, etc.

   Inputs:
      infiles   -  list of input files.  These files will be modified in place.
   """

   print ""
   for f in infiles:
      hdulist = im.open_fits(f,"update")
      hdr = hdulist[0].header

      """ Get useful header information """
      pixscale = hdr['cdelt2']*3600.
      rot = hdr['ccrot']

      """ Delete the extraneous headers if they exist """
      clist = ['wcsname', 'wcsaxes', 'pc1_1', 'pc1_2', 'pc2_1', 'pc2_2',
               'cdelt1', 'cdelt2', 'crota1', 'crota2']
      for i in clist:
         try:
            del hdr[i]
         except:
            print "No %s header card found for %s" % (i,f)

      """ Generate a CD matrix based on the inputs """
      cdmatx = coords.rscale_to_cdmatrix(pixscale, rot, verbose=False)

      """ Transfer new CD matrix to fits header cards and save """
      hdr.update('cd1_1',cdmatx[0,0])
      hdr.update('cd1_2',cdmatx[0,1])
      hdr.update('cd2_1',cdmatx[1,0])
      hdr.update('cd2_2',cdmatx[1,1])

      """ Write out result and exit """
      hdulist.flush()
      del hdulist, hdr
      print "Cleaned WCS information for %s" % f

#------------------------------------------------------------------------------

def run_astrom(listfile, astcat='usno-b1'):
   """
   Sets up for scamp and then runs it on the list.
   """

   """ Make the LDAC fits catalogs as inputs to scamp """
   os.system('rm tmptmp.cat')
   os.system('touch tmptmp.cat')
   for line in open(listfile,'r'):
      f = line.split()[0]
      tmpfile = f.replace('ff_wirc0','tmp')
      tmpcat = tmpfile.replace('fits','cat')
      os.system('echo %s >> tmptmp.cat' % tmpcat)
      ast.remove_nans(f,tmpfile)
      ast.make_fits_cat(tmpfile,gain=5.467,ncoadd=5,outcat=tmpcat)
   
   """ Run scamp on the list of SExtractor catalogs  """
   scamplist = '@tmptmp.cat'
   ast.run_scamp(scamplist,astcat)

   """ Rename the output *.head files produced by scamp """
   for line in open('tmptmp.cat','r'):
      f = line.split()[0]
      tmphead = f.replace('.cat','.head')
      outhead = tmphead.replace('tmp','ff_wirc0')
      os.rename(tmphead,outhead)

#------------------------------------------------------------------------------

def apply_astrom(listfile):
   """
   Once the solutions are acceptable, create new versions of the original
   input files (the flat-fielded ones) and copy the new wcs header information
   from the appropriate *.head file into the file
   """

   """ 
   Create wcs_*fits copies of the input file list and lists of wcs and
   ascii header files.
   """
   for line in open(listfile,'r'):
      ff = line.split()[0]
      wcsfile = ff.replace('ff','wcs_ff')
      headfile = ff.replace('fits','head')
      os.system('cp -v %s %s' % (ff,wcsfile))
      ast.import_ascii_header(wcsfile,headfile)
   
#------------------------------------------------------------------------------

def fixphot(inlist, photcat, edgedist=50., dpixmax=6.):
   """
   A porting of the photometric part of the wirc_fixwcs.pro procedure in 
   Code/IDL/wircsoft.
   The reason for doing this is that the astrometry part of the IDL procedure
   wasn't completely working, and so the astrometry is now being done
   with scamp instead.

   Inputs:
      inlist - list of data files
      astcat - photometric catalog (acquired from 2MASS) in
                an earlier part of the pipeline.
   """

   from scipy import stats

   """ Read astrometric and photometric information """
   astra0,astdec0,j0,h0,k0 = n.loadtxt(photcat,unpack=True)

   """ Loop on images """
   for f in inlist:

      """ Find the stars in the catalog that lie within the image """
      hdr = pf.getheader(f)
      nx = hdr['naxis1']
      ny = hdr['naxis2']

      mx0,my0 = wcsmwa.sky2pix(hdr,astra0,astdec0)
      goodmask = (mx0>edgedist) & (mx0<nx-edgedist) & \
          (my0>edgedist) & (my0<ny-edgedist)
      mra  = astra0[goodmask]
      mdec = astdec0[goodmask]
      mx   = mx0[goodmask]
      my   = my0[goodmask]
      j    = j0[goodmask]
      h    = h0[goodmask]
      k    = k0[goodmask]

      """ Determine the filter for the data """
      filter2 = hdr['aft']
      if filter2[0:2] == 'Ks':
         compmag = k
         print "Comparing %s vs. 2MASS Ks" % f
      elif filter2[0] == 'H':
         compmag = h
         print "Comparing %s vs. 2MASS H" % f
      elif filter2[0] == 'J':
         compmag = j
         print "Comparing %s vs. 2MASS J" % f
      else:
         compmag = k
         print "No filter info found.  Defaulting to K"

      """ 
      Run SExtractor on the image and read in positions, mag, fwhm, and flags
      """
      ast.remove_nans(f,'tmp_sext.fits')
      outcat = f.replace('fits','cat')
      ast.make_fits_cat('tmp_sext.fits',format='ascii',outcat=outcat,
                        logfile='sext_wcs.log')
      x,y,mag,fwhm,flag = n.loadtxt(outcat,usecols=(6,7,17,19,20),unpack=True)

      """
      Estimate FWHM from good detections in image
      """
      mask = (x>200.) & (x<1800.) & (y>200.) & (y<1800.) & (flag==0)
      if mask.sum() >0.:
         goodfwhm = 0.248 * n.median(fwhm[mask])
      else:
         goodfwhm = -99.

      """ 
      Do pixel-based distance calculations, which should be good enough 
      since the images have had scamp run on them.
      """
      magnum = mx.size
      matchmag = n.zeros(magnum)
      for ii in range(magnum):
         dist = n.sqrt((x-mx[ii])**2 + (y-my[ii])**2)
         dind = n.argsort(dist)
         if (dist[dind[0]])>dpixmax:
            matchmag[ii] = n.nan
         else:
            matchmag[ii] = compmag[ii] - mag[dind[0]]

      medzp = stats.stats.nanmedian(matchmag)
      hdulist = im.open_fits(f,'update')
      hdr = hdulist[0].header
      hdr.update('seeing',goodfwhm,'Median FWHM in arcsec')
      hdr.update('magnum',magnum,'Number of 2MASS stars in image')
      hdr.update('magzp',medzp,'Magnitude zeropoint relative to 2MASS')
      hdulist.flush()
      os.remove('tmp_sext.fits')
      print ""
      print "Average seeing estimate: %6.2f" % goodfwhm
      print "Number of 2MASS stars in image: %d" % magnum
      print "Median zero point: %+7.2f" % medzp
      print ""
      print "-----------------------------------------------------------------"
      print ""
