"""
 add_rough_wcs.py

 Usage: python add_rough_wcs.py [pixel_scale] [fits_file(s)]
  Inputs:
   pixel_scale  - a numerical value giving the rough pixel scale in arcsec/pix
   fits_file(s) - one or more fits files on which the program will be run
                  NOTE: Wild cards are fine
  Examples:
   python add_rough_wcs.py 0.211 lris0234.fits
   python add_rough_wcs.py 0.211 lris0234.fits lris0235.fits lris0236.fits
   python add_rough_wcs.py 0.211 lris02*.fits

 Description: Uses the telescope pointing information encoded in the fits
  header cards to add a rough WCS solution to a fits file that is missing
  WCS information.  The user passes the rough pixel scale (in arcsec/pix)
  as the first argument.  The remaining argument(s) is/are the name(s) of
  the fits file(s) to update.  Wildcards can be used.
"""

import pyfits,scipy
from pyfits import getval
import sys,os
import wcs

def update_wcs_from_pointing(files,pixscale):
   for f in files:
      if f[-5:]!=".fits":
         continue
      print "%s ....Updating values" % f
      fullfits = pyfits.open(f,mode='update')
      sci = fullfits[0]
      shape = sci.data.shape
      ra = wcs.ra2deg(sci.header['RA'].strip())
      dec = wcs.dec2deg(sci.header['DEC'].strip())
      sci.header.update('CRPIX1',shape[1]/2.)
      sci.header.update('CRVAL1',ra)
      sci.header.update('CRPIX2',shape[0]/2.)
      sci.header.update('CRVAL2',dec)
      sci.header.update('CD1_1',-pixscale/3600.)
      sci.header.update('CD1_2',0.)
      sci.header.update('CD2_1',0.)
      sci.header.update('CD2_2',pixscale/3600.)
      sci.header.update('CTYPE1','RA---TAN')
      sci.header.update('CTYPE2','DEC--TAN')
      sci.header.update('EQUINOX',2000.0)
      sci.header.update('RADESYS','FK5')
      del sci.header['CDELT1']
      del sci.header['CDELT2']
      fullfits.flush()
      print "  n   CTYPEn  CRPIXn    CRVALn       CDn_1        CDn_2"
      print " --- -------- ------- ---------- ------------- -------------"
      print "  1  %s %7.2f %10.6f %13.6e %13.6e" \
       % (getval(f,'CTYPE1'),getval(f,'CRPIX1'),getval(f,'CRVAL1'), \
       getval(f,'CD1_1'),getval(f,'CD1_2'))
      print "  2  %s %7.2f %10.6f %13.6e %13.6e" \
       % (getval(f,'CTYPE2'),getval(f,'CRPIX2'),getval(f,'CRVAL2'), \
       getval(f,'CD2_1'),getval(f,'CD2_2'))


"""------------------------Main program -----------------------------------"""
def main():
   if len(sys.argv)<3:
      print ""
      print "Usage: python add_rough_wcs.py [pixel_scale] [fits_file(s)]"
      print " Inputs:"
      print "  pixel_scale  - a numerical value giving the rough pixel scale in arcsec/pix"
      print "  fits_file(s) - one or more fits files on which the program will be run"
      print "                 NOTE: Wild cards are fine"
      print " Examples:"
      print "  python add_rough_wcs.py 0.211 lris0234.fits"
      print "  python add_rough_wcs.py 0.211 lris0234.fits lris0235.fits lris0236.fits"
      print "  python add_rough_wcs.py 0.211 lris02*.fits"
      print ""
      print "Description: Uses the telescope pointing information encoded in the fits"
      print " header cards to add a rough WCS solution to a fits file that is missing"
      print " WCS information.  The user passes the rough pixel scale (in arcsec/pix)"
      print " as the first argument.  The remaining argument(s) is/are the name(s) of"
      print " the fits file(s) to update.  Wildcards can be used."
      print ""
      return
   else:
      pixscale = float(sys.argv[1])
      files = sys.argv[2:]
   print ""
   update_wcs_from_pointing(files,pixscale)



if __name__ == "__main__":
   main()
