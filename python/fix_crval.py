"""
 fix_crval.py

 Usage: python fix_crval.py [delta_ra] [delta_dec] [fits_files]
  Inputs:
   delta_ra     - offset to apply to RA, in arcsec
   delta_dec    - offset to apply to Dec, in arcsec
   fits_file(s) - one or more fits files on which the program will be run
                  NOTE: Wild cards are fine
  Examples:
   python fix_crval.py 680.1 71.4 0.211 lris0234.fits
   python fix_crval.py 680.1 71.4 lris0234.fits lris0235.fits lris0236.fits
   python fix_crval.py 680.1 71.4 lris02*.fits

 Description: Updates the CRVAL header cards for one or more fits files.
  The users passes the offsets, in ARCSEC, to be applied to the current CRVAL 
  values as the first two arguments.  
  The remaining argument(s) is/are the name(s) of the fits file(s) to update.
  Wildcards can be used.
"""

import pyfits,scipy
from pyfits import getval
from math import cos,pi
import sys,os
import wcs

def fix_wcs_crval(files,delta_ra,delta_dec):
   for f in files:
      if f[-5:]!=".fits":
         continue
      print "%s ....Updating CRVAL values" % f
      fullfits = pyfits.open(f,mode='update')
      sci = fullfits[0]
      ra  = getval(f,'CRVAL1')
      dec = getval(f,'CRVAL2')
      ranew  = ra - delta_ra/(3600.0*cos(dec*pi/180.0))
      decnew = dec - delta_dec/3600.0
      sci.header.update('CRVAL1',ranew)
      sci.header.update('CRVAL2',decnew)
      print " %f %f --> %f %f" % (ra,dec,ranew,decnew)
      fullfits.flush()


"""------------------------Main program -----------------------------------"""
def main():
   if len(sys.argv)<4:
      print ""
      print "fix_crval.py"
      print ""
      print "Usage: python fix_crval.py [delta_ra] [delta_dec] [fits_files]"
      print " Inputs:"
      print "  delta_ra     - offset to apply to RA, in arcsec"
      print "  delta_dec    - offset to apply to Dec, in arcsec"
      print "  fits_file(s) - one or more fits files on which the program will"
      print "                  be run. NOTE: Wild cards are fine"
      print " Examples:"
      print "  python fix_crval.py 680.1 71.4 lris0234.fits"
      print "  python fix_crval.py 680.1 71.4 lris0234.fits lris0235.fits"
      print "  python fix_crval.py 680.1 71.4 lris02*.fits"
      print ""
      print "Description: Updates the CRVAL header cards for one or more fits"
      print " files.  The users passes the offsets, in ARCSEC, to be applied"
      print " to the current CRVAL values as the first two arguments."
      print " The remaining argument(s) is/are the name(s) of the fits file(s)"
      print "  to update.  Wildcards can be used."
      print ""
      return
   else:
      dalpha = float(sys.argv[1])
      ddelta = float(sys.argv[2])
      files = sys.argv[3:]
   print ""
   fix_wcs_crval(files,dalpha,ddelta)



if __name__ == "__main__":
   main()
