"""
Code to adjust the WCS of one or more fits files by assigning the
requested RA and Dec (given on the command line) to a selected pixel in 
each fits file.  The pixel is selected interactively.

A typical usage would be in fixing a set of files where the WCS does
not quite agree.  An object that is in all of the images with know (RA,Dec)
coordinates is selected and the user interactively marks the position
of the star in each file.

Usage: python apply_crval_interactive [ra] [dec] [fitsfile(s)]

Inputs:
    1. ra  - RA in decimal degrees
    2. dec - Dec in decimal degrees
    3. fitsfile(s) - Either a single filename or a wildcard expression
       e.g.,  m13*fits
"""

try:
    from astropy.io import fits as pf
except:
    try:
        import pyfits as pf
    except:
        print ''
        print 'ERROR: Could not import pyfits or astropy.io.fits'
        print ''
        exit()
from matplotlib import pyplot as plt
import numpy as n
import imfuncs as imf
import sys

""" Check command line syntax """
if len(sys.argv)<4:
    print ''
    print 'Usage: python apply_crval_interactive.py [ra] [dec] [fitsfile(s)]'
    print ''
    print 'Inputs:'
    print '  1. ra  - RA in decimal degrees'
    print '  2. dec - Dec in decimal degrees'
    print '  3. fitsfile(s) - Either a single filename or a wildcard expression'
    print '     e.g.,  m13*fits'
    print ''
    exit()

""" Set up variables for later use """
ra  = float(sys.argv[1])
dec = float(sys.argv[2])
if len(sys.argv)>4:
   files = sys.argv[3:]
else:
    files = [sys.argv[3],]
subimsize = 21
crpix1 = n.zeros(len(files))
crpix2 = n.zeros(len(files))

""" Loop through the input files, marking the object in each one """
for i in range(len(files)):
    im1 = imf.Image(files[i])
    im1.zoomsize = subimsize
    im1.display(sighigh=30.)
    im1.keypress_info()
    im1.start_interactive()
    plt.show()
    crpix1[i] = im1.xmark
    crpix2[i] = im1.ymark

""" 
Second pass through the fits files, assigning the RA and Dec to the
appropriate CRPIX
"""
print ''
print 'File                        CRVAL1      CRVAL2     CRPIX1   CRPIX2 '
print '------------------------- ----------- ----------- -------- --------'
for i in range(len(files)):
    hdu = pf.open(files[i],mode='update')
    hdr = hdu[0].header
    hdr['crval1'] = ra
    hdr['crval2'] = dec
    hdr['crpix1'] = crpix1[i]
    hdr['crpix2'] = crpix2[i]
    hdu.flush()
    f = files[i][:-5]
    print '%-24s %11.7f %+11.7f %8.2f %8.2f' % (f,ra,dec,crpix1[i],crpix2[i])
