"""
Code to adjust the WCS of one or more fits files by assigning the
requested RA and Dec (given on the command line) to a selected pixel in 
each fits file.  The pixel is selected interactively.

A typical usage would be in fixing a set of files where the WCS does
not quite agree.  An object that is in all of the images with know (RA,Dec)
coordinates is selected and the user interactively marks the position
of the star in each file.

Usage: python wcs_1click [ra] [dec] [fitsfile(s)]

Inputs:
    1. ra  - RA in decimal degrees
    2. dec - Dec in decimal degrees
    3. fitsfile(s) - Either a single filename or a wildcard expression
       e.g.,  m13*fits
"""

import sys
import numpy as np
try:
    from astropy.io import fits as pf
except ImportError:
    try:
        import pyfits as pf
    except ImportError:
        print ''
        print 'ERROR: Could not import pyfits or astropy.io.fits'
        print ''
        exit()
from matplotlib import pyplot as plt
try:
    from SpecIm import imfuncs as imf
except ImportError:
    import imfuncs as imf

""" Check command line syntax """
if len(sys.argv)<4:
    print('')
    print('Usage:')
    print(' python wcs_1click.py [ra] [dec] (-p [pixscale]) [fitsfile(s)]')
    print('              --- or ---')
    print(' python wcs_1click.py [ra] [dec] (-p [pixscale]) -f [flatfile]'
          ' [fitsfile(s)]')
    print('')
    print('Inputs:')
    print(' ra  - RA in decimal degrees')
    print(' dec - Dec in decimal degrees')
    print(' pixscale [OPTIONAL] - pixel scale in arcsec/pix')
    print(' flatfile [OPTIONAL] - flat-field file to be applied to the input'
          ' fits files')
    print(' fitsfile(s) - Either a single filename or a wildcard expression')
    print('  e.g.,  m13*fits')
    print('')
    exit()

""" Set up variables for later use """
filestart = 3
pixscale = None
flat = None
flatfile = None
start_files = False
subimsize = 21
no_error = True

""" Parse the command line """
ra  = float(sys.argv[1])
dec = float(sys.argv[2])
while start_files is False and no_error:
    if sys.argv[filestart] == '-p':
        try:
            pixscale = float(sys.argv[filestart+1])
        except ValueError:
            msg = 'ERROR: pixel scale is not a floating point number'
            no_error = False
        except IndexError:
            msg = 'ERROR: -p used but no pixel scale given'
            no_error = False
        filestart += 2
    elif sys.argv[filestart] == '-flat':
        try:
            flatfile = sys.argv[filestart+1]
        except IndexError:
            msg = 'ERROR: -flat used but no flat-field file is given'
            no_error = False
        filestart += 2
    else:
        start_files = True

if no_error is not True:
    print('')
    print('%s' % msg)
    print('')
    exit()

""" Create the input file list """
if len(sys.argv) > filestart + 1:
   files = sys.argv[filestart:]
else:
    files = [sys.argv[filestart],]
crpix1 = np.zeros(len(files))
crpix2 = np.zeros(len(files))

""" Read in the flat-field data """
if flatfile is not None:
    flat = pf.getdata(flatfile)
    print('')
    print('Using flat-field file: %s' % flatfile)

""" Loop through the input files, marking the object in each one """
for i in range(len(files)):
    im1 = imf.Image(files[i])
    if flat is not None:
        im1.hdu[0].data /= flat
    im1.zoomsize = subimsize
    im1.display(fmax=10., mode='xy', title=im1.infile)
    im1.start_interactive()
    plt.show()
    if im1.xmark is not None:
        crpix1[i] = im1.xmark
    if im1.ymark is not None:
        crpix2[i] = im1.ymark
    im1.close()
    del(im1)

""" 
Second pass through the fits files, assigning the RA and Dec to the
appropriate CRPIX
"""
print ''
print 'File                        CRVAL1      CRVAL2     CRPIX1   CRPIX2 '
print '------------------------- ----------- ----------- -------- --------'
for i in range(len(files)):
    hdu = pf.open(files[i], mode='update')
    hdr = hdu[0].header
    hdr['crval1'] = ra
    hdr['crval2'] = dec
    hdr['crpix1'] = crpix1[i]
    hdr['crpix2'] = crpix2[i]
    try:
        foo = hdr['ctype1']
    except:
        hdr['ctype1'] = 'RA---TAN'
    try:
        foo = hdr['ctype2']
    except:
        hdr['ctype2'] = 'DEC--TAN'
    hdu.flush()
    f = files[i][:-5]
    print '%-24s %11.7f %+11.7f %8.2f %8.2f' % (f,ra,dec,crpix1[i],crpix2[i])

