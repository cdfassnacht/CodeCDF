"""
Code to adjust the WCS of one or more fits files by assigning the
requested RA and Dec (given on the command line) to a selected pixel in 
each fits file.  The pixel is selected interactively.

A typical usage would be in fixing a set of files where the WCS does
not quite agree.  An object that is in all of the images with know (RA,Dec)
coordinates is selected and the user interactively marks the position
of the star in each file.

Usage: python wcs_1click [ra] [dec] (-p [pixscale]) [fitsfile(s)]

             -- or --

   python wcs_1click.py [ra] [dec] (-p [pixscale]) -f [flatfile] [fitsfile(s)]

Input parameters:
    ra          - RA in decimal degrees
    dec         - Dec in decimal degrees
    pixscale    - [OPTIONAL] pixel scale in arcsec/pix
    fmax        - [OPTIONAL] maximum display value, in sigma above mean
                  Default value if this is not set is fmax=10
    flatfile    - [OPTIONAL] flat-field file to apply to the file(s) before
                   displaying it
    fitsfile(s) - Either a single filename or a wildcard expression
                   e.g.,  m13*fits
"""

import sys
import numpy as np
from astropy import wcs
from astropy.io import fits as pf
from matplotlib import pyplot as plt
from specim import imfuncs as imf

""" Check command line syntax """
if len(sys.argv)<4:
    print('')
    print('Usage:')
    print(' python wcs_1click.py [ra] [dec] (flag1 flag1val'
          ' flag2 flag2val...) [fitsfile(s)]')
    print('')
    print('Required inputs:')
    print(' ra  - RA in decimal degrees')
    print(' dec - Dec in decimal degrees')
    print(' fitsfile(s) - Either a single filename or a wildcard expression')
    print('  e.g.,  m13*fits')
    print('')
    print('OPTIONAL FLAGS and associated parameters')
    print('  -p [pixscale]    - pixel scale in arcsec/pix')
    print('  -flat [flatfile] - flat-field file to be applied to the input'
          ' fits files')
    print('  -fmax [fmax]     - maximum flux value, in sigma above mean,'
          ' for display')
    print('                     Default value: 10')
    print('')
    exit()

""" Set up variables for later use """
filestart = 3
pixscale = None
fmax = 10.
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
    elif sys.argv[filestart] == '-fmax':
        try:
            fmax = float(sys.argv[filestart+1])
        except ValueError:
            msg = 'ERROR: fmax is not a floating point number'
            no_error = False
        except IndexError:
            msg = 'ERROR: -fmax used but no fmax value is given'
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

"""
Set up the pixel scale to use
The default is to use the WCS information in the file header, but if
 the pixscale parameter has been set then its value overrides any
 pixel scale information in the header
"""
if pixscale is not None:
    pixscale /= 3600.

""" Loop through the input files, marking the object in each one """
for infile, crp1, crp2 in zip(files, crpix1, crpix2):
    """ Open and display the image """
    im1 = imf.Image(infile)
    if flat is not None:
        im1.data /= flat
    im1.zoomsize = subimsize
    im1.display(fmax=fmax, mode='xy', title=im1.infile)

    """ Run the interactive zooming and marking """
    im1.start_interactive()
    plt.show()

    """ Set the crpix values to the marked location """
    if im1.dispim.xmark is not None:
        crp1 = im1.dispim.xmark + 1
    if im1.dispim.ymark is not None:
        crp2 = im1.dispim.ymark + 1

    """
    If there is no WCS information in the input file, create a base version
    to be filled in later
    """
    if im1['input'].wcsinfo is None:
        im1['input'].wcsinfo = wcs.WCS(naxis=2)
        im1['input'].wcsinfo.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    im1['input'].update_crpix((crp1, crp2), verbose=False)
    im1['input'].update_crval((ra, dec), verbose=False)
    im1.wcsinfo = im1['input'].wcsinfo
    if flat is not None:
        im1.data *= flat
    im1.save(verbose=False)
    del(im1)

"""
Report on the updated values
"""
print('')
print('File                       CRVAL1      CRVAL2     CRPIX1   CRPIX2 ')
print('------------------------ ----------- ----------- -------- --------')
for infile in files:
    hdr = pf.getheader(infile)
    f = infile[:-5]
    print('%-24s %11.7f %+11.7f %8.2f %8.2f' % (f, ra, dec, hdr['crpix1'],
                                                hdr['crpix2']))
    del hdr

# """ 
# Second pass through the fits files, assigning the RA and Dec to the
# appropriate CRPIX, and setting the pixel scale if requested (pix scale
# not yet implemented)
# """
# if pixscale is not None:
#     pixscale /= 3600.
# print ''
# print 'File                       CRVAL1      CRVAL2     CRPIX1   CRPIX2 '
# print '------------------------ ----------- ----------- -------- --------'
# for i in range(len(files)):
#     hdu = pf.open(files[i], mode='update')
#     hdr = hdu[0].header
#     hdr['crval1'] = ra
#     hdr['crval2'] = dec
#     hdr['crpix1'] = crpix1[i]
#     hdr['crpix2'] = crpix2[i]
#     try:
#         foo = hdr['ctype1']
#     except KeyError:
#         hdr['ctype1'] = 'RA---TAN'
#     try:
#         foo = hdr['ctype2']
#     except KeyError:
#         hdr['ctype2'] = 'DEC--TAN'
#     hdu.flush()
#     f = files[i][:-5]
#     print '%-24s %11.7f %+11.7f %8.2f %8.2f' % (f,ra,dec,crpix1[i],crpix2[i])

