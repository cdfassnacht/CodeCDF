"""
Does exactly what the name says, adds a pixel scale (given in arcsec/pix) to
a fits file.  This code does _not_ change the crpix or crval -- only the
pixel scale.
"""

import sys
import glob
from astropy.io import fits as pf

""" Check the command line syntax """
if len(sys.argv) < 3:
    print('')
    print('Usage: python add_pixscale.py [pixscale] '
          '[file / file list / wildcard]')
    print('')
    print('Examples:')
    print('  python add_pixscale 0.211 myobj.fits')
    print('  python add_pixscale 0.211 my*fits')
    print('')
    print('NOTE: the pixel scale should be in arcsec/pix')
    print('')
    exit()

""" Get the information from the command line """
try:
    pixscale = float(sys.argv[1]) / 3600.
except:
    print('')
    print('ERROR: Expected the first passed parameter to be a number')
    print('')

infiles = sys.argv[2:]

""" Add the pixel scale to the input files """
print('')
for i in infiles:
    hdu = pf.open(i, mode='update')
    hdr = hdu[0].header
    """ Delete any old CD matrix that may exist """
    for j in ['cd1_1', 'cd1_2', 'cd2_1', 'cd2_2']:
        if j in hdr.keys():
            del hdr[j]
    """ Add in the new pixel information """
    hdr['cdelt1'] = -1.*pixscale
    hdr['cdelt2'] = pixscale
    hdu.flush()
    print('Added pixel scale information to %s' % i)
print('')



