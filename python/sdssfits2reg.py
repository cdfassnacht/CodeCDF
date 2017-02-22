"""
Takes a SDSS catalog that is stored in a binary fits table format and converts
it to a ds9 region file.  This type of catalog is produced by selecting the
FITS output option from the SDSS imaging query web site
"""

import astropy
from astropy import units as u
if astropy.__version__[:3] == '0.3':
   from astropy.coordinates import ICRS as SkyCoord
else:
   from astropy.coordinates import SkyCoord
try:
    from astropy.io import fits as pf
except:
    import pyfits as pf
import numpy as np
import sys

""" Test the command line for proper number of arguments """
if len(sys.argv)<4:
    print ''
    print 'USAGE: python sdssfits2reg.py [input_file] [output_file] [rcirc]'
    print ''
    exit()

""" Get info from the command line """
infile  = sys.argv[1]
outfile = sys.argv[2]
rcirc   = float(sys.argv[3])

""" Get the table data """
try:
    hdu = pf.open(infile)
except:
    print ''
    print 'ERROR: Could not read input file %s' % infile
    print ''
    exit()
tdat = hdu[1].data

""" Get the ra and dec information from the catalog """
ra  = tdat['ra']
dec = tdat['dec']
#radec = SkyCoord(gra,gdec,unit=(u.deg,u.deg))

""" Write the output region file """
f = open(outfile,'w')
f.write('global color=green\n')
for i in range(ra.size):
    f.write('fk5;circle(%f,%f,%.1f")\n' % (ra[i],dec[i],rcirc))
f.close()

""" Clean up and exit """
hdu.close()
del ra,dec
#del radec
