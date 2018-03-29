"""
Takes a SDSS catalog that is stored in a binary fits table format and converts
it to a ds9 region file.  This type of catalog is produced by selecting the
FITS output option from the SDSS imaging query web site


NOTE: This should really incorporate the code in catfuncs.py, which already
does this.

"""

import sys
import numpy as np
from astropy.table import Table
import catfuncs as cf

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

""" Read in the data """
scat = cf.Secat(infile, catformat='sdssfits')
tdat = scat.data

""" Get the ra and dec information from the catalog """
ragal  = tdat['ra'][scat.galmask]
decgal = tdat['dec'][scat.galmask]
if scat.starmask is not None:
   rastar  = tdat['ra'][scat.starmask]
   decstar = tdat['dec'][scat.starmask]

""" Write the output region file """
f = open(outfile,'w')
f.write('global color=green\n')
for i in range(ragal.size):
   f.write('fk5;circle(%f,%f,%.1f")\n' % (ragal[i], decgal[i], rcirc))
if scat.starmask is not None:
   f.write('global color=red\n')
   for i in range(rastar.size):
      f.write('fk5;box(%f,%f,%.1f",%.1f",45d)\n' % \
                 (rastar[i], decstar[i], 1.4*rcirc, 1.4*rcirc))
f.close()

""" Clean up and exit """
del(ragal, decgal)
if scat.starmask is not None:
   del(rastar, decstar)

#del radec
