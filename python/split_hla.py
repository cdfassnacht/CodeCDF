"""
Code to take an input fits file that has been downloaded from the HST
Legacy Archve (HLA) and split it into a science fits file and a weight
fits file.

The HLA file is expected to have the following structure:

 No.    Name         Type      Cards   Dimensions   Format
 0    PRIMARY     PrimaryHDU     580   ()              
 1    SCI         ImageHDU       215   (630, 630)   float32   
 2    WHT         ImageHDU        71   (630, 630)   float32   
 3    CTX         ImageHDU        71   (630, 630)   int32   

so the science file is in extension 1 and the weight file is in extension 2

"""

import sys
from astropy.io import fits as pf

""" Check command line """
if len(sys.argv)<3:
    print ''
    print 'Usage:'
    print 'python %s [infile] [outfile_root]' % sys.argv[0]
    print ''
    exit()

""" Get the file names """
infile = sys.argv[1]
outroot = sys.argv[2]
outsci = '%s.fits' % outroot
outwht = '%s_wht.fits' % outroot

""" Open data file """
try:
    hdu = pf.open(infile)
except:
    print ''
    print 'ERROR: Could not open %s' % infile
    print ''
    exit()

""" Get the data and header """
scidat = hdu[1].data
whtdat = hdu[2].data
prihdr = hdu[0].header
scihdr = hdu[1].header
whthdr = hdu[2].header

""" List of important header cards """
prilist = ['telescop','instrume','targname','camera','focus','aperture',
           'exptime','date-obs','proposid'
           ]
scilist = ['wcsaxes','ctype1','ctype2','crpix1','crpix2','crval1','crval2'
           ]

""" Populate the output file headers """
scihdu = pf.PrimaryHDU(scidat)
whthdu = pf.PrimaryHDU(whtdat)
for i in prilist:
    try:
        scihdu[i] = prihdu[i]
    except:
        print 'Warning: %s keyword not found in %s' % (i.upper(),infile)
    try:
        

""" Create the output files """

""" Clean up and close """
hdu.close()
