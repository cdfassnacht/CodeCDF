"""
Estimates the rms noise level in an image (stored as a fits file), given the
center and dimensions of a region to use for the estimation.  This region
should be representative of the sky level in the image and should not
contain any bright sources.

Usage: python get_rms.py [filename] [xcent] [ycent] [xsize] [ysize]

 filename - name of fits file containing the image data
 xcent    - x coordinate of the center of the region used to calculate the
            image statistics
 ycent    - y coordinate of the center of the region used to calculate the
            image statistics
 xsize    - size of the region in the x direction
 ysize    - size of the region in the y direction
"""

import sys
import numpy as np
from astropy.io import ascii
import imfuncs as imf

# """ Check the command line syntax """
# if len(sys.argv)<6:
#     print ''
#     print 'get_rms.py requires 5 input parameters.  Only %d were given' \
#         % (len(sys.argv) - 1)
#     print ''
#     print 'Usage: python get_rms.py [filename] [xcent] [ycent] [xsize] [ysize]'
#     print ''
#     print 'filename - name of fits file containing the image data'
#     print 'xcent    - x coordinate of the center of the region used to '
#     print '           calculate the image statistics'
#     print 'ycent    - y coordinate of the center of the region used to'
#     print '           calculate the image statistics'
#     print 'xsize    - size of the region in the x direction'
#     print 'ysize    - size of the region in the y direction'
#     print ''
#     exit()
# 
# """ Get the information from the command line """
# infile = sys.argv[1]
# xcent = float(sys.argv[2])
# ycent = float(sys.argv[3])
# xsize = float(sys.argv[4])
# ysize = float(sys.argv[5])
# 
# """ Get the rms from the selected region in the image"""
# rms = imf.get_rms(infile, xcent, ycent, xsize, ysize, verbose=False)
# print ''
# print '%s: In the requested region, RMS = %f' % (infile,rms)
# print ''

infile = sys.argv[1]
intab = ascii.read(infile, guess=False, format='commented_header')
rms = np.zeros(len(intab))
print '    File               RMS'
print '--------------- ---------------'
for i in range(len(intab)):
    trow = intab[i]
    f = '%s/%s.fits' % (trow['Folder'], trow['Filename'])
    rms[i] = imf.get_rms(f, trow['xcent'], trow['ycent'], trow['xsize'],
                         trow['ysize'], verbose=False)
    print '%-13s %f'


