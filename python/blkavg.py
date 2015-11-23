"""
blkavg.py

Description: This bit of code replicates (more or less) the blkavg task
 within iraf/pyraf.  The purpose is to take an input fits file and to
 create an output fits file that is smaller by an integer factor (N).
 The code takes NxN blocks of pixels from the input file and creates 1
 pixel in the output file.  Therefore, unlike ndimage.zoom or 
 ndimage.map_coordinates, there is no interpolation and therefore no
 introduction of correlated noise between the pixels.

NOTE: This code has taken the resamp function directly from Matt
Auger's indexTricks library.  Right now there is NO ERROR CHECKING
in that part of the code.  User beware!

Usage: python blkavg.py [input_fits] [output_fits] [blkavg_factor]
 Example:  python blkavg.py big_file.fits ave_file.fits 3
"""

import sys
import numpy as n
try:
    from astropy.io import fits as pf
except:
    try:
        import pyfits
    except:
        print ''
        print 'ERROR. Neither astropy.io.fits nor pyfits found'
        print ''
        exit()

#---------------------------------------------------------------------------

def resamp(a,factor,add=False):
    """
    Code to do the block averaging, taken from Matt Auger's indexTricks library

    Inputs:
       a      - input data array
       factor - block averaging / summing factor
       add    - set to True if desired output is sum rather than average
                default = False

    Output:
       o      - block-averaged or block-summed array
    """
    arr = a.copy()
    """ 
    Cut off rows and columns to get an integer multiple of the factor
    in each dimension
    """
    dx = arr.shape[1] % factor
    dy = arr.shape[0] % factor
    if dx>0:
        arr = arr[:,:-dx]
    if dy>0:
        arr = arr[:-dy,:]
    x = arr.shape[1]/factor
    y = arr.shape[0]/factor
    o = n.zeros((y,x))
    for i in range(factor):
        for j in range(factor):
            o += arr[i::factor,j::factor]
    if add==True:
        return o
    return o/factor**2

#-----------------------------------------------------------------------

""" Main program """

""" Check command-line format """
if len(sys.argv)<4:
    print ''
    print 'Usage: python blkavg.py [input_fits] [output_fits] [blkavg_factor]'
    print '  Example:  python blkavg.py big_file.fits ave_file.fits 3'
    print ''
    exit()

"""
Get information from the command line
"""

infits  = sys.argv[1]
outfits = sys.argv[2]
try:
    blkfact = int(sys.argv[3])
except:
    print ''
    print 'ERROR. Third command-line input cannot be converted to an integer'
    print ''
    exit()

""" 
Do the averaging.
NOTE: For now assume data are in HDU 0
"""

try:
    (indat,hdr) = pf.getdata(infits,header=True)
except:
    print ''
    print 'ERROR.  Could not read data from input file %s' % infits
    print ''
    exit()
outdat = resamp(indat,blkfact)

""" Save the output file, clean up, and exit """
pf.PrimaryHDU(outdat,header=hdr).writeto(outfits)
del indat,outdat,hdr
