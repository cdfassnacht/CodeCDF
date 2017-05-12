#!/usr/bin/env/python

import sys
import numpy as n
from astropy.io import fits as pf
from scipy.ndimage.filters import minimum_filter
import glob
import ccdredux as ccd
from math import sqrt as msqrt

""" Get command line parameters """
medfile = sys.argv[1]
resampdir = sys.argv[2]
outfile = sys.argv[3]

# list of input files
# inlist=glob.glob("check/check_??.fits")
inlist = glob.glob("%s/resamp_??.fits" % resampdir)
print ''
print 'Input files'
print '-------------------'
for i in inlist:
    print i


#check that files (images) are same size
hdr0 = pf.getheader(inlist[0])
x0 = hdr0["naxis1"]
y0 = hdr0["naxis2"]
print ''
print 'Check that images are all the same size: %dx%d' % (x0, y0)
print '-------------------------------------------------------'

for i in inlist:
    hdr=pf.getheader(i)
    x=hdr["naxis1"]
    y=hdr["naxis2"]
    if x != x0 or y != y0:
        print ""
        print "Error: shape of %s does not match shape of %s" %(i, inlist[0]) 
        print ""
        exit()
    else: 
        print "%s: Shape matches" %i


#get the median files
meddata = pf.getdata(medfile)
#print meddata.shape
meanmed, noisemed = ccd.sigma_clip(meddata)
#print meanmed, noisemed

""" initialize 3d arrays with zeros """
check  = n.zeros((len(inlist), y0, x0))
wht    = n.zeros((len(inlist), y0, x0))
newwht = n.zeros((len(inlist), y0, x0))

""" Find and mask the bad pixels, and the area around them """
print ''
print 'Masking the bad pixels...'
for i in range(len(inlist)):
    '''Put data in arrays'''
    infile = inlist[i]
    check[i, :, :] = pf.getdata(infile)
    checki = check[i, :, :].copy()
    tmp = pf.getdata(infile.replace(".fits", "_wht.fits"))
    tmpwht = n.ones(tmp.shape)
    wht[i, :, :] = tmp.copy()

    '''Mask pixels associated with weight of zero'''
    mask = tmp !=0
    newwht[i, :, :][mask] = 1
    newwhti = newwht[i, :, :].copy()
    data = checki[mask]

    """ 
    Find bad pixels 
    NOTE: right now this does not account for Poisson noise where there are
    real objects.  This needs to be fixed!
    """
    meani, noisei = ccd.sigma_clip(data)
    noise = msqrt(noisei**2 + noisemed**2)
    diffi = (checki - meddata) * newwhti
    mask2 = diffi > (3*noise)

    ''' Grow the masked region around the bad pixels '''
    tmpwht[mask2] = 0
    tmp2 = minimum_filter(tmpwht,size=3)
    # outfile="tmp%d.fits" %i
    wht[i,:,:] *= tmp2
    newwht[i,:,:] *= tmp2

""" Do the weighted sum and also create the sum of the weights """
print ''
print 'Creating the weighted average image (may take a while)'
datasum = (check*wht).sum(axis=0) #2d array
whtsum = wht.sum(axis=0) #2d array

""" Avoid devision by zero errors """
mask = whtsum==0 
tmpwht = whtsum.copy()
tmpwht[mask] = 1 
datasum[mask] = 0

"""' Finally, define the weighted average """
whtedav = datasum / tmpwht
print '  ...Done'

""" Save the outputs """
print 'Writing output to %s' % outfile
pf.PrimaryHDU(whtedav).writeto(outfile, clobber=True)
outwht = outfile.replace('.fits', '_wht.fits')
pf.PrimaryHDU(whtsum).writeto(outwht, clobber=True)
print ''
