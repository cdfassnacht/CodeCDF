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

# list of input files
# inlist=glob.glob("check/check_??.fits")
inlist=glob.glob("%/resamp_??.fits" % resampdir)
print inlist


#check that files (images) are same size
hdr0=pf.getheader(inlist[0])
x0=hdr0["naxis1"]
y0=hdr0["naxis2"]
print x0,y0

for i in inlist:
    hdr=pf.getheader(i)
    x=hdr["naxis1"]
    y=hdr["naxis2"]
    if x != x0 or y != y0:
        print ""
        print "Error: shape of %s does not match shape of %s" %(i, inlist[0]) 
        print ""
    else: 
        print "%s matches" %i


#get the median files
meddata=pf.getdata(medfile)
#print meddata.shape
meanmed, noisemed = ccd.sigma_clip(meddata)
#print meanmed, noisemed

#initialize 3d arrays with zeros
check  = n.zeros((len(inlist), y0, x0))
wht    = n.zeros((len(inlist), y0, x0))
newwht = n.zeros((len(inlist), y0, x0))

#Find and mask the bad pixels, and the area around them
for i in range(len(inlist)):
    '''Put data in arrays'''
    infile=inlist[i]
    check[i,:,:]=pf.getdata(infile)
    checki = check[i,:,:].copy()
    tmp=pf.getdata(infile.replace(".fits", "_wht.fits"))
    tmpwht = n.ones(tmp.shape)
    wht[i,:,:]=tmp.copy()
    '''Mask pixels associated with weight of zero'''
    mask = tmp !=0
    newwht[i,:,:][mask] = 1
    newwhti=newwht[i,:,:].copy()
    data=checki[mask]
    '''Find bad pixels '''
    meani,noisei=ccd.sigma_clip(data)
    noise= msqrt(noisei**2 + noisemed**2)
    diffi=(checki - meddata)*newwhti
    mask2 = diffi > (3*noise)
    ''' Grow the masked region around the bad pixels '''
    tmpwht[mask2] = 0
    tmp2 = minimum_filter(tmpwht,size=3)
    outfile="tmp%d.fits" %i
    wht[i,:,:] *= tmp2
    newwht[i,:,:] *= tmp2

#define weighted data sum and sum of weights
datasum = (check*wht).sum(axis=0) #2d array
whtsum = wht.sum(axis=0) #2d array

#avoid devision by zero errors
mask=whtsum==0 
whtsum[mask]=1 
datasum[mask]=0

#define weighted average
whtedav = datasum/whtsum

#save to fits file
pf.PrimaryHDU(whtedav).writeto("whtedav234.fits", clobber=True)
