import sys
import numpy as n
from astropy.io import fits as pf
from scipy.ndimage.filters import minimum_filter
import glob
from math import sqrt as msqrt
import imfuncs as imf

""" Get command line parameters """
medfile = sys.argv[1]
resampdir = sys.argv[2]
outfile = sys.argv[3]

# list of input files
# inlist = glob.glob("check/check_??.fits")
inlist = glob.glob("%s/resamp_??.fits" % resampdir)
print ''
print 'Input files'
print '-------------------'
for i in inlist:
    print i


# check that files (images) are same size
hdr0 = pf.getheader(inlist[0])
x0 = hdr0["naxis1"]
y0 = hdr0["naxis2"]
print ''
print 'Check that images are all the same size: %dx%d' % (x0, y0)
print '-------------------------------------------------------'

for i in inlist:
    hdr = pf.getheader(i)
    x = hdr["naxis1"]
    y = hdr["naxis2"]
    if x != x0 or y != y0:
        print ""
        print "Error: shape of %s does not match shape of %s"  % (i, inlist[0])
        print ""
        del hdr
        exit()
    else:
        print "%s: Shape matches" % i
        del hdr
del hdr0

""" Get the median file data """
med = imf.Image(medfile, verbose=False)
meddata = med.hdu[0].data
med.sigma_clip()
noisemed = med.rms_clip

""" initialize 3d arrays with zeros """
insci  = n.zeros((len(inlist), y0, x0))
wht    = n.zeros((len(inlist), y0, x0))

""" Find and mask the bad pixels, and the area around them """
print ''
print 'Masking the bad pixels...'
print '-------------------------'
for i in range(len(inlist)):
    print '  %s' % inlist[i]
    '''Put data in arrays'''
    infile = inlist[i]
    whtfile = infile.replace('.fits', '_wht.fits')
    imi = imf.Image(infile, verbose=False)
    insci[i, :, :] = imi.hdu[0].data.copy()
    insci_i = insci[i, :, :].copy()
    whti = imf.Image(whtfile, verbose=False)
    wht[i, :, :] = whti.hdu[0].data.copy()
    tmpwht = n.ones(whti.hdu[0].data.shape)

    '''Mask pixels associated with weight of zero'''
    mask = whti.hdu[0].data != 0
    newwhti = np.zeros(y0, x0)
    newwhti[mask] = 1

    """
    Find bad pixels
    NOTE: right now this does not account for Poisson noise where there are
    real objects.  This needs to be fixed!
    """
    imi.sigma_clip(mask=mask)
    noisei = imi.rms_clip
    noise = msqrt(noisei**2 + noisemed**2)
    diffi = (insci[i, :, :] - meddata) * newwhti
    mask2 = diffi > (3 * noise)

    ''' Grow the masked region around the bad pixels '''
    tmpwht[mask2] = 0
    tmp2 = minimum_filter(tmpwht, size=3)
    # outfile = "tmp%d.fits" %i
    wht[i, :, :] *= tmp2
    # newwht[i,:,:] *= tmp2

""" Do the weighted sum and also create the sum of the weights """
print ''
print 'Creating the weighted average image'
datasum = (insci * wht).sum(axis=0)  # 2d array
whtsum = wht.sum(axis=0)  # 2d array

""" Avoid devision by zero errors """
mask = whtsum == 0
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
