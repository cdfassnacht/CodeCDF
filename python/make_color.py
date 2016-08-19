# make_color.py
# Takes 3 input images and makes a rgb output.

import sys,pyfits,scipy,pylab
from scipy import ndimage

file1 = sys.argv[1]
file2 = sys.argv[2]
file3 = sys.argv[3]

data1 = pyfits.open(file1)[0].data.astype(scipy.float32)
data2 = pyfits.open(file2)[0].data.astype(scipy.float32)
data3 = pyfits.open(file3)[0].data.astype(scipy.float32)

# Load weight files
file1 = file1.replace(".fits",".weight.fits")
file2 = file2.replace(".fits",".weight.fits")
file3 = file3.replace(".fits",".weight.fits")
w1 = pyfits.open(file1)[0].data.astype(scipy.float32)
w2 = pyfits.open(file2)[0].data.astype(scipy.float32)
w3 = pyfits.open(file3)[0].data.astype(scipy.float32)

# Mask of good points in ALL THREE images
w = w1*w2*w3

# The array bad masks non-image reasons. This is different than the weight
#   array w because w masks parts within the image interior that we would like
#   to keep in the output. The maximum_filter 'fills in' these bad pixels.
#   It looks like the wht map from multidrizzle or swarp is screwing up; I
#   don't know why these pixels are zero (a small number sure, but zero?!)....
w1 = ndimage.maximum_filter(w1,11)
w2 = ndimage.maximum_filter(w2,11)
w3 = ndimage.maximum_filter(w3,11)
bad = w1*w2*w3

cond = w!=0
cond2 = bad==0

# Perform statistics for scaling
tmp = data1[cond]
tmp.sort()
min = tmp[tmp.size*0.1]
max = tmp[tmp.size*0.995]
data1[data1<min] = min
data1[data1>max] = max
data1 -= min
data1 /= data1[cond].max()
data1[cond2] = 0.

tmp = data2[cond]
tmp.sort()
min = tmp[tmp.size*0.1]
max = tmp[tmp.size*0.995]
data2[data2<min] = min
data2[data2>max] = max
data2 -= min
data2 /= data2[cond].max()
data2[cond2] = 0.

tmp = data3[cond]
tmp.sort()
min = tmp[tmp.size*0.1]
max = tmp[tmp.size*0.995]
data3[data3<min] = min
data3[data3>max] = max
data3 -= min
data3 /= data3[cond].max()
data3[cond2] = 0.

# Create the RGB array
data = scipy.empty((data1.shape[0],data1.shape[1],3))
# Note that for some reason the Y-axis is flipped. The ::-1 syntax flips an
#   array about the axis.
data[:,:,0] = data1[::-1].copy()
data[:,:,1] = data3[::-1].copy()
data[:,:,2] = data2[::-1].copy()
pylab.imshow(data)
pylab.show()
