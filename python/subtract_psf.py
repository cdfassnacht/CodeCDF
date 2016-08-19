import os
import pyfits
import scipy
from scipy import ndimage,optimize

# Function poststamp - cuts out a postage stamp from a larger image
#
# Inputs:
#  data   - full image data array
#  cx     - x value of central pixel
#  cy     - y value of central pixel
#  csize  - length of one side of the postage stamp
# Output:
#  cutout - postage stamp data array

def poststamp(data,cx,cy,csize):
   oddeventest = csize / 2.0 - int(csize/2.0)
   if oddeventest==0:
      halfsize = int(csize/2.0)
   else:
      halfsize = int((csize+1)/2.0)

   # Make the cutout.  Remember y coordinates come first
   cutout = data[cy-halfsize:cy+halfsize,cx-halfsize:cx+halfsize].copy()
   return cutout

def eval_psf_match(p,data1,data2):
   amp,x,y = p
   coords = scipy.indices(data1.shape).astype(scipy.float64)
   coords[0] += y
   coords[1] += x
   data = amp * data1

   shift = ndimage.map_coordinates(data,coords,output=scipy.float64)

   return (shift - data2).flatten()

def find_shift(template,target):

   #
   # Initialize guess for amplitude and shift
   #

   p = [(template.max()/target.max()),0.,0.]
   
   #
   # Solve for the shift between the template and target
   #
   
   pfinal,ier = optimize.leastsq(eval_psf_match,p,(template,target))
   print pfinal,ier

   return pfinal

def shift_template(template,tempcore,target,targcore):

   # Get coordinates of each pixel in the full template array

   coords = scipy.indices(template.shape).astype(scipy.float64)

   # Find the shift between template and target, using just the central
   #  core region of each PSF

   pshift = find_shift(tempcore,targcore)
   
   #
   # Shift the template star to match the centering of the target and return
   #  the result
   #
   
   coords[1] += pshift[1]
   coords[0] += pshift[2]
   shiftstar = ndimage.map_coordinates(template,coords,output=scipy.float64)

   return pshift[0] * shiftstar

def putpost(data,cx,cy,cutout):
   csize = cutout.shape()[0]
   oddeventest = csize / 2.0 - int(csize/2.0)
   if oddeventest==0:
      halfsize = int(csize/2.0)
   else:
      halfsize = int((csize+1)/2.0)

   # Make the cutout.  Remember y coordinates come first
   cutout = data[cy-halfsize:cy+halfsize,cx-halfsize:cx+halfsize]
   return cutout

#***********************************************************************
#
# Main program

#
# Read in data file
#

data = pyfits.open("1520_nirc2_ao_K.fits")[0].data.astype(scipy.float64)
cutsize = 20
fullsize = 100
s1x = 349
s1y = 420
ax = 691
ay = 700
bx = 548
by = 634

#
# Define the boxes for star 1
#
# *** NB: x coordinate given first here ***

star1_core = poststamp(data,s1x,s1y,cutsize)
star1 = poststamp(data,s1x,s1y,fullsize)

#
# Define the boxes for quasar image A
#

qa_core = poststamp(data,ax,ay,cutsize)
qa = poststamp(data,ax,ay,fullsize)

#
# Define the boxes for quasar image B
#

qb_core = poststamp(data,bx,by,cutsize)
qb = poststamp(data,bx,by,fullsize)

#
# Create a version of star 1, shifted to match centering of quasar A
#

shift1_2_a = shift_template(star1,star1_core,qa,qa_core)

#
# Create a version of star 1, shifted to match centering of quasar B
#

shift1_2_b = shift_template(star1,star1_core,qb,qb_core)
shifta_2_b = shift_template(qa,qa_core,qb,qb_core)

#
# Make a copy of the original data (just in case we screw up)
#

newdata = data.copy()

#
# Subtract the shifted PSF from the original data and write output to
#  fits file
#

oddeventest = fullsize / 2.0 - int(fullsize/2.0)
if oddeventest==0:
   halfsize = int(fullsize/2.0)
else:
   halfsize = int((fullsize+1)/2.0)

newdata[ay-halfsize:ay+halfsize,ax-halfsize:ax+halfsize] -= shift1_2_a
newdata[by-halfsize:by+halfsize,bx-halfsize:bx+halfsize] -= shift1_2_b
#newdata = qb - shift1_2_b
pyfits.PrimaryHDU(newdata).writeto("test.fits")
