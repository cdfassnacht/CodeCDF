"""
This program is set up to do a very basic reduction of FORS1 imaging
 data.
"""

import pyfits,scipy
import numpy
import sys,os
import wcs
from scipy import stats,ndimage

#######################################################################

def sigma_clip(data,tol=3.):
   d = data.flatten()
   avg = d.mean()
   std = d.std()

   delta = 1
   while delta:
      size = d.size
      d = d[abs(d-avg)<tol*std]
      avg = d.mean()
      std = d.std()
      delta = size-d.size
   return avg,std

#######################################################################

#def make_fringe(frames,bias,flat):
#   y = slice(150,1950,1)
#   x = slice(150,2020,1)
#   fstaty = slice(1500,2500,1)
#   fstatx = slice(1500,2500,1)
#
#   fringe = scipy.empty(bias.shape)
#   tmp = []
#   band = pyfits.open(frames[0])[0].header['DWFILNAM'].strip()
#   print "Creating fringe frame for filter = %s.  Input frames:" % band
#   for f in frames:
#      print "   %s" % f
#      t = (pyfits.open(f)[0].data[y,x].astype(scipy.float32)-bias)/flat
#      t -= stats.stats.median(t[fstaty,fstatx],axis=None)
#      tmp.append(t.copy())
#   #
#   # Now assume that the observer has taken enough exposures
#   #  so that the median is actually a decent estimator of the
#   #  median that would be obtained after bad pixels were rejected
#   #
#   tmp2 = scipy.empty((len(tmp),tmp[0].shape[0],tmp[0].shape[1]))
#   print " Co-adding frames..."
#   for i in range(tmp2.shape[0]):
#      tmp2[i] = tmp[i].copy()
#   tmp3 = stats.stats.median(tmp2,axis=0)
#   #
#   # Now do a small-scale smoothing of the co-added frame
#   #
#   print " Smoothing the co-added frame"
#   fringe = ndimage.median_filter(tmp3,11)
#   pyfits.PrimaryHDU(fringe).writeto('fringe_%s.fits' % band)
#
#   return fringe

#######################################################################

def subtract_overscan(data,x,y):

   """This function finds the median values in each of the four overscan
      regions and subtracts them from the appropriate regions of the
      input data file.  It then converts the results back to electrons
      rather than ADU"""

   # Define bias region limits
   bx1 = slice(0,15,1)
   bx2 = slice(2065,2080,1)
   y1 = slice(0,1024,1)
   y2 = slice(1024,2048,1)

   # Define limits of regions associated with the four amps
   x1 = slice(16,1040)
   x2 = slice(1040,2064)

   # Define median values of overscan regions from appropriate data regions
   newdata = data.astype(scipy.float32)
   overscan = scipy.zeros((4,1))
   overscan[0] = scipy.median(newdata[y1,bx1].ravel())
   overscan[1] = scipy.median(newdata[y2,bx1].ravel())
   overscan[2] = scipy.median(newdata[y1,bx2].ravel())
   overscan[3] = scipy.median(newdata[y2,bx2].ravel())

   # Subtract overscan
   newdata[y1,x1] = newdata[y1,x1] - overscan[0]
   newdata[y2,x1] = newdata[y2,x1] - overscan[1]
   newdata[y1,x2] = newdata[y1,x2] - overscan[2]
   newdata[y2,x2] = newdata[y2,x2] - overscan[3]

   newdata = newdata[y,x]
   return newdata

#######################################################################

def gain_correct(filename,data):

   # Get the conversion from ADU to electrons.  Apparently, the chips
   #  are arranged on the detector as:
   #
   #     C   D
   #     A   B
   #
   # Hmm, my trial and error suggest that C and B need to be switched to
   #  bring the chip levels into (rough) agreement, except for the very
   #  first FORS I-band image of HE 0230, from Oct 1999.

   gain = scipy.zeros((4,1))

   gain[0] = pyfits.getval(filename,'HIERARCH ESO DET OUT1 GAIN')
   gain[1] = pyfits.getval(filename,'HIERARCH ESO DET OUT2 GAIN')
   gain[2] = pyfits.getval(filename,'HIERARCH ESO DET OUT3 GAIN')
   gain[3] = pyfits.getval(filename,'HIERARCH ESO DET OUT4 GAIN')

   for i in range(4):
      if gain[i] == 0:
         print "No gain found for detector %d" % (i+1)
         print "Setting gain to 1"
         gain[i] = 1.
   print "    gains: %f %f %f %f" % (gain[0], gain[1], gain[2], gain[3])

   # Define amplifier region limits
   x1 = slice(0,1024,1)
   x2 = slice(1024,2048,1)
   y1 = slice(0,1024,1)
   y2 = slice(1024,2048,1)

   # Correct data for the gains
   data[y1,x1] /= gain[0]
   data[y1,x2] /= gain[2]
   data[y2,x1] /= gain[1]
   data[y2,x2] /= gain[3]

   return data

#######################################################################

def forsim_redux(object,filter,inbias=None,inflat=None,insky=None,infringe=None):
   dir = os.curdir
   files = os.listdir(dir)

   _bias = []
   _flat = []
   _sky = []
   _sci = []

   # Define good part of data
   # Note that the raw imaging data (FORS1 in 2001, at least) has
   #  non-data (overscan??) regions both to the left and right of
   #  the main data block.  The data are produced by 4 amps, and
   #  it looks as if each overscan region contains the info on
   #  2 of the amps.  Thus, the LHS overscan region is split into
   #  two regions vertically.  The top one probably corresponds to
   #  the top left amp, while the bottom one probably corresponds
   #  to the bottom left amp.
   # The overscan regions are 16 pix wide. 1-15 and 2066-2080
   # The vertical cut in the overscan region occurs btwn y=1024 and 
   #  y=1025
   y = slice(0,2048,1)
   #x = slice(16,2064,1)
   x = slice(16,2050,1)

   # Define region of image to be used to calculate statistics
   staty = slice(700,1300,1)
   statx = slice(700,1300,1)

   # Loop through all files in the directory
   for f in files:
      if f[-5:]!=".fits":
         continue
      if f[0:4]!="FORS":
         continue
      try:
         obj = pyfits.getval(f,'HIERARCH ESO OBS TARG NAME').strip().upper()
      except:
         continue
      if obj.find('BIAS')>=0:
         _bias.append(f)
         continue
      if obj.find('STANDARD')>=0:
         continue
      try:
         filt = pyfits.getval(f,'FILTER1').strip()
      except:
         continue
      if filter!=filt:
         continue
      if obj.find('FLAT')>=0:
         _flat.append(f)
      else:
         _sky.append(f)
      if obj.find(object)>=0:
         _sci.append(f)


   if len(_sky)==0:
      print ""
      print "No science frames.  Stopping here."
      return

#   if filter=="I":
#      if infringe is None:
#         fringe = make_fringe(_sky,bias,flat)
#      else:
#         print "Using %s as fringe frame" % infringe
#         fringe = pyfits.open(infringe)[0].data.copy()
#      totbias = bias + fringe
#   else:
#      totbias = bias

   # Temporarily set holders for bias and flat files, and default for
   #  sky.
   arrx = x.stop - x.start
   arry = y.stop - y.start
   totbias = scipy.zeros((arry,arrx))
   flat = 1. + totbias
   sky = flat

   if insky is None:
      sky = scipy.empty((flat.shape))
      tmp = []
      print "Creating sky superflat frame for filter = %s.  " % filter
      print " Input frames:"
      for f in _sky:
         print "   %s" % f
         stmp = pyfits.open(f)[0]
         #t = subtract_overscan(stmp.data,x,y)
         #t = gain_correct(f,t)
         t = stmp.data
         t /= stats.stats.median(t[staty,statx],axis=None)
         tmp.append(t.copy())
      #
      # Now assume that the observer has taken enough exposures
      #  so that the median is actually a decent estimator of the
      #  median that would be obtained after bad pixels were rejected
      #
      tmp2 = scipy.empty((len(tmp),tmp[0].shape[0],tmp[0].shape[1]))
      print " Co-adding frames..."
      for i in range(tmp2.shape[0]):
         tmp2[i] = tmp[i].copy()
      #tmp3 = stats.stats.median(tmp2,axis=0)
      tmp3 = numpy.minimum(tmp2[0,:,:],tmp2[1,:,:],tmp2[2,:,:])
      #
      # Now do a large-scale smoothing of the co-added frame
      #
      #print " Smoothing the co-added frame. May take a while..."
      #sky = ndimage.median_filter(tmp3,51)
      sky = tmp3
      sky /= sky.mean()
      band = pyfits.open(_sky[0])[0].header['FILTER1'].strip()
      print " Writing skyflat_%s.fits" % band
      pyfits.PrimaryHDU(sky).writeto('skyflat_%s.fits' % band)
   else:
     print "Using %s as skyflat frame" % insky
     sky = pyfits.open(insky)[0].data.copy()

   totalflat = flat*sky

   print ""
   print "Correcting science frames:"
   for f in _sci:
      print "  %s" % f
      sci = pyfits.open(f)[0]
      #sci.data = subtract_overscan(sci.data,x,y)
      #sci.data = gain_correct(f,sci.data) / totalflat
      sci.data /= totalflat
#      shape = sci.data.shape
#      ra = wcs.ra2deg(sci.header['RA'].strip())
#      dec = wcs.dec2deg(sci.header['DEC'].strip())
#      sci.header.update('CRPIX1',shape[1]/2.)
#      sci.header.update('CRVAL1',ra)
#      sci.header.update('CRPIX2',shape[0]/2.)
#      sci.header.update('CRVAL2',dec)
#      sci.header.update('CD1_1',-0.371/3600.)
#      sci.header.update('CD1_2',0.)
#      sci.header.update('CD2_1',0.)
#      sci.header.update('CD2_2',0.370/3600.)
#      sci.header.update('CTYPE1','RA---TAN')
#      sci.header.update('CTYPE2','DEC--TAN')
#      sci.header.update('EQUINOX',2000.0)
#      sci.header.update('RADESYS','FK5')
#      del sci.header['CDELT1']
#      del sci.header['CDELT2']
#
      sci.verify('fix')
      newname = f.replace('FORS','%s_%s_' % (object,filter))
      newname = newname.strip()
      sci.writeto(newname)

#pfcam_redux("B1359","B","../Calib/bias.fits","../Calib/flat_B.fits","skyflat_B.fits")
#pfcam_redux("B1555","B","../Calib/bias.fits","../Calib/flat_B.fits","skyflat_B.fits")
#pfcam_redux("B1359","G","../Calib/bias.fits","../Calib/flat_G.fits")
#pfcam_redux("B1555","G","../Calib/bias.fits","../Calib/flat_G.fits","skyflat_G.fits")
#pfcam_redux("B1359","I","../Calib/bias.fits","../Calib/flat_I.fits","skyflat_I.fits","fringe_I.fits")
#pfcam_redux("B1555","I","../Calib/bias.fits","../Calib/flat_I.fits","skyflat_I.fits","fringe_I.fits")
forsim_redux("0230","I_BESS",None,None,None)
