import pyfits,scipy
import sys,os
import wcs
from scipy import stats,ndimage

def make_fringe(frames,bias,flat):
	#
	# Define the good area of the chip.  Only read in rows and columns
	#  defined by x and y below
	#
	y = slice(150,1950,1)
	x = slice(150,2020,1)
        #
        # Define an area of the chip to use for statistics
        #
	fstaty = slice(1500,2500,1)
	fstatx = slice(1500,2500,1)

	fringe = scipy.empty(bias.shape)
	tmp = []
	band = pyfits.open(frames[0])[0].header['DWFILNAM'].strip()
	print "Creating fringe frame for filter = %s.  Input frames:" % band
	for f in frames:
		print "   %s" % f
		t = (pyfits.open(f)[0].data[y,x].astype(scipy.float32)-bias)/flat
		t -= stats.stats.median(t[fstaty,fstatx],axis=None)
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
	tmp3 = stats.stats.median(tmp2,axis=0)
	#
	# Now do a small-scale smoothing of the co-added frame
	#
	print " Smoothing the co-added frame"
	fringe = ndimage.median_filter(tmp3,11)
	pyfits.PrimaryHDU(fringe).writeto('fringe_%s.fits' % band)

	return fringe

def pfcam_redux(object,filter,inbias=None,inflat=None,insky=None,infringe=None):
	dir = os.curdir
	files = os.listdir(dir)

	_bias = []
	_flat = []
	_sky = []
	_sci = []

	#
	# Define the good area of the chip.  Only read in rows and columns
	#  defined by x and y below
	#
	x = slice(150,2020,1)
	y = slice(150,1950,1)
        #
        # Define an area of the chip to use for statistics
        #
	staty = slice(1500,2500,1)
	statx = slice(1500,2500,1)

	for f in files:
		if f[-5:]!=".fits":
			continue
		if f[0:3]!="pfc":
			continue
		try:
			obj = pyfits.getval(f,'OBJECT').strip().upper()
		except:
			continue
		if obj.find('BIAS')>=0:
			_bias.append(f)
			continue
		if obj.find('STANDARD')>=0:
			continue
		try:
			filt = pyfits.getval(f,'DWFILNAM').strip()
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

	if inbias is None:
		print ""
		print "Creating bias frame.  Input frames:"
		bias = []
		for f in _bias:
			bias.append(pyfits.open(f)[0].data[y,x].astype(scipy.float32))
			print "   %s" % f
		tmp = scipy.empty((len(bias),bias[0].shape[0],bias[0].shape[1]))
		for i in range(tmp.shape[0]):
			tmp[i] = bias[i].copy()
		bias = stats.stats.median(tmp,axis=0)
		pyfits.PrimaryHDU(bias).writeto('bias.fits')
		print " Done: Created file bias.fits"
	else:
		print "Using %s as bias frame" % inbias
		bias = pyfits.open(inbias)[0].data.copy()

	if inflat is None:
		flat = []
		print ""
		print "Creating domeflat frame for filter = %s.  Input frames:" % filter
		for f in _flat:
			flat.append(pyfits.open(f)[0].data[y,x].astype(scipy.float32)-bias)
			print "   %s" % f
		tmp = scipy.empty((len(flat),flat[0].shape[0],flat[0].shape[1]))
		for i in range(tmp.shape[0]):
			tmp[i] = flat[i].copy()
		flat = stats.stats.median(tmp,axis=0)
		flat /= flat.mean()
		band = pyfits.open(_flat[0])[0].header['DWFILNAM'].strip()
		pyfits.PrimaryHDU(flat).writeto("flat_%s.fits" % band)
		print " Done: Created file flat_%s.fits" % band
	else:
		print "Using %s as domeflat frame" % inflat
		flat = pyfits.open(inflat)[0].data.copy()

	if len(_sky)==0:
		print ""
		print "No science frames.  Stopping here."
		return

	if filter=="I":
		if infringe is None:
			fringe = make_fringe(_sky,bias,flat)
		else:
			print "Using %s as fringe frame" % infringe
			fringe = pyfits.open(infringe)[0].data.copy()
		totbias = bias + fringe
	else:
		totbias = bias

	if insky is None:
		sky = scipy.empty(flat.shape)
		tmp = []
		print "Creating sky superflat frame for filter = %s.  Input frames:" % filter
		for f in _sky:
			print "   %s" % f
			t = (pyfits.open(f)[0].data[y,x].astype(scipy.float32)-totbias)/flat
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
		tmp3 = stats.stats.median(tmp2,axis=0)
		#
		# Now do a large-scale smoothing of the co-added frame
		#
		print " Smoothing the co-added frame"
		sky = ndimage.median_filter(tmp3,51)
		sky /= sky.mean()
		band = pyfits.open(_sky[0])[0].header['DWFILNAM'].strip()
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
		sci.data = ((sci.data[y,x].astype(scipy.float32)-totbias)/totalflat).T[::-1,::-1]
		shape = sci.data.shape
		ra = wcs.ra2deg(sci.header['RA'].strip())
		dec = wcs.dec2deg(sci.header['DEC'].strip())
		sci.header.update('CRPIX1',shape[1]/2.)
		sci.header.update('CRVAL1',ra)
		sci.header.update('CRPIX2',shape[0]/2.)
		sci.header.update('CRVAL2',dec)
		sci.header.update('CD1_1',-0.371/3600.)
		sci.header.update('CD1_2',0.)
		sci.header.update('CD2_1',0.)
		sci.header.update('CD2_2',0.370/3600.)
		sci.header.update('CTYPE1','RA---TAN')
		sci.header.update('CTYPE2','DEC--TAN')
		sci.header.update('EQUINOX',2000.0)
		sci.header.update('RADESYS','FK5')
		del sci.header['CDELT1']
		del sci.header['CDELT2']

		sci.verify('fix')
		newname = f.replace('pfc','%s_%s_' % (object,filter))
		newname = newname.strip()
		sci.writeto(newname)

#pfcam_redux("B1359","B","../Calib/bias.fits","../Calib/flat_B.fits","skyflat_B.fits")
#pfcam_redux("B1555","B","../Calib/bias.fits","../Calib/flat_B.fits","skyflat_B.fits")
#pfcam_redux("B1359","G","../Calib/bias.fits","../Calib/flat_G.fits")
#pfcam_redux("B1555","G","../Calib/bias.fits","../Calib/flat_G.fits","skyflat_G.fits")
#pfcam_redux("B1359","I","../Calib/bias.fits","../Calib/flat_I.fits","skyflat_I.fits","fringe_I.fits")
#pfcam_redux("B1555","I","../Calib/bias.fits","../Calib/flat_I.fits","skyflat_I.fits","fringe_I.fits")
