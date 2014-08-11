import numpy as n
import scipy,pyfits,math
from scipy import arccos,arcsin,arctan,arctan2,sin,cos,tan,pi,linalg

raddeg = 180./pi

# Convert pixels to sky coordinates, including projection
def pix2sky(header,x,y):
	hdr_info = parse_header(header)
	x0 = x-hdr_info[1][0]+1.	# Plus 1 python->image
	y0 = y-hdr_info[1][1]+1.
	x0 = x0.astype(scipy.float64)
	y0 = y0.astype(scipy.float64)
	x = hdr_info[2][0,0]*x0 + hdr_info[2][0,1]*y0
	y = hdr_info[2][1,0]*x0 + hdr_info[2][1,1]*y0
	if hdr_info[3]=="DEC":
		a = x.copy()
		x = y.copy()
		y = a.copy()
		ra0 = hdr_info[0][1]
		dec0 = hdr_info[0][0]/raddeg
	else:
		ra0 = hdr_info[0][0]
		dec0 = hdr_info[0][1]/raddeg
	if hdr_info[5]=="TAN":
		r_theta = scipy.sqrt(x*x+y*y)/raddeg
		theta = arctan(1./r_theta)
		phi = arctan2(x,-1.*y)
	elif hdr_info[5]=="SIN":
		r_theta = scipy.sqrt(x*x+y*y)/raddeg
		theta = arccos(r_theta)
		phi = artan2(x,-1.*y)
	ra = ra0 + raddeg*arctan2(-1.*cos(theta)*sin(phi-pi),
				   sin(theta)*cos(dec0)-cos(theta)*sin(dec0)*cos(phi-pi))
	dec = raddeg*arcsin(sin(theta)*sin(dec0)+cos(theta)*cos(dec0)*cos(phi-pi))

	return ra,dec

# Convert sky coordinates to grid (ie ccd) coordinates, including projection
def sky2pix(header,ra,dec):
	hdr_info = parse_header(header)
	if scipy.isscalar(ra):
		ra /= raddeg
		dec /= raddeg
	else:
		ra = ra.astype(scipy.float64)/raddeg
		dec = dec.astype(scipy.float64)/raddeg
	if hdr_info[3]=="DEC":
		ra0 = hdr_info[0][1]/raddeg
		dec0 = hdr_info[0][0]/raddeg
	else:
		ra0 = hdr_info[0][0]/raddeg
		dec0 = hdr_info[0][1]/raddeg

	phi = pi + arctan2(-1*cos(dec)*sin(ra-ra0),sin(dec)*cos(dec0)-cos(dec)*sin(dec0)*cos(ra-ra0))
	argtheta = sin(dec)*sin(dec0)+cos(dec)*cos(dec0)*cos(ra-ra0)
	if scipy.isscalar(argtheta):
		theta = arcsin(argtheta)
	else:
		argtheta[argtheta>1.] = 1.
		theta = arcsin(argtheta)

	if hdr_info[5]=="TAN":
		r_theta = raddeg/tan(theta)
		x = r_theta*sin(phi)
		y = -1.*r_theta*cos(phi)
	elif hdr_info[5]=="SIN":
		r_theta = raddeg*cos(theta)
		x = r_theta*sin(phi)
		y = -1.*r_theta*cos(phi)
	if hdr_info[3]=="DEC":
		a = x.copy()
		x = y.copy()
		y = a.copy()
	inv = linalg.inv(hdr_info[2])
	x0 = inv[0,0]*x + inv[0,1]*y
	y0 = inv[1,0]*x + inv[1,1]*y

	x = x0+hdr_info[1][0]-1
	y = y0+hdr_info[1][1]-1

	return x,y

# Create a wcs header given a central RA/DEC and image size and scale
def make_header(ra,dec,xsize,ysize,xscale,yscale=-1,docdmatx=True):
	"""
	Inputs:
	  ra
	  dec
	  xsize
	  ysize
	  xscale
	  yscale
	  docdmatx - If set to True (the default), then put the image scale
	             in terms of a CD matrix.  If False, then use the
		     CDELT and PC matrix formalism
	"""
	if yscale==-1:
		yscale = xscale
	hdr = pyfits.PrimaryHDU(scipy.empty((ysize,xsize))).header.copy()
	if is_degree(ra)==False:
		ra = ra2deg(ra)
	if is_degree(dec)==False:
		dec = dec2deg(dec)
	hdr.update("RA",ra)
	hdr.update("DEC",dec)
	hdr.update("CTYPE1","RA---TAN")
	hdr.update("CTYPE2","DEC--TAN")
	hdr.update("CRVAL1",ra)
	hdr.update("CRPIX1",xsize/2.)
	hdr.update("CRVAL2",dec)
	hdr.update("CRPIX2",ysize/2.)
	if docdmatx:
		hdr.update('cd1_1',-xscale/3600.)
		hdr.update('cd1_2',0.)
		hdr.update('cd2_1',0.)
		hdr.update('cd2_2',yscale/3600.)
	else:
		hdr.update("CDELT1",xscale/3600.)
		hdr.update("CDELT2",yscale/3600.)
		hdr.update("PC1_1",-1.)
		hdr.update("PC1_2",0.)
		hdr.update("PC2_1",0.)
		hdr.update("PC2_2",1.)

	return hdr

# Apply a rotation to a wcs header
def rotate_header(oheader,angle):
	DEGRAD = pi/180.
	header = oheader.copy()
	header['PC1_1'] = cos(DEGRAD*angle)*oheader['PC1_1']+sin(DEGRAD*angle)*oheader['PC2_1']
	header['PC1_2'] = cos(DEGRAD*angle)*oheader['PC1_2']+sin(DEGRAD*angle)*oheader['PC2_2']
	header['PC2_1'] = -1.*sin(DEGRAD*angle)*oheader['PC1_1']+cos(DEGRAD*angle)*oheader['PC2_1']
	header['PC2_2'] = -1.*sin(DEGRAD*angle)*oheader['PC1_2']+cos(DEGRAD*angle)*oheader['PC2_2']
	return header

#-----------------------------------------------------------------------
# Parse the relevent wcs headers
def parse_header(header):
	crval = scipy.array([header['crval1'],header['crval2']])
	crpix = scipy.array([header['crpix1'],header['crpix2']])
	cd = scipy.zeros((2,2))
	try:
		cd[0,0] = header['cd1_1']
		cd[1,1] = header['cd2_2']
		try:
			cd[0,1] = header['cd1_2']
		except:
			cd[0,1] = 0.0
		try:
			cd[1,0] = header['cd2_1']
		except:
			cd[1,0] = 0.0
	except:
		try:
			cd[0,0] = header['pc1_1']
			cd[0,1] = header['pc1_2']
			cd[1,0] = header['pc2_1']
			cd[1,1] = header['pc2_2']
			cd[0] *= header['cdelt1']
			cd[1] *= header['cdelt2']
		except:
			try:
				rota2 = header['crota2']
			except:
				#rota2 = header['crota1']
				rota2 = 0.0
			cdelt1 = header['cdelt1']
			cdelt2 = header['cdelt2']
			cd[0,0] = cdelt1*cos(rota2)
			cd[0,1] = -1.*cdelt2*sin(rota2)
			cd[1,0] = cdelt1*sin(rota2)
			cd[1,1] = cdelt2*cos(rota2)
	ctype1,projection = parse_ctype(header['ctype1'])
	ctype2,projection = parse_ctype(header['ctype2'])

	return [crval,crpix,cd,ctype1,ctype2,projection]

#-----------------------------------------------------------------------

# Parse the ctype header card
def parse_ctype(ctype):
	tmp = ctype.split("-")
	while 1:
		try:
			tmp.remove("")
		except:
			break
	return tmp[0],tmp[1]

# Check if the coordinate is in degrees
def is_degree(comp):
	if type(comp)==float or type(comp)==n.float64 or type(comp)==n.float32:
		return True
	return False

# Convert string ra to degrees
def ra2deg(ra):
	comp = ra.split(" ")
	if comp[0]==ra:
		comp = comp[0].split(":")
	deg = float(comp[0])+float(comp[1])/60.+float(comp[2])/3600.
	return deg*15.

# Convert string declination to degrees
def dec2deg(dec):
	comp = dec.split(" ")
	if comp[0]==dec:
		comp = comp[0].split(":")
	if comp[0][0]=="-":
		comp[0] = comp[0][1:]
		sign = -1.
	else:
		sign = 1.
	return sign*(float(comp[0])+float(comp[1])/60.+float(comp[2])/3600.)

# Convert decimal ra to HMS format
def deg2ra(ra,sep=" "):
	ra /= 15.
	h = math.floor(ra)
	res = (ra-h)*60.
	m = math.floor(res)
	s = (res-m)*60.
	if sep=="hms":
		sep1 = "h"
		sep2 = "m"
		sep3 = "s"
	else:
		sep1 = sep
		sep2 = sep
		sep3 = ""
	return "%02d%s%02d%s%06.3f%s" % (h,sep1,m,sep2,s,sep3)

# Convert decimal declination to DaMaS format
def deg2dec(dec,sep=" "):
	if dec<0:
		sign = -1.
		dec = abs(dec)
	else:
		sign = 1.
	d = math.floor(dec)
	res = (dec-d)*60.
	m = math.floor(res)
	s = (res-m)*60.
	if sep=="dms":
		sep1 = "d"
		sep2 = "m"
		sep3 = "s"
	else:
		sep1 = sep
		sep2 = sep
		sep3 = ""
	if sign==-1:
		return "-%02d%s%02d%s%06.3f%s" % (d,sep1,m,sep2,s,sep3)
	return "%02d%s%02d%s%06.3f%s" % (d,sep1,m,sep2,s,sep3)
