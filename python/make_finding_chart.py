import wcs,scipy,pyfits,pylab,sys
from scipy import ndimage

def postage_stamp(input,output,locations,xsize,ysize,scale,angle):
	pylab.close()

#	f = open(locations).readlines()[0].split()
#	ra = wcs.ra2deg(f[3]+":"+f[4]+":"+f[5])
#	dec = wcs.dec2deg(f[6]+":"+f[7]+":"+f[8])
#

#	if wcs.is_degree(ra)==False:
#		ra = wcs.ra2deg(ra)
#	if wcs.is_degree(dec)==False:
#		dec = wcs.dec2deg(dec)
	ra = 317.72512
	dec = 21.516883
	outheader = wcs.make_header(ra,dec,xsize,ysize,scale)
	outheader = wcs.rotate_header(outheader,angle)
	coords = scipy.indices((ysize,xsize)).astype(scipy.float32)
	skycoords = wcs.pix2sky(outheader,coords[1],coords[0])

	file = pyfits.open(input)
	inheader = file[0].header.copy()
	data = file[0].data.copy()
	ccdcoords = wcs.sky2pix(inheader,skycoords[0],skycoords[1])
	coords[1] = ccdcoords[0]
	coords[0] = ccdcoords[1]

	image = ndimage.map_coordinates(data,coords,output=scipy.float64)
	bounds = scipy.sort(image.flatten())
	vmin = bounds[bounds.size*0.65]
	vmax = bounds[bounds.size*0.995]
	pylab.imshow(image[::-1],cmap=pylab.cm.gist_yarg,vmin=vmin,vmax=vmax)

	pylab.axis('off')
	title = r"B2108 Finding Chart"
	
	pylab.text(420,320,r"Star TO Target offsets:")
	pylab.text(440,370,r"11.91$^{\prime\prime}$ S, 7.20$^{\prime\prime}$ W")
	pylab.text(440,420,r"Slit PA 60 degrees E from N")
	pylab.title(title)

	length = 10./scale
	dx = length
	dy = length
	pylab.arrow(300,300,-1.*dx,0)
	pylab.arrow(300,300,0,-1*dx)

	pylab.rc('text',usetex=True)
	pylab.text(310,290-dy,'N')
	pylab.text(310.,290-dy/2.,r'10$^{\prime\prime}$')

	ax = pylab.gca()

	ax.figure.set_size_inches((7.5,7.5))
	import matplotlib as mpl
	a = 0
	for i in locations:
		ra = i[0]
		dec = i[1]
		x,y = wcs.sky2pix(outheader,ra,dec)
		y -= ysize/2.
		y *= -1
		y += ysize/2.
		if a==0:
			a = 1
			pylab.text(x+30,y-30,"Target")
		else:
			pylab.text(x+30,y-30,"Offset Star")
		ax.patches.append(pylab.Circle((x,y),25,transform=ax.transData,fill=False,ec='r',lw=1.5))
		
	pylab.savefig(output+".eps")
	import os
	os.system('/usr/bin/convert %s.eps %s.png' % (output,output))

inname = sys.argv[1]
outname = sys.argv[2]

loc = [[317.72512,21.516883],[317.72868,21.518882]]
postage_stamp(inname,outname,loc,80/0.05,80/0.05,0.05,0.)

