"""
Program that takes an input fits file with WCS information and a list of
(RA,Dec) coordinates and then displays a small cutout image around each
(RA,Dec) coordinate.

Inputs:
   fitsfile  -  Input fits file with WCS information
   catfile   -  Input ascii file containing RA and Dec as decimal degrees
                in separate columns
   imsize    -  cutout image size in arcsec
   fmax      -  [OPTIONAL] upper end of display range, in sigma 
                above the mean.  Default=25
   racol     -  [OPTIONAL] Column in catfile containing the RA in decimal
                degrees.  Default value (1) means that RA is in the second
                column (since python is zero-indexed)
   deccol    -  [OPTIONAL] Column in catfile containing the Dec in decimal
                degrees.  Default value (2) means that Dec is in the third
                column (since python is zero-indexed)
   idcol     -  [OPTIONAL] Column in catfile containing the ID. '
                Default value (0) means that the ID '
                is in the first column (since python is zero-indexed)'

Usage: 
 python display_cutouts.py [fitsfile] [catfile] [imsize] ([fmax] 
    [racol] [deccol] [idcol])

"""

import sys
from matplotlib import pyplot as plt
try:
    from SpecIm import imfuncs as imf
except ImportError:
    import imfuncs as imf
import catfuncs as cf

""" Parse the command line """
if len(sys.argv) < 4:
    print ''
    print 'display_cutouts.py'
    print '------------------------'
    print 'Program that takes an input fits file with WCS information and a '
    print 'list of (RA,Dec) coordinates and then displays a small cutout image'
    print 'around each (RA,Dec) coordinate.'
    print ''
    print 'Inputs:'
    print '   fitsfile -  Input fits file with WCS information'
    print '   catfile  -  Input ascii file containing RA and Dec as decimal'
    print '               degrees in separate columns'
    print '   imsize   -  cutout image size in arcsec'
    print '   fmax     -  [OPTIONAL] upper end of display range, in sigma '
    print '               above the mean.  Default=25'
    print '   racol    -  [OPTIONAL] Column in catfile containing the RA in '
    print '               decimal degrees.  Default value (1) means that RA '
    print '               is in the second column (since python is zero-indexed)'
    print '   deccol   -  [OPTIONAL] Column in catfile containing the Dec in '
    print '               decimal degrees.  Default value (2) means that Dec '
    print '               is in the third column (since python is zero-indexed)'
    print '   idcol    -  [OPTIONAL] Column in catfile containing the ID. '
    print '               Default value (0) means that the ID '
    print '               is in the first column (since python is zero-indexed)'
    print ''
    print 'Examples of usage'
    print '------------------------'
    print 'python display_cutouts.py [fitsfile] [catfile] [imsize] [fmax]'
    print ''
    exit()

# Simple mouse click function to move on by closing the window
def onclick(event):
    global ix,iy
    ix, iy = event.xdata, event.ydata

    # assign global variable to access outside of function
    global coords
    coords.append((ix, iy))

    # Disconnect after 1 click
    if len(coords) == 1:
        fig.canvas.mpl_disconnect(cid)
        plt.close(1)
    return

""" Assign the variables """
fitsfile = sys.argv[1]
catfile = sys.argv[2]
imsize = float(sys.argv[3])
if len(sys.argv)>4:
    fmax = float(sys.argv[4])
else:
    fmax = 25.
if len(sys.argv)>5:
    racol = int(sys.argv[5])
else:
    racol = 1
if len(sys.argv)>6:
    deccol = int(sys.argv[6])
else:
    deccol = 2
if len(sys.argv)>7:
    namecol = int(sys.argv[7])
else:
    namecol = 0

print ''
print 'Paramter values'
print '----------------'
print 'imsize:  %4.1f arcsec' % imsize
print 'racol:   %d' % racol
print 'deccol:  %d' % deccol
print 'namecol: %d' % namecol
print 'fmax:    %5.1f' % fmax


""" Open the fits file """
try:
    infits = imf.Image(fitsfile)
except:
    print ''
    print 'ERROR: could not open input fits file: %s' % fitsfile
    print ''
    exit()

""" Open the catalog file """
try:
    incat = cf.Secat(catfile,catformat='ascii',racol=racol,deccol=deccol,
                     namecol=namecol)
except:
    print ''
    print 'ERROR: problem reading input catalog file: %s' % catfile
    print ''
    infits.close()
    exit()

""" Loop through the (RA,Dec) pairs """
incat.get_radec()
print ''
xx = []
yy = []
for i in range(incat.ra.size):
    print ''
    print 'Image Center: %11.7f %+11.7f' % (incat.ra[i],incat.dec[i])
    title = 'ID: %s' % incat.data[incat.namefield][i]
    #fig = plt.figure(1)
    #ax = fig.add_subplot(111)
    infits.display(subimdef='radec',subimcent=(incat.ra[i],incat.dec[i]),
                   subimsize=(imsize,imsize),dispunits='radec',fmax=fmax,
                   title=title,show_xyproj=True)
    fig = plt.gcf()
    coords = []

    """ Call click func """
    print ''
    print 'Click on the figure to move to the next cutout'
    print ''
    cid = fig.canvas.mpl_connect('button_press_event', onclick)

    plt.show()


""" Clean up and exit """
print ''
infits.close()
del infits,incat
    





