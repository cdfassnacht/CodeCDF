try:
    from astropy.io import fits as pf
except:
    try:
        import pyfits as pf
    except:
        print ''
        print 'ERROR: Could not import pyfits or astropy.io.fits'
        print ''
        exit()
from matplotlib import pyplot as plt
import numpy as n
import imfuncs as imf
import sys

subimsize = 21

""" Check command line syntax """
if len(sys.argv)<2:
    print ''
    print 'Usage: python test_mouse_and_keypress.py [fitsfile]'
    print ''
    exit()
infits = sys.argv[1]
print ''
print 'Running code on input fits file %s' % infits
print ''

#def keypress_info():
#    """
#    Prints useful information about what the key presses do
#    """
#    print ''
#    print 'Actions available by pressing a key in the Figure 1 window'
#    print '----------------------------------------------------------'
#    print 'Key  Action'
#    print '---  ---------------------------------'
#    print ' m   Mark the position of an object'
#    print ' q   Quit and save the position of the last marked object'
#    print ''
#
#def onclick(event):
#    """
#    Actions taken if a mouse button is clicked
#    """
#
#    xd,yd = event.xdata,event.ydata
#    print xd,yd
#    return
#
#def keypress(event):
#    """
#    Actions taken if a key on the keyboard is pressed
#    """
#    if event.key == 'm':
#        global xmark, ymark
#        print ''
#        print 'Marking position %8.2f %8.2f' % (event.xdata,event.ydata)
#        print ''
#        xmark = event.xdata
#        ymark = event.ydata
#        fig2 = plt.figure(2,figsize=(10,3))
#        fig2.add_subplot(131)
#        im1.display(subimcent=(xmark,ymark),subimsize=(subimsize,subimsize))
#        fig2.add_subplot(132)
#        xsum = im1.subim.sum(axis=0)
#        plt.plot(xsum)
#        plt.xlabel('Relative x Coord')
#        fig2.add_subplot(133)
#        ysum = im1.subim.sum(axis=1)
#        plt.plot(ysum)
#        plt.xlabel('Relative y Coord')
#        fig2.show()
#    if event.key == 'q':
#        print ''
#        print 'Closing down'
#        print ''
#        fig1.canvas.mpl_disconnect(cid1)
#        fig1.canvas.mpl_disconnect(cid2)
#        plt.close(1)
#        plt.close(2)
#
#    keypress_info()
#    return

im1 = imf.Image(infits)
im1.zoomsize = subimsize
#fig2 = plt.figure(2)
#fig1 = plt.figure(1)
im1.display(sighigh=30.)
im1.keypress_info()
#im1.cid_mouse = fig1.canvas.mpl_connect('button_press_event',im1.onclick)
#im1.cid_keypress = fig1.canvas.mpl_connect('key_press_event',im1.keypress)
im1.start_interactive()

plt.show()

print 'Got final marked position of %8.2f %8.2f' % (im1.xmark,im1.ymark)
print ''
