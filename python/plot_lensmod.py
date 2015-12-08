"""
Functions that are used to plot the output from a lens modeling run.
Right now these functions work with files that are produced by
glafic (Masamune Oguri) and gravlens (Chuck Keeton)
"""

import numpy as n
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib import patches

#-----------------------------------------------------------------------

def plot_critcaust(infile, plottype, ax=None, icolor='b', scolor='r', ils='-',
                   sls='-'):
    """

    Plots either the critical curves or the caustics contained in the input file.
    This file is expected to be in the format produced by both glafic and
    gravlens, namely as a series of line segments defined by

       xi1 yi1 xs1 ys1 xi2 yi2 xs2 ys2

    where xi1, yi1, etc., denote positions in the image plane and are used
    for the critical curve(s), and
    where xs1, ys1, etc., denote positions in the source plane and are used
    for the caustics.

    Note: The simplest coding to draw the curves, once the arrays have been
    loaded, is just: plot([xi1,xi2],[yi1,yi2]).
    However, even if the xi1, etc. are vectors, this is MUCH slower than
     it has to be.
    The code in this function is much faster, since only a single 
    LineCollection object needs to be plotted rather than hundreds of 
    Line2D objects.

    Inputs:
     infile   - input file containing the line segment info
     plottype - which curves to plot, either 'crit' or 'caust'
     ax       - the matplotlib axis to which to add the plotted curve.
                The default (None) should be fine for most applications.
     icolor   - color for the critical curves.  Default is 'b'
     scolor   - color for the caustic curves.   Default is 'r'

    """

    """ Load the data """
    print ""
    try:
        xi1,yi1,xs1,ys1,xi2,yi2,xs2,ys2 = n.loadtxt(infile,unpack=True)
    except:
        print ''
        print "ERROR: plot_crit. Unable to read %s." % infile
        print "Will not plot critical curves and caustics"
        print ""
        return
    print "Read data from input file %s" % infile
   
    """ Set up the line collection object.  Start by creating a 3d array """
    segs_i = n.zeros((xi1.size,2,2))
    segs_s = n.zeros((xs1.size,2,2))
    for i in range(segs_i.shape[0]):
        segs_i[i,0,0] = xi1[i]
        segs_i[i,1,0] = xi2[i]
        segs_i[i,0,1] = yi1[i]
        segs_i[i,1,1] = yi2[i]
        segs_s[i,0,0] = xs1[i]
        segs_s[i,1,0] = xs2[i]
        segs_s[i,0,1] = ys1[i]
        segs_s[i,1,1] = ys2[i]
    line_segs_i = LineCollection(segs_i,colors=icolor,linewidths=2.,
                                 linestyles=ils)
    line_segs_s = LineCollection(segs_s,colors=scolor,linewidths=2.,
                                 linestyles=sls)

    """ Actually do the plotting.  Only plot the requested curve """
    if plottype[:6]=='caust':
        print "Plotting caustics"
        if ax is None:
            ax = plt.axes()
        ax.set_xlim(xs1.min(), xs2.max())
        ax.set_ylim(ys1.min(), ys2.max())
        ax.add_collection(line_segs_s)
        #ax_s.set_xlim(xs1.min(), xs2.max())
        #ax_s.set_ylim(ys1.min(), ys2.max())
        #ax_s.add_collection(line_segs_s)
    else:
        print "Plotting critical curves"
        xmini = min((xi1.min(),xi2.min()))
        xmaxi = max((xi1.max(),xi2.max()))
        xd = xmaxi - xmini
        xmini -= 0.05*xd
        xmaxi += 0.05*xd
        ymini = min((yi1.min(),yi2.min()))
        ymaxi = max((yi1.max(),yi2.max()))
        yd = ymaxi - ymini
        ymini -= 0.05*yd
        ymaxi += 0.05*yd
        if ax is None:
            ax = plt.axes()
        ax.set_xlim(xmini, xmaxi)
        ax.set_ylim(ymini, ymaxi)
        ax.add_collection(line_segs_i)
        #ax_i = plt.axes()
        #ax_i.set_xlim(xmini, xmaxi)
        #ax_i.set_ylim(ymini, ymaxi)
        #ax_i.add_collection(line_segs_i)

