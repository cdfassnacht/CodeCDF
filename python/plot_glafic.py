"""
Functions used to plot the output from glafic as figures
"""

import os
import numpy as n
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib import patches

def read_glafic_point(infile):
   """
   Reads in data from a glafic *point.dat output file.
   Returns two arrays: image information and source information

   Inputs:
      infile  -  input file (should be named *point.dat)
   """

   try:
      data = n.loadtxt(infile)
   except:
      print ""
      print "ERROR: plot_points.  Could not read input file %s" % infile
      print "Will not plot image and source locations"
      print ""
      return None,None

   # Separate the source-plane and lens-plane info
   srcinfo = data[0,:]
   imginfo = data[1:,:]

   return imginfo, srcinfo

#-----------------------------------------------------------------------

def read_glafic_lensmodel(infile):
   """
   Reads input from a glafic optimized lens output file.
   NOTE: For now assumes one SIE and one shear
   Returns a dtype (sort of like a C structure) with model info.
   """

   # Open the input file
   try:
      f = open(infile)
   except:
      print ""
      print "ERROR.  Could not open input file %s" % infile
      return None

   # Get the number of SIEs in the model file
   nsie = 0
   ilines = f.readlines()
   for i in ilines:
      ic = i.split()
      if len(ic)>1:
         if ic[0]=='lens' and ic[1]=='sie':
            nsie += 1
   if nsie == 0:
      print ""
      print "ERROR: read_glafic_lensmod.  No SIE models found in input file %s" \
          % infile
      return None

   # Create a dtype container
   modpars = n.zeros(nsie,dtype={'names':['siex','siey','vdisp','e','theta_e',
                                       'rcore','zs','gamma','theta_gamma',
                                       'kappa'],
                              'formats':[n.float64,n.float64,n.float64,
                                         n.float64,n.float64,n.float64,
                                         n.float64,n.float64,n.float64,
                                         n.float64]})

   # Read the file
   count = 0
   for i in ilines:
      ic = i.split()
      if len(ic)>1:
         if ic[0]=='lens' and ic[1]=='sie':
            modpars[count]['vdisp']   = ic[2]
            modpars[count]['siex']    = ic[3]
            modpars[count]['siey']    = ic[4]
            modpars[count]['e']       = ic[5]
            modpars[count]['theta_e'] = ic[6]
            modpars[count]['rcore']   = ic[7]
            count += 1
         if ic[0]=='lens' and ic[1]=='pert':
            modpars[0]['zs']      = ic[2]
            modpars[0]['gamma']   = ic[5]
            modpars[0]['theta_gamma'] = ic[6]
            modpars[0]['kappa']   = ic[8]
   return modpars

#-----------------------------------------------------------------------

def plot_mesh(prefix, suff='mesh.dat'):
   """
   Plots the mesh pattern contained in prefix_mesh.dat.  

   Each line of the input file contains the start and end points of a line
    segment that, when all are drawn, will show the mesh used for the modeling.

   The format of the input file is:
       xi1 yi1 xs1 ys1 xi2 yi2 xs2 ys2
    where the "i" variables are used to draw the image-plane mesh (i.e., the
    line segments go from (xi1,yi1) to (xi2,yi2)) and the "s" variables are
    used for the source-plane mesh.

   Note: The simplest coding to draw the mesh, once the arrays have been
    loaded, is just:

      plot([xi1,xi2],[yi1,yi2])

    This format works even if the xi1, etc., are arrays rather than scalars.
    However, for the many thousands of line segments that can be produced
    for a mesh, this format can be quite slow since a new Line2D object has
    to get created for each array element (i.e., each line in the input file).
    The code in this function is much faster, since only a single 
    LineCollection object needs to be plotted.

   Inputs:
    prefix - the prefix part of the input file name (i.e. the part that
             preceeds "_mesh.dat"

   """

   # Load the data
   infile = "%s_%s" % (prefix,suff)
   print ""
   try:
      xi1,yi1,xs1,ys1,xi2,yi2,xs2,ys2 = n.loadtxt(infile,unpack=True)
   except:
      print "ERROR: plot_mesh. Unable to read %s." % infile
      print "Will not plot mesh"
      print ""
      return
   print "Read data from input file %s" % infile
   
   # Set up the line collection object.  Start by creating a 3d array
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
   line_segs_i = LineCollection(segs_i,colors='k')
   line_segs_s = LineCollection(segs_s,colors='k')

   # Actually do the plotting
   print "Plotting mesh (can take a while)"
   plt.figure(1)
   ax_i = plt.axes()
   ax_i.set_xlim(xi1.min(), xi2.max())
   ax_i.set_ylim(yi1.min(), yi2.max())
   ax_i.add_collection(line_segs_i)
   plt.figure(2)
   ax_s = plt.axes()
   ax_s.set_xlim(xs1.min(), xs2.max())
   ax_s.set_ylim(ys1.min(), ys2.max())
   ax_s.add_collection(line_segs_s)
   plt.show()

#-----------------------------------------------------------------------

def plot_crit(prefix, icolor='b', scolor='r', suff='crit.dat'):
   """
   Plots the critical curves and caustics contained in prefix_crit.dat.  

   Each line of the input file contains the start and end points of a line
    segment that, when all are drawn, will show the mesh used for the modeling.

   The format of the input file is:
       xi1 yi1 xs1 ys1 xi2 yi2 xs2 ys2
    where the "i" variables are used to draw the image-plane mesh (i.e., the
    line segments go from (xi1,yi1) to (xi2,yi2)) and the "s" variables are
    used for the source-plane mesh.

   Note: The simplest coding to draw the curves, once the arrays have been
    loaded, is just:

      plot([xi1,xi2],[yi1,yi2])

    This format works even if the xi1, etc., are arrays rather than scalars.
    However, even though it is more complex, the code in this function is much
    faster, since only a single LineCollection object needs to be plotted
    rather than hundreds of Line2D objects.

   Inputs:
    prefix - the prefix part of the input file name (i.e. the part that
             preceeds "_mesh.dat"

   """

   # Load the data
   infile = "%s_%s" % (prefix,suff)
   print ""
   try:
      xi1,yi1,xs1,ys1,xi2,yi2,xs2,ys2 = n.loadtxt(infile,unpack=True)
   except:
      print "ERROR: plot_crit. Unable to read %s." % infile
      print "Will not plot critical curves and caustics"
      print ""
      return
   print "Read data from input file %s" % infile
   
   # Set up the line collection object.  Start by creating a 3d array
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
   line_segs_i = LineCollection(segs_i,colors=icolor,linewidths=2.)
   line_segs_s = LineCollection(segs_s,colors=scolor,linewidths=2.)

   # Actually do the plotting
   print "Plotting critical curves"
   plt.figure(1)
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
   ax_i = plt.axes()
   ax_i.set_xlim(xmini, xmaxi)
   ax_i.set_ylim(ymini, ymaxi)
   ax_i.add_collection(line_segs_i)
   print "Plotting caustics"
   plt.figure(2)
   ax_s = plt.axes()
   ax_s.set_xlim(xs1.min(), xs2.max())
   ax_s.set_ylim(ys1.min(), ys2.max())
   ax_s.add_collection(line_segs_s)
   plt.show()

#-----------------------------------------------------------------------

def plot_lensmod(prefix):

   # Load the data
   infile = "%s_optresult.dat" % prefix
   modpars = read_glafic_lensmodel(infile)
   if modpars is None:
      print "Exiting plot_lensmod"
      return

   # Give model info
   print ""
   print "Model Parameters"
   print "----------------"
   for i in range(modpars.size):
      print "x_SIE_%d:      %8.3f"    % ((i+1),modpars[i]['siex'])
      print "y_SIE_%d:      %8.3f"    % ((i+1),modpars[i]['siey'])
      print "vdisp_SIE_%d:  %7.2f"    % ((i+1),modpars[i]['vdisp'])
      print "e_%d:             %5.3f" % ((i+1),modpars[i]['e'])
      print "theta_e_%d:    %7.2f"    % ((i+1),modpars[i]['theta_e'])
   print "gamma:          %.3g"   % modpars[0]['gamma']
   print "theta_gamma:  %7.2f"    % modpars[0]['theta_gamma']
   print "kappa:          %6.3f"  % modpars[0]['kappa']

   # Plot the lens position
   plt.figure(1)
   ax = plt.axes()
   for i in range(modpars.size):
      q = 1. - modpars[i]['e']
      a = 0.5
      b = q * a
      lens = patches.Ellipse((modpars[i]['siex'],modpars[i]['siey']),a,b,
                             modpars[i]['theta_e']+90.,color='r',alpha=0.3)
      ax.add_patch(lens)

#-----------------------------------------------------------------------

def plot_points(prefix, title=None, show_src_in_implane=True, scale=0.03,
                xmin=None, xmax=None, ymin=None, ymax=None, dtboxloc='tr'):

   # Load the data
   infile = "%s_point.dat" % prefix
   imginfo,srcinfo = read_glafic_point(infile)
   if imginfo is None:
      print "Exiting plot_points"
      return

   # Plot the source location
   plt.figure(2)
   xs = srcinfo[2]
   ys = srcinfo[3]
   plt.plot(xs,ys,'go',ms=10)
   #plt.scatter(xs,ys,s=100,c='r',marker='^')

   # Generate information about the lens-plane images
   xi  = imginfo[:,0]
   yi  = imginfo[:,1]
   mui = imginfo[:,2]
   ri  = n.sqrt(n.abs(mui))*scale
   td  = imginfo[:,3]

   # Sort according to the time delays
   ind = n.argsort(td)
   xx = xi[ind]
   yy = yi[ind]
   tt = td[ind]
   mm = mui[ind]
   rr = ri[ind]

   # Assign labels to the images
   if xx.size == 2:
      labi = ['M1','S1']
   elif xx.size == 4:
      labi = ['M1','M2','S1','S2']
   else:
      labi = n.arange(xx.size) + 1

   # Plot the image locations, with the size of the image (ri) determined
   #  by its magnification.
   plt.figure(1)
   plt.tick_params(labelsize=14)
   ax = plt.axes()
   for i in range(xi.size):
      im = patches.Circle((xi[i],yi[i]),ri[i],color='c',alpha=0.8,
                          label='t(%s): %6.2f' % (labi[i],tt[i]))
      ax.add_patch(im)

   # Make calculations that affect sizes
   xdiff = xi.max() - xi.min()
   ydiff = yi.max() - yi.min()
   dm = max(xdiff,ydiff)

   # Draw arrows between the image locations
   for i in range(xx.size-1):
      xxx = xx[i:]
      yyy = yy[i:]
      ttt = tt[i:]
      for j in range(1,xxx.size):
         dx = xxx[j] - xxx[0]
         dy = yyy[j] - yyy[0]
         plt.arrow(xxx[0],yyy[0],dx,dy,color='k',head_width=0.04*dm,
                   length_includes_head=True)

   xm = xx.mean()
   ym = yy.mean()
   dpl = 1.15
   for i in range(xx.size):
      if xx[i]<xm:
         lsx = -1.
      else:
         lsx = 1.
      if yy[i]<ym:
         lsy = -1.
      else:
         lsy = 1.
      laboff = max((dpl*rr[i], 0.04*dm))
      xl = xx[i] + lsx * laboff
      yl = yy[i] + lsy * laboff
      plt.text(xl,yl,labi[i],ha='center',va='center',fontsize=14)

   # Plot source position, if desired
   if show_src_in_implane:
      plt.plot(xs,ys,'go',ms=10)

   # Print out image info
   print ""
   print "Image     x       y        mu    Delta_t"
   print "----- -------- -------- -------- --------"
   for i in range(xx.size):
      print "%3s   %8.4f %8.4f %8.4f %8.3f" % \
          (labi[i],xx[i],yy[i],mm[i],tt[i])

   # Finally, set plot limits
   imsize = dm * 1.9
   xminp = (xi.max() + xi.min())/2.0 - 0.5*imsize
   xmaxp = xminp + imsize
   yminp = (yi.max() + yi.min())/2.0 - 0.5*imsize
   ymaxp = yminp + imsize
   if xmin is not None:
      xminf = min(xmin,xminp)
   else:
      xminf = xminp
   if ymin is not None:
      yminf = min(ymin,yminp)
   else:
      yminf = yminp
   if xmax is not None:
      xmaxf = max(xmax,xmaxp)
   else:
      xmaxf = xmaxp
   if ymax is not None:
      ymaxf = max(ymax,ymaxp)
   else:
      ymaxf = ymaxp
   plt.xlim(xminf,xmaxf)
   plt.ylim(yminf,ymaxf)
   plt.xlabel(r"$\Delta \alpha$ (arcsec)",fontsize=14)
   plt.ylabel(r"$\Delta \delta$ (arcsec)",fontsize=14)
   if title is not None:
      plt.title(title)
   #plt.legend(markerscale=0.)
   labtxt = ""
   for i in range(tt.size):
      if i < tt.size-1:
         labtxt += "$\Delta$t(%s): %6.2f d\n" % (labi[i],tt[i])
      else:
         labtxt += "$\Delta$t(%s): %6.2f d" % (labi[i],tt[i])
   if dtboxloc == 'tl':
      xloc = 0.03
      yloc = 0.97
      valoc = 'top'
      haloc = 'left'
   elif dtboxloc == 'br':
      xloc = 0.97
      yloc = 0.03
      valoc = 'bottom'
      haloc = 'right'
   else:
      xloc = 0.97
      yloc = 0.97
      valoc = 'top'
      haloc = 'right'
   plt.text(xloc, yloc, labtxt, va=valoc, ha=haloc, fontsize=14,
            transform=ax.transAxes, bbox=dict(facecolor='none'))

#-----------------------------------------------------------------------

def clear_all():
   plt.figure(1)
   plt.clf()
   plt.figure(2)
   plt.clf()

#-----------------------------------------------------------------------

def make_model_plot(lensroot, lensname, dtboxloc='tr', 
                    moddir='/Users/cdf/Projects/Active/Master_database/Models/omega'):
   """
   A top-level function to call plot_crit, plot_lensmod, and plot_points
   in order to make a nice plot.

   For now, assumes that all of the models will have prefixes of the form
    out_[lensroot]_best

   Inputs:
      lensroot  - root name used for file and directory names (e.g., he0435)
      lensname  - full name of lens, used for title (e.g. HE0435-1223)
      moddir    - directory containing all of the glafic models
   """

   # Start by defining the prefix for the input files

   prefix = "%s/%s/out_%s_best" % (moddir,lensroot,lensroot)

   # Make the plots

   plt.figure(1)
   plt.axes().set_aspect('equal')
   plot_crit(prefix)
   plt.figure(1)
   xmin,xmax = plt.xlim()
   ymin,ymax = plt.ylim()
   plot_lensmod(prefix)
   plot_points(prefix,lensname,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,
               dtboxloc=dtboxloc)
