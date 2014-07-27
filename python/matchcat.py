"""
A set of functions to match sources in two input catalogs
"""

import numpy as n
import coords
from matplotlib import pyplot as plt

#------------------------------------------------------------------------------

class matchcat:

   """ A class contaiing matching methods """

   def __init__(self, catfile1, catfile2):
      """
      Sets up the matchcat container
      """

      self.catfile1 = catfile1
      self.catfile2 = catfile2

   """ TO BE CONTINUED """

#-----------------------------------------------------------------------------

def match_coords(ra1, dec1, ra2, dec2, rmatch, dra2=0., ddec2=0., doplot=True):
   """
   The main function to match coordinates.  
   Inputs:
      ra1      - RA (decimal degrees) for first catalog
      dec1     - Dec (decimal degrees) for first catalog
      ra2      - RA (decimal degrees) for second catalog
      dec2     - Dec (decimal degrees) for second catalog
      rmatch   - max distance for a valid match (arcsec)
      dra2     - optional offset in ARCSEC to apply to ra2, if there is a known
                 offset between the catalogs (default=0.0)
      ddec2    - optional offset in ARCSEC to apply to dec2, if there is a known
                 offset between the catalogs (default=0.0)
   """

   print ""
   print "Catalog info"
   print "--------------------------------------------"
   print " Catalog 1: %d coordinates" % ra1.size
   print " Catalog 2: %d coordinates" % ra2.size

   """ Initialize containers for output information """
   nmatch   = n.zeros(ra1.size,dtype=int)
   dxmatch  = n.zeros(ra1.size)
   dymatch  = n.zeros(ra1.size)
   ramatch  = n.zeros(ra1.size)
   decmatch = n.zeros(ra1.size)
   indmatch = n.ones(ra1.size,dtype=int) * -1

   """ Correct for known shifts """
   #dec2 += 1.39e-4 temporary kludge for fixing Cl1604 matches
   ra2 = ra2.copy() + dra2/(3600.*n.cos(dec2))
   dec2 = dec2.copy() + ddec2/3600.

   """ Loop over catalog1 """
   print ""
   print "Searching for matches..."
   print "------------------------------"
   for i in range(ra1.size):
      dx,dy = coords.sky_to_darcsec(ra1[i],dec1[i],ra2,dec2)
      dpos = n.sqrt(dx**2 + dy**2)
      isort = n.argsort(dpos)
      if dpos[isort[0]]<=rmatch:
         dxmatch[i]  = dx[isort[0]]
         dymatch[i]  = dy[isort[0]]
         ramatch[i]  = ra1[i]
         decmatch[i] = dec1[i]
         nmatch[i] = dpos[dpos<=rmatch].size
         indmatch[i] = isort[0]
      del dx,dy,dpos
   print " Sources in catalog 1 with matches in catalog 2:  %d" % \
       (nmatch>0).sum()
   mdx = dxmatch[nmatch>0]
   mdy = dymatch[nmatch>0]
   mra = ramatch[nmatch>0]
   mdec = decmatch[nmatch>0]
   print " Median offset for matches (RA):  %+6.2f arcsec" % n.median(mdx)
   print " Median offset for matches (Dec): %+6.2f arcsec" % n.median(mdy)

   """ Plot up some offsets, if desired """
   if doplot:
      plt.figure(1)
      plt.scatter(mdx,mdy)
      plt.axis('scaled')
      plt.xlabel(r'$\Delta \alpha$ (arcsec)')
      plt.ylabel(r'$\Delta \delta$ (arcsec)')
      plt.title('Offsets between matched sources (rmatch = %5.2f)' % rmatch)
      plt.axvline(0.0,color='r')
      plt.axhline(0.0,color='r')
      plt.xlim(-1.1*rmatch,1.1*rmatch)
      plt.ylim(-1.1*rmatch,1.1*rmatch)

      plt.figure(2)
      #
      ax1 = plt.subplot(221)
      plt.scatter(mra,mdy)
      plt.setp(ax1.get_xticklabels(), visible=False)
      plt.ylabel(r'$\Delta \delta$ (arcsec)')
      plt.axhline(0.0,color='r')
      #
      ax2 = plt.subplot(223, sharex=ax1)
      plt.scatter(mra,mdx)
      plt.xlabel(r'$\alpha$')
      plt.ylabel(r'$\Delta \alpha$ (arcsec)')
      plt.axhline(0.0,color='r')
      #
      ax3 = plt.subplot(222, sharey=ax1)
      plt.scatter(mdec,mdy)
      plt.axhline(0.0,color='r')
      plt.setp(ax3.get_xticklabels(), visible=False)
      plt.setp(ax3.get_yticklabels(), visible=False)
      #
      ax4 = plt.subplot(224)
      plt.scatter(mdec,mdx)
      plt.xlabel(r'$\delta$')
      plt.axhline(0.0,color='r')
      plt.setp(ax4.get_yticklabels(), visible=False)

   """ Clean up """
   #del ra1,dec1,ra2,dec2
   del ramatch,decmatch

   return indmatch,nmatch,dxmatch,dymatch

#------------------------------------------------------------------------------

def match_xy(x1, y1, x2, y2, rmatch, dx2=0., dy2=0., doplot=True):
   """
   The main function to match coordinates.  
   Inputs:
      x1       - x coordinate in first catalog
      y1       - y coordinate in first catalog
      x2       - x coordinate in second catalog
      y2       - y coordinate in second catalog
      rmatch   - max distance for a valid match (pixels)
      dx2      - optional offset to apply to x2, if there is a known offset
                 between the catalogs (default=0.0)
      dy2      - optional offset to apply to y2, if there is a known offset
                 between the catalogs (default=0.0)
   """

   print ""
   print "Catalog info"
   print "--------------------------------------------"
   print " Catalog 1: %d coordinates" % x1.size
   print " Catalog 2: %d coordinates" % x2.size

   """ Initialize containers for output information """
   nmatch   = n.zeros(x1.size,dtype=int)
   dxmatch  = n.zeros(x1.size)
   dymatch  = n.zeros(x1.size)
   xmatch  = n.zeros(x1.size)
   ymatch = n.zeros(x1.size)
   indmatch = n.ones(x1.size,dtype=int) * -1

   """ Correct for known offsets """
   x2 = x2.copy() + dx2
   y2 = y2.copy() + dy2

   """ Loop over catalog1 """
   print ""
   print "Searching for matches..."
   print "------------------------------"
   for i in range(x1.size):
      dx = x2 - x1[i]
      dy = y2 - y1[i]
      dpos = n.sqrt(dx**2 + dy**2)
      isort = n.argsort(dpos)
      if dpos[isort[0]]<=rmatch:
         dxmatch[i]  = dx[isort[0]]
         dymatch[i]  = dy[isort[0]]
         xmatch[i]   = x1[i]
         ymatch[i]   = y1[i]
         nmatch[i]   = dpos[dpos<=rmatch].size
         indmatch[i] = isort[0]
      del dx,dy,dpos
   print " Sources in catalog 1 with matches in catalog 2:  %d" % \
       (nmatch>0).sum()
   mdx = dxmatch[nmatch>0]
   mdy = dymatch[nmatch>0]
   mx = xmatch[nmatch>0]
   my = ymatch[nmatch>0]
   print " Median offset for matches (X):  %+6.2f pixels" % n.median(mdx)
   print " Median offset for matches (Y): %+6.2f pixels" % n.median(mdy)

   """ Plot up some offsets, if desired """
   if doplot:
      plt.figure(1)
      plt.scatter(mdx,mdy)
      plt.xlabel(r'$\Delta x$ (pixels)')
      plt.ylabel(r'$\Delta y$ (pixels)')
      plt.title('Offsets between matched sources (rmatch = %5.2f)' % rmatch)
      plt.axvline(0.0,color='r')
      plt.axhline(0.0,color='r')
      plt.xlim(-1.1*rmatch,1.1*rmatch)
      plt.ylim(-1.1*rmatch,1.1*rmatch)

      plt.figure(2)
      plt.subplot(221)
      plt.scatter(mx,mdy)
      plt.ylabel(r'$\Delta y$ (pixels)')
      plt.axhline(0.0,color='r')
      plt.subplot(223)
      plt.scatter(mx,mdx)
      plt.xlabel('x')
      plt.ylabel(r'$\Delta x$ (pixels)')
      plt.axhline(0.0,color='r')
      plt.subplot(222)
      plt.scatter(my,mdy)
      plt.ylabel(r'$\Delta y$ (pixels)')
      plt.axhline(0.0,color='r')
      plt.subplot(224)
      plt.scatter(my,mdx)
      plt.xlabel('y')
      plt.ylabel('$\Delta x$ (pixels)')
      plt.axhline(0.0,color='r')

   """ Clean up """
   del x1,y1,x2,y2
   del xmatch,ymatch

   return indmatch,nmatch,dxmatch,dymatch

#-----------------------------------------------------------------------

def find_match(catfile1, catfile2, rmatch, catformat1='ascii', 
               catformat2='ascii', 
               racol1=None, deccol1=None, racol2=None, deccol2=None,
               namecol1=None, namecol2=None, dra2=0., ddec2=0., doplot=True):
   """
   The main function to match catalogs contained in two input files.  The
   input files are expected (for now) to have RA and Dec in decimal degrees
   Inputs:
      catfile1    - file containing the first catalog
      catfile2    - file containing the second catalog
      rmatch      - max distance for a valid match
      catformat1  - Format for first catalog file: 'ascii' (default) or 'ldac'
      catformat2  - Format for second catalog file: 'ascii' (default) or 'ldac'
      racol1      - column containing RA in the first file
      deccol1     - column containing Dec in the first file
      racol2      - column containing RA in the second file
      deccol2     - column containing Dec in the second file
      dra2        - Offset in ARCSEC to add to ra2 in case of known shifts 
                    between cats
      ddec2       - Offset in ARCSEC to add to dec2 in case of known shifts
                    between cats
   """

   """ Read inputs """
   import astrom_simple as astsimp
   try:
      cat1 = astsimp.Secat(catfile1,catformat=catformat1,racol=racol1,
                           deccol=deccol1,namecol=namecol1)
      cat1.get_radec()
   except:
      print ""
      print "ERROR: Could not read RA and Dec from %s" % catfile1
      return
   try:
      cat2 = astsimp.Secat(catfile2,catformat=catformat2,racol=racol2,
                           deccol=deccol2,namecol=namecol2)
      cat2.get_radec()
   except:
      print ""
      print "ERROR: Could not read RA and Dec from %s" % catfile2
      return
   print ""
   
   #dec2 += 1.39e-4 temporary kludge for fixing Cl1604 matches
   
   """ Do the matching """
   cat1.indmatch,cat1.nmatch,cat1.matchdx,cat1.matchdy = \
       match_coords(cat1.ra,cat1.dec,cat2.ra,cat2.dec,
                    rmatch,dra2,ddec2,doplot)
   
   return cat1,cat2

#--------------------------------------------------------------------------

def write_matchcat(cat1,cat2,outfile,rmatch,c1fluxcol,c2fluxcol):
   """
   Writes an output file in the format of the file produced by catcomb.c.
   *** NOTE *** The catalog matching (with the find_match function in
     matchcat.py) has to have been run before running this code.
   
   Inputs:
      cat1    - first catalog used in the matching
      cat2    - second catalog used in the matching (needs to be fixed for more)
      outfile - output file
      rmatch  - match radius used for the matching - used only to put info
                into the output file.
   """

   """ Get the mask for the matches """
   mask = cat1.indmatch>-1

   """ Get info on the matched objects """
   c1d  = cat1.data
   c1id = n.arange(1,c1d.size+1)
   c1mi = cat1.indmatch.copy()
   indm = cat1.indmatch[mask]
   ra1  = cat1.ra
   dec1 = cat1.dec[mask]
   ct1  = (n.arange(1,cat1.data.size+1))[mask]
   c1m  = cat1.data[mask]
   c2m  = cat2.data[indm]
   dx   = cat1.matchdx[mask]
   dy   = cat1.matchdy[mask]
   dpos = n.sqrt(dx**2 + dy**2)

   """ Write match info to output file """
   ofile = open(outfile,'w')
   #
   # Need to fix format here to match matchcat.c
   #
   for i in range(cat1.ra.size):
      c1dat = cat1.data[i]
      c1flux = c1dat['f%d'% c1fluxcol] 
      if c1mi[i]>-1:
         c2dat = cat2.data[c1mi[i]]
         ind2  = c2dat['f0']
         flux2 = c2dat['f%d' % c2fluxcol]
         dpos  = n.sqrt(cat1.matchdx[i]**2 + cat1.matchdy[i]**2)
      else:
         ind2  = 0
         flux2 = 0.
         dpos  = 0.
      ofile.write('%05d %11.7f %+11.7f     0.00     0.00 %5d     0.00 %7.4f ' % \
                     (c1id[i],cat1.ra[i],cat1.dec[i],c1id[i],c1flux))
      ofile.write('%5d %8.2f %6.2f' % (ind2,dpos,flux2))
      ofile.write('\n')
   ofile.close()

