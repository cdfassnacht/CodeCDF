""" 
A python program to plot a grade histogram 

Usage: python plot_grade_hist.py [filename] [nscorecols] [maxy]

Inputs:
 filename   - text file containing either two columns (name total) or
              three columns (name, score_multiple-choice, score_short-answer)
 nscorecols - number of score columns (1 for 2-column input, 2 for 3-column
              input).
 maxy       - maximum value for y axis 
"""

import numpy as n
from matplotlib import pyplot as p
import sys

if len(sys.argv) < 4:
   print ''
   print 'ERROR: This program requires 3 input parameters:'
   print '   1. infile   - name of the input file containing scores'
   print '   2. ncolumns - number of columns with scores (1 or 2)'
   print '   3. maxy     - maximum y value for plot'
   print ''
   print 'Format: python plot_grade_hist.py infile ncolumns maxy'
   print ''
   sys.exit()

infile = sys.argv[1]
nscorecols = int(sys.argv[2])
maxy = float(sys.argv[3])
binsize = 3

namecol = 0
mccol  = 1
sacol  = 2
totcol = 1

print ""
print "Reading %d columns of data from %s" % (nscorecols,infile)

if nscorecols == 1:
   tot = n.loadtxt(infile,usecols=(totcol,),unpack=True)
else:
   mc,sa = n.loadtxt(infile,usecols=(mccol,sacol),unpack=True)
   tot = mc + sa

mn = tot.mean()
med = n.median(tot)
mp = tot.mean() + tot.std()
mm = tot.mean() - tot.std()

print ""
print "Statistics for %s" % infile
print "---------------------------------"
print "  Mean:         %5.1f" % mn
print "  Median:       %5.1f" % med
print "  Sigma:        %5.1f" % tot.std()
print "  Mean - 1 sig: %5.1f" % mm
print "  Mean + 1 sig: %5.1f" % mp
print ""

binhist = range(int(tot.min())-1,int(tot.max())+3,binsize)
p.hist(tot,binhist,histtype='step',ec='k')
p.ylim(0,maxy)
p.axvline(x=mn, ymin=0, ymax=maxy, c='r', lw=3)
p.axvline(x=mm, ymin=0, ymax=maxy, c='b', lw=3)
p.axvline(x=mp, ymin=0, ymax=maxy, c='b', lw=3)
p.title("Distribution of scores for %s" % infile)
p.xlabel("Scores")
p.ylabel("N")
p.show()

#totdat = gradedat[:,0:1].copy()
#totdat[:,1] = tot.copy()
#print totdat
