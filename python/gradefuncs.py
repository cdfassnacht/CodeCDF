""" 
Functions that are useful in plotting a grade histogram 

"""

import numpy as np
from matplotlib import pyplot as plt
from astropy.io import ascii

#---------------------------------------------------------------------------

def read_table(infile, colname):
   """
   This new code can be used to read both the old-school text files (if they
    had the two-column format - see help for the read_text function) but
    also (and more importantly) the information directly from a CSV table
    of the form that is exported from smartsite or canvas.

   Inputs:
    infile  - input file name
    colname - the name of the column containing the score of interest.
              NOTE: for the old-school text files, this will be 'col2'
              while for the CSV files it could be something like
              'Midterm 2 (32620)' or 'MT2' or 'Final Score'
   """

   """ Read in the table """
   try:
      tab = ascii.read(infile, guess=False, format='csv')
   except:
      tab = ascii.read(infile)
   print(tab.colnames)

   """ Get the relevant information """
   try:
      tot = tab[colname].copy()
   except KeyError:
      print ''
      print 'Could not find a column matching %s in %s' % (colname,infile)
      tot = None
   return tot

#---------------------------------------------------------------------------

def read_text(infile, nscorecols=1):
   """
   Function to read in the scores from the old-school text files that were
   created by modifying the csv files that came from downloading the gradebook
   from smartsite or canvas.

   There are two expected input formats:
     Two-column, designated by setting nscorecols=1
       Name  total_score
     Three-column, designated by setting nscorecols=2
       Name  multiple_choice_score  short_answer_score

   The old code used the numpy loadtxt function to load the data, but this
   new code uses the astropy.io.ascii read function
   """

   """ Read the data into an astropy Table structure """
   tab = ascii.read(infile)

   """ Generate the total score array """
   if nscorecols == 1:
      tot = tab['col2'].copy()
   else:
      tot = tab['col2'] + tab['col3']

   return tot

#---------------------------------------------------------------------------

def plot_tothist(infile, tot, maxy, binsize=3):
   """
   Plot the total-score histogram, where the total score (tot) has been
   previous calculated or read-in by the input functions
   """

   """ Calculate moments of the distribution """
   mn = tot.mean()
   med = np.median(tot)
   mp = tot.mean() + tot.std()
   mm = tot.mean() - tot.std()

   """ Report on the properties of the distibution """
   print ""
   print "Statistics for %s" % infile
   print "---------------------------------"
   print "  Mean:         %5.1f" % mn
   print "  Median:       %5.1f" % med
   print "  Sigma:        %5.1f" % tot.std()
   print "  Mean - 1 sig: %5.1f" % mm
   print "  Mean + 1 sig: %5.1f" % mp
   print ""

   """ Plot the distribution """
   binhist = range(int(tot.min())-1,int(tot.max())+3,binsize)
   plt.hist(tot,binhist,histtype='step',ec='k')
   plt.ylim(0,maxy)
   plt.axvline(x=mn, ymin=0, ymax=maxy, c='r', lw=3)
   plt.axvline(x=mm, ymin=0, ymax=maxy, c='b', lw=3)
   plt.axvline(x=mp, ymin=0, ymax=maxy, c='b', lw=3)
   plt.title("Distribution of scores for %s" % infile)
   plt.xlabel("Scores")
   plt.ylabel("N")
   plt.show()

