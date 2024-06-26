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
      print('')
      print('Could not find a column matching %s in %s' % (colname,infile))
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

def plot_tothist(infile, tot, maxy, binsize=3, colname=None):
   """
   Plot the total-score histogram, where the total score (tot) has been
   previously calculated or read-in by the input functions
   """

   """ Calculate moments of the distribution """
   mn = tot.mean()
   med = np.median(tot)
   mp = tot.mean() + tot.std()
   mm = tot.mean() - tot.std()
   tot.sort()
   m25 = tot[int(0.25*tot.size)]
   m75 = tot[int(0.75*tot.size)]
   m16 = tot[int(0.16*tot.size)]
   m84 = tot[int(0.84*tot.size)]
   m80 = tot[int(0.80*tot.size)]

   """ Report on the properties of the distibution """
   print('')
   print("Statistics for %s" % infile)
   print("---------------------------------")
   print("  Mean:                 %5.1f" % mn)
   print("  Sigma:                %5.1f" % tot.std())
   print("  Mean - 1 sig:         %5.1f" % mm)
   print("  Mean + 1 sig:         %5.1f" % mp)
   print("  Median:               %5.1f" % med)
   print('  Lowest quartile:      %5.1f' % m25)
   print('  Highest quartile:     %5.1f' % m75)
   print('  Low end of 68%% range: %5.1f' % m16)
   print('  Top end of 68%% range: %5.1f' % m84)
   print('  Top 20%%:              %5.1f' % m80)
   print('  Maximum:              %5.1f' % tot.max())
   print('  Minimum:              %5.1f' % tot.min())
   print('')

   """ Plot the distribution """
   binhist = range(int(tot.min())-1,int(tot.max())+3,binsize)
   plt.hist(tot,binhist,histtype='step',ec='k')
   plt.ylim(0,maxy)
   plt.axvline(x=med, ymin=0, ymax=maxy, c='r', lw=3)
   plt.axvline(x=m16, ymin=0, ymax=maxy, c='b', lw=3)
   plt.axvline(x=m84, ymin=0, ymax=maxy, c='b', lw=3)
   if colname is not None:
      plt.title("Distribution of scores for %s" % colname)
   else:
      plt.title("Distribution of scores for %s" % infile)
   plt.xlabel("Scores")
   plt.ylabel("N")
   plt.show()

