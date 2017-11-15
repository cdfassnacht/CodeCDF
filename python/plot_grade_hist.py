""" 
A python program to plot a grade histogram 

Usage: python plot_grade_hist.py [filename] [nscorecols] [maxy]

Required inputs:
 filename   - text file containing either two columns (name total) or
              three columns (name, score_multiple-choice, score_short-answer)
 colname    - the name of the column containing the score of interest.
              NOTE: for the old-school text files, this will be 'col2'
               _unless_ the old-school file also is in 3-column format, in which
               case this parameter is ignored and the optional nscorecols
               parameter should be set to 2.
              If the input file is in CSV format, the colname parameter could
              be something like 'Midterm 2 (32620)' or 'MT2' or 'Final Grade'
 maxy       - maximum value for y axis 

Optional input:
 nscorecols - number of score columns (1 for 2-column input, 2 for 3-column
              input). ONLY set this if the input file is in the old-school
              text format AND it is in 3-column format (i.e., with 
              nscorecols=2).  If it is in the old-school text format but is
              in the 2-column input, then DO NOT set this keyword, but just
              set the colname variable above to 'col2'
"""

import numpy as n
from matplotlib import pyplot as p
import sys
import gradefuncs as gf

if len(sys.argv) < 4:
   print ''
   print 'ERROR: This program requires at least 3 input parameters:'
   print '   1. infile   - name of the input file containing scores'
   print '   2. colname  - name of column containing the relevant score if the'
   print '       input file is in csv format produced by smartsite or canvas or'
   print '       if it is in old-school text format with one total-score column'
   print '       In the second case (text format with one column of scores) the'
   print '       colname parameter should be set to "col2"'
   print '   3. maxy     - maximum y value for plot'
   print 'It may also take an optional fourth parameter, which should ONLY BE'
   print ' USED if the file is BOTH in the old-school text format and has'
   print ' two columns with scores (one for multiple-choice and one for short'
   print ' answer), in which case, this parameter should be used and set to 2.'
   print ''
   print 'Format: python plot_grade_hist.py infile colname maxy'
   print '               --- or ---'
   print 'Format: python plot_grade_hist.py infile colname maxy 2'
   print ''
   sys.exit()

if len(sys.argv) == 5:
   old_3col = True
else:
   old_3col = False

infile = sys.argv[1]
colname = sys.argv[2]
maxy = float(sys.argv[3])

if old_3col:
   tot = gf.read_text(infile,2)
else:
   tot = gf.read_table(infile, colname)

if tot is None:
   print 'Could not plot histogram'
   print ''
   sys.exit()

binsize = 3
gf.plot_tothist(infile,tot,maxy,binsize)

