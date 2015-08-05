""" 
This code is just an example of how to read in columns from a textfile.
 This can be useful in reading in a SExtractor catalog.  It also involves
 selecting a subset of the rows that satisfy some condition(s).

To execute this program, just type "python example_renata.py [filename]"
where [filename] is the file that you want to read in.

Note that triple quotes denote a comment block.
You can also comment out a line by starting it with the hash sign (#).

Python is a bit odd in that if-then statements, etc., don't have "endif"
 statements or don't define the executable block with {}.  Instead,
 python defines these blocks by the indentation.  

Version history:
 2013_11_18: CDF, Initial version
"""

# At the beginning load the python libraries that you are going to use.
# This program may not use the math library, but it is included to give
#  an example of the "from [library] import [function(s)]" syntax.
# If you don't use a line like this, then you would have to call, for example,
#  sqrt as math.sqrt

import math
from math import fabs,sqrt
import numpy as n  # Use this syntax to create a short-cut name to a library
import sys,os

"""
Here we can define functions.
"""

def read_inputs(file):
   print ""
   print "Reading in %s" % file
   data = n.loadtxt(file)
   # Note that loadtxt has lots of options.  Google "numpy loadtxt" to see.
   
   return data

"""
Here is the main program.  Note that it does not start with a def statement

Start by getting the name of the text file to read from the command line.
#  Note that argv is the list of things that follow "python" on the command
   line.  So if you typed "python read_textfile_example.py CFHT1.cat" then
   argv[0] = read_textfile_example.py
   argv[1] = CFHT1.cat
"""

if len(sys.argv)>1:
   file = sys.argv[1]
else:
   print ""
   print "ERROR: Expected an input filename but there was none."
   print ""
   sys.exit(1)

""" Read in the file.  The output (data) is a 2D array """

data = read_inputs(file)

""" 
Extract the columns that you want.  Remember that the first column is 0.
 Also, python flips the normal array convention for columns and rows, so that
 the first number is the row and the second number is the column.  To
 get all rows associated with the first column, take data[:,0], etc.
 To get all columns associated with the second row, take data[1,:], etc.
"""

id = data[:,0]
x  = data[:,1]
y  = data[:,2]
sgclass = data[:,4]

"""
Select objects that are likely to be galaxies (sgclass<0.5) and also have
You can do this in a 2-step or a 1-step way.  Both are shown below, with
the two-step way first
"""
mask = sgclass<0.5
galx = x[mask]
# Now the 1-step way
galy = y[sgclass<0.5]

"""
You can also select on multiple conditions
"""
mask2 = (sgclass<0.5) & (x<512)
galx_left = x[mask2]
galy_left = y[mask2]

""" 
Make a temporary container for our new variables.  It will only have 2 columns  
Note that galx.shape[0], used below, gives the number of rows in the galx vector
"""

outdat = n.zeros((galx.shape[0],2))
outdat[:,0] = galx
outdat[:,1] = galy

""" Write the output to a file """

outfile = 'foo.txt'
print "Writing output file %s" % outfile
n.savetxt(outfile,outdat,fmt='%8.2f')
