"""
example_2.py - This is just the same as example.py, but is meant to be
run from within an interactive python environment such as ipython, as
opposed to example.py, which is meant to be run from the command line.

To run this, first start up ipython with the pylab capability.  Then
type the following two lines:

  import example_2
  example_2.run_the_example()

This is an example python program.  The triple quotation marks indicate
a comment block.  It is good practice, especially for code that you submit
in this class, to have lots of documentation comments in your code.
"""

# Comments on a single line can also be designated by the leading pound sign
#  (hash mark).


"""
Start the code by importing the libraries that you are going to need.
For this class, numpy, scipy, pyfits, and matplotlib.pyplot are
probably the most common libraries that we will use.  You can also
write your own libraries.  You can import a library (import numpy),
import a function or sub-library from a library (from numpy import
sqrt), and import a library while giving it a nickname (import numpy
as n).
"""

import numpy
from matplotlib import pyplot as plt

#-----------------------------------------------------------------------

"""
You may want to define functions that will do important things.  Here is
a function definition (of a somewhat silly function).  Note that in python
the indentation tells you what is associated with what.  Thus, you should
make sure that everything that belongs together (as in a function definition,
a for loop, an if statement, etc., has the same indentation.  The end of
the indentation indicates the end of the block
"""

def plot_my_data(x,y,verbose=True):
   """
   This function will plot y vs x.  It has two required parameters, x and y.
   It also has an optional parameter, called 'verbose'.  The default value
   of the verbose parameter is True.  Thus, if you want verbose to be true,
   then you don't have to pass it to the function.  All optional parameters
   must have a default value.  This can be a logical (True or False), an
   integer, a float, a string, etc.
   """

   # Note the indentation for the if block
   if (verbose):
      print ""
      print "Now plotting the data."

   # In the plot command below, the 'bo' tells the program to plot blue ('b')
   #  circles ('o')
   print ""
   #print "Close the plot window to continue..."
   plt.plot(x,y,'bo')

#-----------------------------------------------------------------------



"""
Here is the main part of the example, which for an interactive environment
is also defined as a function


One very useful task is that which reads in a text file containing data
(as floating point numbers) in columns.  The data are put into a 2-dimensional
numpy array.  Note that in this example, polyfit_data.txt contains 50 lines
of data and 2 columns.  Therefore, in python-speak, it is a 50x2 array.
"""

def run_the_example():

   infile = 'polyfit_data.txt'
   data = numpy.loadtxt(infile)

   print ""
   print "Dimensions of input data are %d x %d" % (data.shape[0],data.shape[1])

   """
   Once the data are loaded in, you can separate the columns.  Note that in
   python, (1) arrays are zero-indexed, i.e., the first element of an array
   has an index of 0, and not an index of 1, and (2) for 2-d arrays, the first
   index refers to y and the second refers to x.
   """

   x = data[:,0]  # Extracts the first column
   y = data[:,1]  # Extracts the second column

   """
   In the above definitions, x and y are also numpy arrays, having sizes 50x1.
   numpy arrays have all sorts of nice functions associated with them, such
   as mean and standard deviation, but not median (which is a separate numpy
   function).
   """

   ymean = y.mean()
   yrms = y.std()
   ymed = numpy.median(y)
   
   print ""
   print "The y vector has mean=%f, rms=%f, and median=%f." % (ymean,yrms,ymed)


   """
   Note that for data that are more than one dimension, you can take the overall
   mean or the mean of the columns or of the rows.  You can do those things
   by using the optional axis parameter.
   """

   meanall = data.mean()
   meancols = data.mean(axis=0)
   meanrows = data.mean(axis=1)

   print "The overall data set has mean = %5.2f" % meanall
   print "The means of the columns are:",meancols

   """
   You can create new arrays from existing arrays by testing for some
   condition.
   """

   # The following line creates a mask that is true where y is greater than 1.5
   # Note that this mask is just an array of True or False that has the same
   #  size as y.  Therefore, you can use it to select the corresponding members
   #  of the x array as well.
   mask = y>1.5
   print ""
   print "mask = ",mask
   newx1 = x[mask]
   newy1 = y[mask]
   print ""
   print "newy1 = ",newy1

   # You could have taken a shortcut for the above step.  For example:

   newy2 = y[y>1.5]
   print "newy2 = ",newy2

   """
   Finally, call the function to plot the data
   """

   plt.figure(1)
   plot_my_data(x,y)
   plt.figure(2)
   plot_my_data(newx1,newy1)
   plt.show()

