#!/bin/tcsh
#
#

# Check the command line.
if ( $# < 1 ) then
   echo ""
   echo "*** run_pdflatex command line requires 1 input. ***"
   echo "    1. Root name of latex file. (i.e., [root].tex)"
   echo ""
   echo "*** Exiting. ***"
   echo ""
   exit
endif

# Run pdflatex 3 times on the file to get all cross-references correct

pdflatex $1.tex
pdflatex $1.tex
pdflatex $1.tex

# Clean up and display the final result

echo "Cleaning up...."
echo ""
rm *.log *.dvi *.aux *~ *.bbl *.blg
echo "Displaying $1.pdf..."
echo ""
#acroread $1.pdf &
open $1.pdf
