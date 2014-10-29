#!/bin/tcsh
#
#

# Check the command line.
if ( $# < 1 ) then
   echo ""
   echo "*** latex2pdf command line requires at least 1 input. ***"
   echo "    1. Root name of latex file. (i.e., [root].tex)"
   echo "    2. [OPTIONAL] put anything as a second argument to run bibtex"
   echo ""
   echo "*** Exiting. ***"
   echo ""
   exit
endif

# Run latex 3 times on the file (to get all cross-references correct)

latex $1
if ( $# > 1) then
   bibtex $1
endif
latex $1
latex $1

# Run dvips on the file

dvips -t letter -Ppdf -G0 -o $1.ps $1

# Run ps2pdf14 on the postscript file

echo ""
echo "Starting conversion: $1.tex --> $1.pdf"
echo ""
ps2pdf14 $1.ps $1.pdf
echo "Finished converting $1.tex --> $1.pdf"
echo ""


# Clean up and display the final result

echo "Cleaning up...."
echo ""
rm $1.ps
rm *.log *.dvi *.aux *~
echo "Displaying $1.pdf..."
echo ""
#acroread $1.pdf &
open $1.pdf
