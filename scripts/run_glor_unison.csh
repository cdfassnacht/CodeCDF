#!/bin/tcsh
#
# Runs unison to back up the requested directory
#

# Check the command line.
if ( $# < 1 ) then
   echo ""
   echo "run_unison"
   echo "  Usage: run_unison directory_name_to_be_backed_up"
   echo "  Description: runs unison on the selected directory."
   echo ""
   echo "*** Exiting. ***"
   echo ""
   exit
endif

echo
echo
echo Running unison on directory: $1
echo
unison -ui text $1 ssh://glorfindel.physics.ucdavis.edu/$1
echo
echo ------------------------------------------------
