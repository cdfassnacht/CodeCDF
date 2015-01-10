#!/bin/tcsh
#
# Runs unison to back up selected data directories
#

foreach i (Data2/Lenses Data2/SHARP)
    run_unison.csh $i
    end
