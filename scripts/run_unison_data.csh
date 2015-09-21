#!/bin/tcsh
#
# Runs unison to back up selected data directories
#

foreach i (Data/Lenses Data/SHARP)
    run_unison.csh $i
    end
