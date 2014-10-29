#!/usr/bin/sh
echo "Enter the name of the file containing the source names:  "
read filename
mapscript_2 $filename
chmod 755 automap_*.sh
