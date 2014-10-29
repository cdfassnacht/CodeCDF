#!/usr/local/bin/perl
##########################################
# gscbatch.pl is a front end to webquery.pl 
# to return lists of gsc stars from the
# ESO GSC search web page.
#
# Calling sequence and parameter information
# can be found on the SkyView web site:
#     http://skyview.gsfc.nasa.gov/skyview.html
#
##########################################

# Loop over the argument list and look to see
# if file has been defined 
$i=0;
foreach (@ARGV) {
   last if (/^file=/i);
   $i++;
}

# if a filename has been specified, make it the filename
# and splice it out of the argument list
$filename = splice(@ARGV,$i,1);

# Read in the argument list and put the
# values in quotes
$new = 'webquery.pl url="/gsc/gsc" host="archive.eso.org" ' . '\'' . join('\' \'',@ARGV) . '\'';

# Exexcute webquery.pl
@arr = `$new`;

# Put all the blocks togehter
$data = join('',@arr);

# Read past the HTTP header to get to the real data
$locat = index($data, "\n\n") + 2;

# Either write to a file or standard output
if (length($filename)) {
   # get the filename from the argument
   ($key,$truename) = split('=',$filename);

   # Open the file
   open(WWWOUT,">$truename");

   # Print out only from SIMPLE on
   print(WWWOUT substr($data, $locat));
} else {
   # Print to standard outuput
   print substr($data, $locat);
}

# And that's all there is


