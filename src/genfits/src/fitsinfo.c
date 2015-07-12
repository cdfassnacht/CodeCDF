/* fitsinfo.c
 *
 * Usage: fitsinfo filename
 *
 * Description:  Takes an image in FITS format and prints out basic 
 *  information about the file from the header.
 *
 * Revision history:
 * -----------------
 *  2006Sep01 - Chris Fassnacht (CDF) - First working version
 *  2008Jul10 - CDF - Changed main function call to be to the new
 *   print_fits_info, which is more streamlined than the old image_info
 *   function.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "libfits.h"
#include "structdef.h"
#include "fitswcs.h"
#include "fitsfuncs.h"

/*.......................................................................
 *
 * Main program
 *
 */

int main(int argc, char *argv[])
{
  int no_error=1;         /* Flag set to 0 on error */
  char file_name[MAXC];   /* The name of the FITS file to load */
  Image *image=NULL;      /* Fits file container */

  /*
   * Check the command line
   */

  if(argc < 2) {
    fprintf(stderr,"\nfitsinfo.c is a program that prints out some basic");
    fprintf(stderr," information about a fits\n");
    fprintf(stderr," file from the fits header cards.\n");
    fprintf(stderr,"\n   Usage: fitsinfo [fits_file_name]\n\n");
    fprintf(stderr,"For a full printout of all the header cards, use");
    fprintf(stderr," fitshead.c.\n\n");
    return 1;
  }

  /*
   * Get the name of the file.
   */

  strcpy(file_name,argv[1]);

  /*
   * Read the image information.  With the new version of new_Image, this
   *  automatically prints out the information that we want.
   */

  printf("\nReading header from %s.....\n\n",file_name);
  if(!(image = new_Image(file_name))) {
    fprintf(stderr,"\n\nERROR: fitsinfo.  Could not read fits file %s\n\n",
	    file_name);
    return 1;
  }

  /*
   * Clean up and exit
   */

  image = del_Image(image);

  if(no_error) {
    printf("\nProgram fitsinfo.c completed.\n\n");
    return 0;
  }
  else {
    fprintf(stderr,"\nERROR. Exiting fitsinfo.c\n\n");
    return 1;
  }
}

