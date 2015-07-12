/* fitswcs2xy.c
 *
 * Usage: fitswcs2xy filename
 *
 * Description:  Takes as inputs two files, a fits file and a text file
 *  containing a list of (RA,Dec) coordinates.  Givenn the WCS information
 *  in the fits header, converts the coordinates to the appropriate (x,y)
 *  positions in the image.
 *
 * Usage: fitswcs2xy [fits_file] [catfile] [catformat] [output_file]
 *
 * Revision history:
 *  2008Jul07: Chris Fassnacht (CDF), First working version
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "structdef.h"
#include "libfits.h"
#include "dataio.h"
#include "coords.h"
#include "fitswcs.h"

/*.......................................................................
 *
 * Function declaration
 *
 */


/*.......................................................................
 *
 * Main program
 *
 */

int main(int argc, char *argv[])
{
  int i;                /* Looping variable */
  int catformat;        /* Format of input catalog */
  int ncat=0;           /* Number of members of the catalog */
  int no_error=1;       /* Flag set to 0 on error */
  char fitsname[MAXC];  /* The name of the FITS file to load */
  char catname[MAXC];   /* The name of the catalog file to load */
  char outname[MAXC];   /* The name of the output file */
  char line[MAXC];      /* General string for getting input */
  Fits *fits=NULL;      /* Fits file */
  Secat *secat=NULL;    /* Catalog container */

  /*
   * Check the command line
   */

  if(argc < 5) {
    fprintf(stderr,"\n");
    fprintf(stderr,
	    "fitswcs2xy: Given a fits image and a catalog containing \n");
    fprintf(stderr,
	    " (RA,Dec) coordinates, either in HMS or decimal degrees,\n");
    fprintf(stderr,
	    " converts those coordinates into (x,y) positions on the image\n");
    fprintf(stderr," by using the WCS information in the input\n");
    fprintf(stderr," fits file.\n\n");
    fprintf(stderr,
	    "Usage: fitswcs2xy [fitsfile] [catfile] [catformat] [outfile]\n");
    fprintf(stderr,"\nHit return to see the options for catformat: ");
    fgets(line,MAXC,stdin);
    secat_format();
    return 1;
  }

  /*
   * Get the name of the input files, input format, and output file
   */

  strcpy(fitsname,argv[1]);
  strcpy(catname,argv[2]);
  strcpy(outname,argv[4]);

  if((sscanf(argv[3],"%d",&catformat) != 1)) {
    fprintf(stderr,"ERROR: Bad input value for catformat\n\n");
    fprintf(stderr,"Hit return to see the valid options for catformat: ");
    fgets(line,MAXC,stdin);
    secat_format();
    return 1;
  }

  /*
   * Attempt to open the fits file file.
   */

  fits = new_Fits(fitsname, 1, 1, 1, 0);
  if(!fits) {
    fprintf(stderr,"\n\nERROR: fitswcs2xy.  Could not read fits file\n\n");
    return 1;
  }

  /*
   * Read the catalog.  read_secat should automatically convert any
   *  RA,Dec positions to decimal degrees, if they are not already in that
   *  format
   */

  if(!(secat = read_secat(catname,'#',&ncat,catformat)))
    no_error = 0;

  /*
   * Convert the catalog (RA,Dec) positions to (x,y) positions
   */

  if(no_error)
    if(radec_to_ccdxy(fits,secat,ncat) == ERROR)
      no_error = 0;

  /*
   * Print out results
   */

  if(no_error)
    if(write_secat(secat,ncat,outname,0))
      no_error = 0;

  /*
   * Clean up and exit
   */

  fits = del_Fits(fits);
  secat = del_secat(secat);

  if(no_error) {
    /* printf("\nProgram fitswcs2xy.c completed.\n\n"); */
    return 0;
  }
  else {
    fprintf(stderr,"\nERROR. Exiting fitswcs2xy.c\n\n");
    return 1;
  }
}


int MAIN_(void)
{
  return 0;
}


