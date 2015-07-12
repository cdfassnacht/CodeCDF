/* fitsxy2wcs.c
 *
 * Usage: fitsxy2wcs filename
 *
 * Description:  Reads the WCS information from the header cards
 *  in an input fits file and then converts the input (x,y) coordinates
 *  into (RA,Dec) coordinates.
 *
 * Usage: fitsxy2wcs [fits_file] [x] [y]
 *
 * Revision history:
 *  2006Dec19: Chris Fassnacht (CDF), First working version
 *  2007Jun20: Moved WCSinfo definition into fitswcs.c
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "libfits.h"
#include "coords.h"
#include "structdef.h"
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
  int no_error=1;       /* Flag set to 0 on error */
  double inpos[2];      /* Input x and y */
  double outpos[2];     /* Output position, in alpha and delta */
  double dxp,dyp;       /* Pixel offsets between input and reference pixels */
  double dpos[2];       /* Arcsec offsets between input and reference pixels */
  char *headstr;        /* String for header cards */
  char file_name[MAXC]; /* The name of the FITS file to load */
  char *setupname;      /* Name of the optional setup file */
  char line[MAXC];      /* General string for getting input */
  Fits *fits=NULL;      /* Fits file */
  WCSinfo wcs;          /* Container for WCS information */
  Skypos skypos;        /* Outpos in hms format */

  /*
   * Check the command line
   */

  if(argc < 4) {
    fprintf(stderr,"\n");
    fprintf(stderr,
	    "fitsxy2wcs: Converts an input (x,y) coordinate into a (RA,Dec)\n");
    fprintf(stderr," coordinate by using the WCS information in the input\n");
    fprintf(stderr," fits file.\n\n");
    fprintf(stderr,"\nUsage: fitsxy2wcs [filename] [x] [y]\n\n");
    return 1;
  }

  /*
   * Get the name of the file and input x and y.
   */

  strcpy(file_name,argv[1]);
  if((sscanf(argv[2],"%lf",&inpos[0]) != 1)) {
    fprintf(stderr,"ERROR: Bad input value for x\n");
    return 1;
  }
  if((sscanf(argv[3],"%lf",&inpos[1]) != 1)) {
    fprintf(stderr,"ERROR: Bad input value for y\n");
    return 1;
  }

  /*
   * Attempt to open the file.
   */

  fits = new_Fits(file_name, 1, 1, 1, 0);
  if(!fits) {
    fprintf(stderr,"\n\nERROR: fitsxy2wcs.  Could not read fits file\n\n");
    return 1;
  }

  /*
   * Read the WCS information from the fits header
   */

  if(read_wcs_info(fits,&wcs))
    no_error = 0;
  if(no_error) {
    printf("\n");
    print_wcs_info(wcs);
  }

  /*
   * Convert (x,y) into (alpha,delta)
   */

  if(no_error) {
    dxp = inpos[0] - wcs.crpix[0];
    dyp = inpos[1] - wcs.crpix[1];
    for(i=0; i<2; i++) {
      dpos[i] = dxp * wcs.cdcol1[i] + dyp * wcs.cdcol2[i];
    }

    /*
     * Correct for declination effects (crude approach)
     */

    dpos[0] /= cos(PI * wcs.crval[1] / 180.0);
    for(i=0; i<2; i++)
      outpos[i] = wcs.crval[i] + dpos[i];
    deg2spos(outpos[0],outpos[1],&skypos);

    /*
     * Print out results
     */
    printf("# x(pix)   y(pix)   RA(deg)   Dec(deg)      RA            Dec\n");
    printf("%8.2f %8.2f %13.9lf %13.10lf ",
	   inpos[0],inpos[1],outpos[0],outpos[1]);
    print_skypos(stdout,skypos);
    printf("\n");
  }

  /*
   * Clean up and exit
   */

  fits = del_Fits(fits);

  if(no_error) {
    /* printf("\nProgram fitsxy2wcs.c completed.\n\n"); */
    return 0;
  }
  else {
    fprintf(stderr,"\nERROR. Exiting fitsxy2wcs.c\n\n");
    return 1;
  }
}


int MAIN_(void)
{
  return 0;
}


