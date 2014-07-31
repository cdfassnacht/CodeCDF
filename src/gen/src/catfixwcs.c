/*
 * catfixwcs.c
 *
 * Usage: catfixwcs [catalog_filename] [catalog_format] [output_filename]
 *
 * This program reads in a catalog and then fixes the RA and Dec values
 *  by adding some delta_RA and delta_Dec offsets.
 *
 * 01Aug2007 CDF, First working version
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "structdef.h"
#include "coords.h"
#include "dataio.h"
#include "catlib.h"


/*.......................................................................
 *
 * Main program
 *
 */

int main(int argc, char *argv[])
{
  int i;                /* Looping variable */
  int no_error=1;       /* Flag set to 0 on error */
  int ncat;             /* Number of lines in the input catalog */
  int format;           /* Format of input/output files */
  double alpha0,delta0; /* Temporary holders for positions, in _radians_ */
  double alpha,delta;   /* Temporary holders for positions, in _radians_ */
  char line[MAXC];      /* Generic string variable */
  char catfile[MAXC];   /* Filename for input catalog */
  char outfile[MAXC];   /* Filename for distcalc-like output file */
  Pos offset;           /* Offset in arcsec to add to the coordinates */
  Secat *secat=NULL;    /* Data array from catalog */
  Secat *sptr;          /* Pointer to navigate secat */

  /*
   * Check the command line invocation
   */

  if(argc < 3) {
    fprintf(stderr,"\nUsage: catfixwcs [catfile] [format] [outfile]\n\n");
    fprintf(stderr," The format flag indicates the format of the");
    fprintf(stderr," input catalog file:\n\n");
    fprintf(stderr," Hit return to see formats.\n");
    fgets(line,MAXC,stdin);
    secat_format();
    return 1;
  }
  printf("\n");

  /*
   * Get the inputs from the command line.
   */

  strcpy(catfile,argv[1]);
  strcpy(outfile,argv[3]);
  if(sscanf(argv[2],"%d",&format) != 1) {
    fprintf(stderr,"ERROR: Bad value for format.\n");
    return 1;
  }

  /*
   * Fill the data structure
   */

  if(!(secat = read_secat(catfile,'#',&ncat,format))) {
    fprintf(stderr,"ERROR.  Exiting catfixwcs.c\n");
    secat = del_secat(secat);
    return 1;
  }

  /*
   * Get the offsets
   */

  printf("\nThis program takes an input catalog and applies an offset to \n");
  printf(" the RA, Dec coordinates in the catalog.\n\n");
  printf("Enter offsets to apply, in ARCSEC, in the format ");
  printf("delta_RA delta_Dec: ");
  fgets(line,MAXC,stdin);
  while(sscanf(line,"%lf %lf",&offset.x,&offset.y) != 2) {
    fprintf(stderr,"ERROR: Bad format for offsets.  Enter offsets again: ");
    fgets(line,MAXC,stdin);
  }

  /*
   * Apply the offsets to the RA, Dec coordinates
   */

  for(i=0,sptr=secat; i<ncat; i++,sptr++) {
    alpha0 = sptr->alpha * PI / 180.0;
    delta0 = sptr->delta * PI / 180.0;
    if(offset2rad(alpha0,delta0,offset,&alpha,&delta))
      no_error = 0;
    sptr->alpha = alpha * 180.0 / PI;
    sptr->delta = delta * 180.0 / PI;
    deg2spos(sptr->alpha,sptr->delta,&sptr->skypos);
  }

  /*
   * Write output file
   */

  if(no_error)
    if(write_secat(secat,ncat,outfile,format))
      no_error = 0;

  /*
   * Clean up and exit
   */

  secat = del_secat(secat);

  if(no_error) {
    printf("\nProgram catfixwcs finished.\n\n");
    return 0;
  }
  else {
    fprintf(stderr,"\nERROR.  Exiting catfixwcs.\n\n");
    return 1;
  }
}

