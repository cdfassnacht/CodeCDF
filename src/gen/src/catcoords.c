/*
 * catcoords.c
 *
 * Usage: catcoords [catalog] [output_filename] ([format])
 *
 * This program extracts the (x,y) and (RA,Dec) coordinates for each
 *  object in a SExtractor catalog.
 *
 *  The input catalog is probably produced by the run_sext_astrom.sh script.
 *
 * 11Aug03 CDF,  Modification of catdistcalc.c
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
  char catfile[MAXC];   /* Filename for input catalog */
  char outfile[MAXC];   /* Filename for distcalc-like output file */
  Secat *secat=NULL;    /* Data array from catalog */
  Secat *sptr;          /* Pointer to navigate secat */
  FILE *ofp=NULL;       /* Output file pointer */

  /*
   * Check the command line invocation
   */

  if(argc < 3) {
    fprintf(stderr,"\nUsage: catcoords [catfile] [outfile] ");
    fprintf(stderr,"([format])\n\n");
    return 1;
  }
  printf("\n");

  /*
   * Get the inputs from the command line.
   */

  strcpy(catfile,argv[1]);
  strcpy(outfile,argv[2]);
  if(argc == 4) {
    if(sscanf(argv[3],"%d",&format) != 1)
      no_error = 0;
  }
  else
    format = 7;

  /*
   * Fill the data structure
   */

  if(no_error)
    if(!(secat = read_secat(catfile,'#',&ncat,format)))
      no_error = 0;

  /*
   * Write output file
   */

  if(no_error)
    if(!(ofp = open_writefile(outfile)))
      no_error = 0;

  /*
   * Load the skypos values into the skypos array
   */

  if(no_error)
    for(i=0,sptr=secat; i<ncat; i++,sptr++)
      fprintf(ofp,"%04d %7.2f %7.2f  %02d %02d %07.4f %+03d %02d %06.3f\n",
	      sptr->id,sptr->x,sptr->y,
	      sptr->skypos.hr,sptr->skypos.min,sptr->skypos.sec,
	      sptr->skypos.deg,sptr->skypos.amin,sptr->skypos.asec);

  /*
   * Clean up and exit
   */

  secat = del_secat(secat);
  if(ofp)
    fclose(ofp);

  if(no_error) {
    printf("\nProgram catcoords finished.\n\n");
    return 0;
  }
  else {
    fprintf(stderr,"\nERROR.  Exiting catcoords.\n\n");
    return 1;
  }
}

