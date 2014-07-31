/*
 * catdistcalc.c
 *
 * Usage: catdistcalc [position_filename] [catalog] [output_filename] ([format])
 *
 * This program calculates offsets between a given input position
 *  and all the positions of objects in a SExtractor catalog, based 
 *  on the skypos values for each object.
 *
 *  The input catalog is probably produced by the run_sext_astrom.sh script.
 *
 * The optional format flag describes the format of the input file
 *     0 ==> ID, x, y
 *     1 ==> Format produced by run_sext.sh in
 *        Lenses/SExtractor directory
 *     2 ==> Temporary run_sext.sh format during
 *        format change
 *     3 ==> Format produced by matchcat.c
 *     4 ==> Format returned by 2MASS server (among
 *        others?). alpha, delta, ...many other cols...
 *        Just take alpha and delta
 *     6 ==> Another hardwired SExtractor format
 *     7 ==> Another hardwired SExtractor format
 *     8 ==> Format from USNO star web page
 *        hr, min, sec, deg, amin, asec, mag
 *
 * 29Jul03 CDF
 * v01Aug03 CDF, Changed central position from first object in the
 *                SExtractor catalog to an arbitrary position, given by
 *                the user in the form of another input file.
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
  int ncent;            /* Number of lines in posfile (should be 1) */
  int ncat;             /* Number of lines in the input catalog */
  int format;           /* Format of input/output files */
  char posfile[MAXC];   /* Filename for input position file */
  char catfile[MAXC];   /* Filename for input catalog */
  char outfile[MAXC];   /* Filename for distcalc-like output file */
  Pos *offsets=NULL;    /* Array to hold offsets */
  Skypos *skypos=NULL;  /* Array of sky positions */
  Skypos *skptr;        /* Pointer to navigate skypos */
  Secat *centpos=NULL;  /* Central position for distance calculations */
  Secat *secat=NULL;    /* Data array from catalog */
  Secat *sptr;          /* Pointer to navigate secat */
  FILE *ofp=NULL;       /* Output file pointer */

  /*
   * Check the command line invocation
   */

  if(argc < 4) {
    fprintf(stderr,"\nUsage: catdistcalc [posfile] [catfile] [outfile] ");
    fprintf(stderr,"([format])\n\n");
    fprintf(stderr," The optional format flag indicates the format of the");
    fprintf(stderr," input file:\n");
    fprintf(stderr,"   0 ==> ID, x, y\n");
    fprintf(stderr,"   1 ==> Format produced by run_sext.sh in\n");
    fprintf(stderr,"      Lenses/SExtractor directory\n");
    fprintf(stderr,"   2 ==> Temporary run_sext.sh format during\n");
    fprintf(stderr,"      format change\n");
    fprintf(stderr,"   3 ==> Format produced by matchcat.c\n");
    fprintf(stderr,"   4 ==> Format returned by 2MASS server, etc.\n");
    fprintf(stderr,"      alpha, delta, ...many other cols...\n");
    fprintf(stderr,"      Just take alpha and delta\n");
    fprintf(stderr,"   6 ==> Another hardwired SExtractor format\n");
    fprintf(stderr,"   7 ==> Another hardwired SExtractor format\n");
    fprintf(stderr,"   8 ==> Format returned by USNO star web page\n");
    fprintf(stderr,"      hr, min, sec, deg, amin, asec, mag\n\n");
    return 1;
  }
  printf("\n");

  /*
   * Get the inputs from the command line.
   */

  strcpy(posfile,argv[1]);
  strcpy(catfile,argv[2]);
  strcpy(outfile,argv[3]);
  if(argc == 5) {
    if(sscanf(argv[4],"%d",&format) != 1)
      no_error = 0;
  }
  else
    format = 7;

  /*
   * Read the input position
   */

  if(no_error)
    if(!(centpos = read_distcalc(posfile,'#',&ncent,0)))
      no_error = 0;

  /*
   * Fill the data structure
   */

  if(no_error)
    if(!(secat = read_secat(catfile,'#',&ncat,format)))
      no_error = 0;

  /*
   * Allocate memory for the skypos array
   */

  if(no_error)
    if(!(skypos = new_skypos(ncat)))
      no_error = 0;

  /*
   * Load the skypos values into the skypos array
   */

  if(no_error)
    for(i=0,sptr=secat,skptr=skypos; i<ncat; i++,sptr++,skptr++) {
      *skptr = sptr->skypos;
      if(format == 9)
	sprintf(skptr->label,"%s",sptr->name);
      else
	sprintf(skptr->label,"%04d",sptr->id);
    }

  /*
   * Calculate the offsets
   */

  if(no_error)
    if(!(offsets = dspos2xy(centpos->skypos,skypos,ncat)))
      no_error = 0;

  /*
   * Write output file
   */

  if(no_error)
    if(!(ofp = open_writefile(outfile)))
      no_error = 0;

  if(no_error)
    if(print_offsets(centpos->skypos,skypos,offsets,ncat,ofp))
       no_error = 0;

  /*
   * Clean up and exit
   */

  secat = del_secat(secat);
  centpos = del_secat(centpos);
  skypos = del_skypos(skypos);
  offsets = del_pos(offsets);
  if(ofp)
    fclose(ofp);

  if(no_error) {
    printf("\nProgram catdistcalc finished.\n\n");
    return 0;
  }
  else {
    fprintf(stderr,"\nERROR.  Exiting catdistcalc.\n\n");
    return 1;
  }
}

