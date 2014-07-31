/*
 * catsort.c
 *
 * Usage: catsort [position_filename] [input_catalog] [output_filename] 
 *   [calcmethod] ([format])
 *
 * This program sorts a SExtractor or SDSS catalog in terms of increasing 
 *  distance from some position, which is provided via the posfile parameter.
 * The calcmethod parameter sets whether offsets are calculated based
 *  on (x,y) coordinates (calcmethod=xy) or (RA,Dec) coordinates
 *  (calcmethod=radec). 
 *
 * Revision history:
 *  2003Jul29, Chris Fassnacht (CDF) - First working version
 *  2003Aug13, CDF - Added purge_cat call to eliminate catalog members with
 *                    SExtractor fit flag values > MASTERLIM
 *  2005Feb24, CDF - Changed approach for (RA,Dec) calculation.  Now use an
 *                    externally provided RA and Dec for the central position
 *                    instead of the RA and Dec associated with the (x,y)
 *                    position of the lens.
 *  2007Jan09, CDF - Added a parallel version for SDSS-format catalogs.
 *                   Improved documentation.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "structdef.h"
#include "dataio.h"
#include "catlib.h"

/*.......................................................................
 *
 * Function declarations
 *
 */

void catsort_help();
int sort_secat(char *catfile, char *posfile, char *outfile, char *calcmethod, 
	       int format);
int sort_sdss(char *catfile, char *posfile, char *outfile, int format);


/*.......................................................................
 *
 * Main program
 *
 */

int main(int argc, char *argv[])
{
  int no_error=1;          /* Flag set to 0 on error */
  int format;              /* Format of input/output files */
  char catfile[MAXC];      /* Filename for input catalog */
  char posfile[MAXC];      /* File containing central position */
  char outfile[MAXC];      /* Filename for input catalog */
  char calcmethod[MAXC];   /* Sets method of calculating offsets */

  /*
   * Check the command line invocation
   */

  if(argc < 4) {
    catsort_help();
    return 1;
  }
  printf("\n");

  /*
   * Get the inputs from the command line.
   */

  strcpy(posfile,argv[1]);
  strcpy(catfile,argv[2]);
  strcpy(outfile,argv[3]);
  strcpy(calcmethod,argv[4]);
  if(argc == 6) {
    if(sscanf(argv[5],"%d",&format) != 1)
      no_error = 0;
  }
  else
    format = 6;

  /*
   * Do the sorting
   */

  if(no_error) {
    switch(format){
    case 11: case 12:
      if(sort_sdss(catfile,posfile,outfile,format))
	no_error = 1;
      break;
    default:
      if(sort_secat(catfile,posfile,outfile,calcmethod,format))
	no_error = 0;
    }
  }

  /*
   * Exit
   */

  if(no_error) {
    printf("\nProgram catsort finished.\n\n");
    return 0;
  }
  else {
    fprintf(stderr,"\nERROR.  Exiting catsort.\n\n");
    return 1;
  }
}


/*.......................................................................
 *
 * Function sort_secat
 *
 * Here, the input file is a SExtractor output file, or at least a file that
 *  has a format compatible with a SExtractor format, in contrast to a
 *  SDSS-like file that has magnitudes for multiple bands all in one file.
 *  This function does the sorting and creates the output for a Secat 
 *  structure.
 *
 * Inputs: 
 *  char *catfile         filename for input catalog
 *  char *posfile         file containing central position
 *  char *outfile         output filename
 *  char *calcmethod      method of calculating offsets
 *  int format            format of input catalog file
 */

int sort_secat(char *catfile, char *posfile, char *outfile, char *calcmethod, 
	       int format)
{
  int i;                   /* Looping variable */
  int no_error=1;          /* Flag set to 0 on error */
  int ninit;               /* Number of lines in the input catalog */
  int ncat;                /* Number of lines in the final catalog */
  int ncent;               /* Number of lines in posfile (should be 1) */
  int centindex;           /* Catalog index of closest match to central pos */
  Pos *offsets=NULL;       /* Array to hold offsets */
  Pos *pptr;               /* Pointer to navigate offsets */
  Skypos *skypos=NULL;     /* Array of sky positions */
  Skypos *skptr;           /* Pointer to navigate skypos */
  Secat *centpos=NULL;     /* Central position for distance calculations */
  Secat *initcat=NULL;     /* Data array from input catalog */
  Secat *secat=NULL;       /* Data array from catalog, after purging */
  Secat *sptr;             /* Pointer to navigate secat */

  /*
   * Fill the data structure
   */

  if(no_error)
    if(!(initcat = read_secat(catfile,'#',&ninit,format)))
      no_error = 0;

  /*
   * Purge input catalog of entries with SExtractor flags greater than 
   *  MASTERLIM to create the final catalog.
   */

  if(no_error)
    if(!(secat = purge_cat(initcat,ninit,catfile,&ncat,MASTERLIM)))
      no_error = 0;

  /*
   * Calculate (RA,Dec) offsets if calcmethod == radec
   */

  if(strcmp(calcmethod,"radec") == 0 && no_error) {

    /*
     * Read the central position
     */

    if(no_error)
      if(!(centpos = read_distcalc(posfile,'#',&ncent,0)))
	no_error = 0;

    /*
     * Allocate memory for the skypos array
     */

    if(!(skypos = new_skypos(ncat)))
      no_error = 0;

    /*
     * Load the skypos values into the skypos array
     */

    else {
      for(i=0,sptr=secat,skptr=skypos; i<ncat; i++,sptr++,skptr++) {
	*skptr = sptr->skypos;
	sprintf(skptr->label,"%04d",sptr->id);
      }

    /*
     * Calculate the offsets
     */

      if(!(offsets = dspos2xy(centpos->skypos,skypos,ncat)))
	no_error = 0;
    }

    /*
     * Put the offsets into the catalog structure
     */

    if(no_error) {
      for(i=0,sptr=secat,pptr=offsets; i<ncat; i++,sptr++,pptr++) {
	sptr->dx = pptr->x;
	sptr->dy = pptr->y;
	sptr->dpos = sqrt(pptr->x * pptr->x + pptr->y * pptr->y);
      }
    }
  }

  /*
   * Otherwise, find lens system in catalog and calculate 
   *  (x,y) distances from it.
   */

  else {
    if(no_error) {
      if(find_lens(secat,ncat,&centindex))
	no_error = 0;
      else {
	printf("centindex = %d, pos = %7.2f %7.2f\n",centindex,
	       (secat+centindex)->x,(secat+centindex)->y);
      }
    }
  }

  /*
   * Sort the catalog in order of
   *  increasing distance from the central position.
   */

  if(no_error) {
    printf("\nSorting the catalog in order of increasing distance ");
    printf("from the central position...");
    qsort(secat,ncat,sizeof(secat[0]),dposcmp);
    printf(" Done.\n");
  }

  /*
   * Renumber the sources
   */

  if(no_error) {
    for(i=0,sptr=secat; i<ncat; i++,sptr++)
      sptr->id = i + 1;
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

  initcat = del_secat(initcat);
  secat = del_secat(secat);
  centpos = del_secat(centpos);
  skypos = del_skypos(skypos);
  offsets = del_pos(offsets);

  if(no_error)
    return 0;
  else {
    fprintf(stderr,"ERROR: sort_secat\n");
    return 0;
  }
}

/*.......................................................................
 *
 * Function sort_sdss
 *
 * Here, the input file is a SExtractor output file, or at least a file that
 *  has a format compatible with a SExtractor format, in contrast to a
 *  SDSS-like file that has magnitudes for multiple bands all in one file.
 *  This function does the sorting and creates the output for a Secat 
 *  structure.
 *
 * Inputs: 
 *  char *catfile         filename for input catalog
 *  char *posfile         file containing central position
 *  char *outfile         output filename
 *  char *calcmethod      method of calculating offsets
 *  int format            format of input catalog file
 */

int sort_sdss(char *catfile, char *posfile, char *outfile, int format)
{
  int i;                   /* Looping variable */
  int no_error=1;          /* Flag set to 0 on error */
  int ncat;                /* Number of lines in the final catalog */
  int ncent;               /* Number of lines in posfile (should be 1) */
  Pos *offsets=NULL;       /* Array to hold offsets */
  Pos *pptr;               /* Pointer to navigate offsets */
  Skypos *skypos=NULL;     /* Array of sky positions */
  Skypos *skptr;           /* Pointer to navigate skypos */
  Secat *centpos=NULL;     /* Central position for distance calculations */
  SDSScat *scat=NULL;      /* Data array from catalog, after purging */
  SDSScat *sptr;           /* Pointer to navigate secat */

  /*
   * Fill the data structure
   */

  if(no_error)
    if(!(scat = read_sdss(catfile,'#',&ncat,format)))
      no_error = 0;

  /*
   * Read the central position
   */

  if(no_error)
    if(!(centpos = read_distcalc(posfile,'#',&ncent,0)))
      no_error = 0;

  /*
   * Allocate memory for the skypos array
   */

  if(!(skypos = new_skypos(ncat)))
    no_error = 0;

  /*
   * Load the skypos values into the skypos array
   */

  else {
    for(i=0,sptr=scat,skptr=skypos; i<ncat; i++,sptr++,skptr++) {
      *skptr = sptr->skypos;
      sprintf(skptr->label,"%04d",sptr->id);
    }

    /*
     * Calculate the offsets
     */

    if(!(offsets = dspos2xy(centpos->skypos,skypos,ncat)))
      no_error = 0;
  }

  /*
   * Put the offsets into the catalog structure
   */

  if(no_error) {
    for(i=0,sptr=scat,pptr=offsets; i<ncat; i++,sptr++,pptr++) {
      sptr->dx = pptr->x;
      sptr->dy = pptr->y;
      sptr->dpos = sqrt(pptr->x * pptr->x + pptr->y * pptr->y);
    }
  }

  /*
   * Sort the catalog in order of
   *  increasing distance from the central position.
   */

  if(no_error) {
    printf("\nSorting the catalog in order of increasing distance ");
    printf("from the central position...");
    qsort(scat,ncat,sizeof(scat[0]),dposcmp_sdss);
    printf(" Done.\n");
  }

  /*
   * Renumber the sources
   */

  if(no_error) {
    for(i=0,sptr=scat; i<ncat; i++,sptr++)
      sptr->id = i + 1;
  }

  /*
   * Write output file
   */

  if(no_error)
     if(write_sdss(scat,ncat,outfile,format))
       no_error = 0;

  /*
   * Clean up and exit
   */

  scat = del_sdsscat(scat);
  centpos = del_secat(centpos);
  skypos = del_skypos(skypos);
  offsets = del_pos(offsets);

  if(no_error)
    return 0;
  else {
    fprintf(stderr,"ERROR: sort_secat\n");
    return 0;
  }
}

/*.......................................................................
 *
 * Function catsort_help
 *
 * Prints useful information for running catsort.
 *
 * Inputs: (none)
 *
 * Output: (none)
 */

void catsort_help()
{
  char line[MAXC];  /* General string */

  fprintf(stderr,"\nUsage: catsort [posfile] [catfile] [outfile] ");
  fprintf(stderr,"[calcmethod] ([format])\n\n");
  fprintf(stderr,"  catfile is the file containing the catalog\n");
  fprintf(stderr,"  posfile is the file containing the RA, Dec of the");
  fprintf(stderr," central object.\n");
  fprintf(stderr,"   Format: id RA_hr RA_min RA_sec Dec_deg Dec_min");
  fprintf(stderr," Dec_sec\n");
  fprintf(stderr,"  calcmethod sets the method for calculating the offsets.");
  fprintf(stderr,"\n");
  fprintf(stderr,"    For calcmethod = xy distances are based on");
  fprintf(stderr," (x,y) coordinates.\n");
  fprintf(stderr,"    For calcmethod = radec distances are based on");
  fprintf(stderr," (RA,Dec) coordinates.\n\n");
  fprintf(stderr," The optional format flag indicates the format of the");
  fprintf(stderr," input file:\n");
  fprintf(stderr,"Hit return to see the format options: ");
  fgets(line,MAXC,stdin);
  secat_format();
}
