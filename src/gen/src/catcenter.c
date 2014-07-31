/*
 * catcenter.c
 *
 * Usage: catcenter [input_catalog] [output_filename] 
 *   [calcmethod] ([format])
 *
 * This program finds the median position in a catalog.
 * The calcmethod parameter sets whether offsets are calculated based
 *  on (x,y) coordinates (calcmethod=xy) or (RA,Dec) coordinates
 *  (calcmethod=radec). 
 *
 * Revision history:
 *  2007Feb09, Chris Fassnacht (CDF) - First working version
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "structdef.h"
#include "dataio.h"
#include "coords.h"
#include "catlib.h"

/*.......................................................................
 *
 * Function declarations
 *
 */

void catcenter_help();
int med_secat(char *catfile, char *outfile, char *calcmethod, int format);
int med_sdss(char *catfile, char *outfile, int format);


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

  if(argc < 3) {
    catcenter_help();
    return 1;
  }
  printf("\n");

  /*
   * Get the inputs from the command line.
   */

  strcpy(catfile,argv[1]);
  strcpy(outfile,argv[2]);
  strcpy(calcmethod,argv[3]);
  if(argc == 5) {
    if(sscanf(argv[4],"%d",&format) != 1)
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
      if(med_sdss(catfile,outfile,format))
	no_error = 1;
      break;
    default:
      if(med_secat(catfile,outfile,calcmethod,format))
	no_error = 0;
    }
  }

  /*
   * Exit
   */

  if(no_error) {
    printf("\nProgram catcenter finished.\n\n");
    return 0;
  }
  else {
    fprintf(stderr,"\nERROR.  Exiting catcenter.\n\n");
    return 1;
  }
}


/*.......................................................................
 *
 * Function med_secat
 *
 * Here, the input file is a SExtractor output file, or at least a file that
 *  has a format compatible with a SExtractor format, in contrast to a
 *  SDSS-like file that has magnitudes for multiple bands all in one file.
 *  This function does the sorting and finds the median position
 *
 * Inputs: 
 *  char *catfile         filename for input catalog
 *  char *outfile         output filename
 *  char *calcmethod      method of calculating offsets
 *  int format            format of input catalog file
 */

int med_secat(char *catfile, char *outfile, char *calcmethod, int format)
{
  int i;                   /* Looping variable */
  int no_error=1;          /* Flag set to 0 on error */
  int ncat;                /* Number of lines in the input catalog */
  int centindex;           /* Catalog index of closest match to central pos */
  double *ra=NULL;         /* Vector for RA (or x) coordinate */
  double *dec=NULL;        /* Vector for Dec (or y) coordinate */
  double *dptr;            /* Pointer to navigate ra and dec */
  double *tptr;            /* Pointer to navigate ra and dec */
  Pos centpos;             /* Output central position */
  Skypos centskypos;       /* Central position in hms format */
  Secat *secat=NULL;       /* Data array from catalog */
  Secat *sptr;             /* Pointer to navigate secat */

  /*
   * Fill the data structure
   */

  if(no_error)
    if(!(secat = read_secat(catfile,'#',&ncat,format)))
      no_error = 0;

  /*
   * Allocate memory for RA and Dec (or x and y) arrays
   */

  if(no_error)
    if(!(ra = new_doubarray(ncat)))
      no_error = 0;
  if(no_error)
    if(!(dec = new_doubarray(ncat)))
      no_error = 0;

  /*
   * Sort (RA,Dec) coordinates if calcmethod == radec
   */

  if(strcmp(calcmethod,"radec") == 0 && no_error) {

    /*
     * Fill the RA and Dec arrays
     */

    for(i=0,sptr=secat,dptr=ra; i<ncat; i++,sptr++,dptr++)
      *dptr = sptr->alpha;
    for(i=0,sptr=secat,dptr=dec; i<ncat; i++,sptr++,dptr++)
      *dptr = sptr->delta;
  }

  /*
   * Otherwise, fill the ra and dec arrays with the x and y positions from
   *  the input catalog.
   */

  else {
    if(no_error) {

      /*
       * Fill the RA and Dec arrays
       */

      for(i=0,sptr=secat,dptr=ra; i<ncat; i++,sptr++,dptr++)
	*dptr = sptr->x;
      for(i=0,sptr=secat,dptr=dec; i<ncat; i++,sptr++,dptr++)
	*dptr = sptr->y;

    }
  }

  /*
   * Sort the RA and Dec arrays
   */

  qsort(ra,ncat,sizeof(ra[0]),dcmp);
  qsort(dec,ncat,sizeof(dec[0]),dcmp);

  /*
   * Find median position by taking middle values of the ra and dec arrays.  
   * The definition of "middle value" depends on whether the number of
   *   elements in the array (ncat) is even or odd.
   */

  if((ncat/2.0) - (int) (ncat/2.0) == 0) {
    centindex = ncat / 2;
    centpos.x = (ra[centindex-1] + ra[centindex])/2.0;
    centpos.y = (dec[centindex-1] + dec[centindex])/2.0;
  }
  else {
    centindex = (ncat - 1) / 2;
    centpos.x = ra[centindex];
    centpos.y = dec[centindex];
  }

  /*
   * Print out the central position.  If calcmethod is radec, also express
   *  the answer in hms format.
   */

  if(no_error) {
    printf("\n");
    for(i=0,dptr=ra,tptr=dec; i<ncat; i++,dptr++,tptr++) {
      deg2spos(*dptr,*tptr,&centskypos);
      print_skypos(stdout,centskypos);
      printf("\n");
    }
    printf("\nMedian position in the catalog is %lf %lf\n",centpos.x,centpos.y);
    if(strcmp(calcmethod,"radec") == 0 && no_error) {
      deg2spos(centpos.x,centpos.y,&centskypos);
      print_skypos(stdout,centskypos);
    }
  }


  /*
   * Clean up and exit
   */

  secat = del_secat(secat);
  ra = del_doubarray(ra);
  dec = del_doubarray(dec);

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
 *  This function does the sorting and finds the median position.
 *
 * Inputs: 
 *  char *catfile         filename for input catalog
 *  char *outfile         output filename
 *  char *calcmethod      method of calculating offsets
 *  int format            format of input catalog file
 */

int med_sdss(char *catfile, char *outfile, int format)
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
#if 0
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
#endif
}

/*.......................................................................
 *
 * Function catcenter_help
 *
 * Prints useful information for running catcenter.
 *
 * Inputs: (none)
 *
 * Output: (none)
 */

void catcenter_help()
{
  fprintf(stderr,"\nUsage: catcenter [catfile] [outfile] ");
  fprintf(stderr,"[calcmethod] ([format])\n\n");
  fprintf(stderr,"  catfile is the file containing the catalog\n");
  fprintf(stderr,"  calcmethod sets the method for calculating the offsets.");
  fprintf(stderr,"\n");
  fprintf(stderr,"    For calcmethod = xy distances are based on");
  fprintf(stderr," (x,y) coordinates.\n");
  fprintf(stderr,"    For calcmethod = radec distances are based on");
  fprintf(stderr," (RA,Dec) coordinates.\n\n");
  fprintf(stderr," The optional format flag indicates the format of the");
  fprintf(stderr," input file:\n");
  secat_format();
}
