/*
 * catcompare.c
 *
 * Usage: catcompare [catfile1] [format1] [catfile2] [format2] 
 *
 * This program reads in two catalog files and finds matches based on the
 *  positions in each catalog.  The output is a file containing the offsets
 *  between each detected pair.
 * The input catalogs are probably produced by the run_sext_*.sh scripts.
 *
 * Revision history:
 *  2009May05, Chris Fassnacht (CDF) - A modified version of catcomb.c
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
 * Function declarations
 *
 */

void help_catcompare();

/*.......................................................................
 *
 * Main program
 *
 */

int main(int argc, char *argv[])
{
  int i,j,k;               /* Looping variables */
  int no_error=1;          /* Flag set to 0 on error */
  int calc_offsets=0;      /* Set to 1 to calculate offsets */
  int nfiles;              /* Number of input catalogs */
  int firstfile=1;         /* argv index of first file */
  int fileindex;           /* Index of current file being read */
  int *nlines=NULL;        /* Number of lines in the input catalog */
  int ninit=0;             /* Number of lines in initial catalog */
  int nmatch=0;            /* Number of additional lines when combinining */
  int ncent;               /* Number of lines in cent pos file (should be 1) */
  int *format=NULL;        /* Format of input/output files */
  int mindex;              /* Catalog index of closest match */
  int *iptr;               /* Pointer to navigate int arrays */
  int **id={NULL};         /* Output IDs */
  float *fptr;             /* Pointer for navigating float arrays */
  double tmpd;             /* Coordinate offset */
  double dmatch;           /* Cutoff distance for matching */
  char **infiles={NULL};   /* Input file names */
  char outfile[MAXC];      /* Filename for output file */
  char posfile[MAXC];      /* Filename for optional central position file */
  char line[MAXC];         /* General string for reading variables */
  Pos *dpos=NULL;          /* Array to hold offsets */
  Skypos currpos;          /* Central position for distance calculations */
  Skypos *skypos=NULL;     /* Array of sky positions */
  Skypos *skptr;           /* Pointer to navigate skypos */
  Secat mindist;           /* Catalog entry for minimum offset */
  Secat **dp={NULL};       /* Offsets between catalogs */
  Secat *dptr;             /* Pointer for navigating dp */
  Secat **incat={NULL};    /* Input catalogs */
  Secat *centpos=NULL;     /* Central position, if offsets are requested */
  Secat *sptr1,*sptr2;     /* Pointers to navigate catalogs */
  Secat *sptr3;            /* Pointers to navigate catalogs */
  FILE *ofp=NULL;          /* Output file pointer */

  /*
   * Check the command line invocation
   */

  if(argc < 5) {
    help_catcompare();
    return 1;
  }
  printf("\n");

  /*
   * Need a check on even number of passed arguments, too
   */

  if(((argc-1)/2.0 - floor((argc-1)/2))>0){
    fprintf(stderr,"\nError. Odd number of arguments!\n\n");
    help_catcompare();
    return 1;
  }
  else {
    nfiles = (argc - 1)/2;
  }

  /*
   * Allocate first level for input catalogs and output arrays
   */

  infiles = (char **) malloc(sizeof(char *) * nfiles);
  if(!infiles) {
    fprintf(stderr,"ERROR:  Insufficient memory for input filename array.\n");
    return 1;
  }
  incat = (Secat **) malloc(sizeof(Secat *) * nfiles);
  if(!incat) {
    fprintf(stderr,"ERROR:  Insufficient memory for input catalog array.\n");
    return 1;
  }
  id = (int **) malloc(sizeof(int *) * nfiles);
  if(!id) {
    fprintf(stderr,"ERROR:  Insufficient memory for output id array.\n");
    return 1;
  }
  dp = (Secat **) malloc(sizeof(Secat *) * nfiles);
  if(!dp) {
    fprintf(stderr,"ERROR:  Insufficient memory for output dpos array.\n");
    return 1;
  }

  if(!(nlines = new_intarray(nfiles,1)))
    no_error = 0;

  if(!(format = new_intarray(nfiles,1)))
    no_error = 0;

  /*
   * Read in the catalogs
   */

  i = 0;
  iptr = format;
  while(no_error && i<nfiles) {
    if(!(infiles[i] = new_string(80)))
      no_error = 0;
    else {
      fileindex = 2*i + firstfile;
      strcpy(infiles[i],argv[fileindex]);
      if((sscanf(argv[fileindex+1],"%d",iptr)) != 1) {
	fprintf(stderr,"ERROR. Expected an integer format after name of ");
	fprintf(stderr,"catalog number %d\n",i+1);
	no_error = 0;
      }
      if(no_error) {
	switch(*iptr) {
	case 11: 
	  *iptr = 1;
	  break;
	case 12:
	  *iptr = 10;
	  break;
	default:
	  break;
	}
	if(!(incat[i] = read_secat(infiles[i],'#',&nlines[i],*iptr)))
	  no_error = 0;
	else
	  ninit += nlines[i];
      }
    }
    i++;
    iptr++;
  }

  /*
   * Allocate memory for output arrays
   */

  i=0;
  while(no_error && i<nfiles) {
    if(!(id[i] = new_intarray(ninit,1)))
      no_error = 0;
    if(!(dp[i] = new_secat(ninit)))
      no_error = 0;
    i++;
  }

  /*
   * Allocate memory for the skypos array
   */

  if(!(skypos = new_skypos(1)))
    no_error = 0;
  else 
    skptr = skypos;

  /*
   * Loop over other catalogs, finding matches in each one
   */

  i=1;
  while(no_error && i<nfiles) {

    nmatch = 0;
    iptr = id[i];

    /*
     * Open output file
     */

    sprintf(outfile,"matchcat1%d.txt",i+1);
    if(!(ofp=open_writefile(outfile)))
      no_error = 0;
    else {
      fprintf(ofp,"#  1 CATID_01    ID in first catalog\n");
      fprintf(ofp,"#  2 DRA_ASEC    RA offset in arcsec\n");
      fprintf(ofp,"#  3 DDEC_ASEC   Dec offset in arcsec\n");
      fprintf(ofp,"#  4 DTOT_ASEC   Total offset in arcsec\n");
      fprintf(ofp,"#  5 CATID_%02d    ID in catalog %d\n",i+1,i+1);
	      
    }

    /*
     * Get maximum distance for matching
     */

    printf("\nEnter maximum offset (in arcsec) for a match ");
    printf("between catalog %d and\n the first catalog(s): ",i+1);
    fgets(line,MAXC,stdin);
    while(sscanf(line,"%lf",&dmatch) != 1) {
      fprintf(stderr," ERROR: bad value.  Try again: ");
      fgets(line,MAXC,stdin);
    }

    /*
     * Loop through the next input catalog
     */

    for(j=0,sptr2=incat[i]; j<nlines[i]; j++,sptr2++) {

      /*
       * Get the position of the current object
       */

      currpos = sptr2->skypos;

      /*
       * Calculate the offset between the current and the first
       *  entry in the first catalog.
       */

      *skptr = incat[0]->skypos;
      if(!(dpos = dspos2xy(currpos,skypos,1)))
	no_error = 0;
      else {
	mindex = 0;
	mindist = *incat[0];
	mindist.dx = dpos->x;
	mindist.dy = dpos->y;
	mindist.dpos = sqrt(dpos->x * dpos->x + dpos->y * dpos->y);
      }
      dpos = del_pos(dpos);

      /*
       * Loop through the first catalog, and find the closest match
       */

      if(no_error) {
	for(k=0,sptr1=incat[0]; k<nlines[0]; k++,sptr1++) {
	  *skptr = sptr1->skypos;
	  if(!(dpos = dspos2xy(currpos,skypos,1)))
	    no_error = 0;
	  else if((tmpd = sqrt(dpos->x * dpos->x + dpos->y * dpos->y)) <
		  mindist.dpos) {
	    mindist = *sptr1;
	    mindist.dx = dpos->x;
	    mindist.dy = dpos->y;
	    mindist.dpos = tmpd;
	    mindex = k;
	  }
	  dpos = del_pos(dpos);
	}
      }

      /*
       * Check to see if closest match is within cutoff.  If it is, add it
       *  to the list of matches and to the output file.
       */

      if(mindist.dpos < dmatch) {
	sptr1 = incat[0] + mindex;
	sptr1->id[0] = i * 10000 + sptr2->id;
	iptr = id[i] + mindex;
	*iptr = sptr2->id;
	dptr = dp[i] + mindex;
	dptr->dpos = mindist.dpos;
	nmatch++;
	fprintf(ofp,"%05d %8.2f %8.2f %8.2f %5d\n",sptr1->id,mindist.dx,
		mindist.dy,mindist.dpos,sptr2->id);
      }
    }
    printf("Between catalog %d and the first catalog, there are %d matches\n",
	   i+1,nmatch);
    if(ofp)
      fclose(ofp);
    i++;
  }

  /*
   * Clean up and exit
   */

  skypos = del_skypos(skypos);
  dpos = del_pos(dpos);
  for(i=0; i<nfiles; i++) {
    incat[i] = del_secat(incat[i]);
    id[i] = del_intarray(id[i]);
    dp[i] = del_secat(dp[i]);
    infiles[i] = del_string(infiles[i]);
  }
  if(incat)
    free(incat);
  if(infiles)
    free(infiles);
  nlines = del_intarray(nlines);
  format = del_intarray(format);
  if(no_error) {
    printf("\nProgram catcompare finished.\n\n");
    return 0;
  }
  else {
    fprintf(stderr,"\nERROR.  Exiting catcompare.\n\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function help_catcompare
 *
 * Prints useful information for running catcompare
 *
 * Inputs: (none)
 *
 * Output: (none)
 */

void help_catcompare()
{
  char line[MAXC];  /* General string */

  fprintf(stderr,"\nUsage: \n");
  fprintf(stderr,"  catcompare [catfile1] [format1] ... [catfileN] ");
  fprintf(stderr,"[formatN]\n\n");
  fprintf(stderr," At least 2 input catalogs are required\n\n");
  fprintf(stderr,"The format flags indicate the formats of the input ");
  fprintf(stderr,"catalogs.\n");
  fprintf(stderr,"Hit return to see format options: ");
  fgets(line,MAXC,stdin);
  secat_format();

}
