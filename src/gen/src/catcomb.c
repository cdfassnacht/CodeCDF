/*
 * catcomb.c
 *
 * Usage: catcomb [catfile1] [format1] [catfile2] [format2] ....
 *
 * This program combines multiple catalogs of the same field, based
 *  on the positions of the catalog objects.
 * The input catalogs are probably produced by the run_sext_*.sh scripts.
 *
 * Revision history:
 *  2003Jul29, Chris Fassnacht (CDF) - First working version
 *  2003Aug13 CDF - Added purge_cat call to eliminate catalog members with
 *                   SExtractor fit flag values > MASTERLIM
 *  2007Jan09 CDF - Made more general by adding format flags to the command
 *                   line.
 *  2007Jul30 CDF - Much more informative output format
 *  2007Aug02 CDF - Added check for input SDSS format
 *  2007Jul04 CDF - Added option to calculate offsets from a fiducial
 *                   central position.  
 *                  Put help info into the new help_catcomb function.
 *  2009Jan24 CDF - Modified so that format 10 now includes magnitude info
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

void help_catcomb();
int write_matchcat(Secat *matchcat, int ncat, int **id, double **dp,
		   float **outmag, char **infiles, int *format, int nfiles, 
		   char *outfile);


/*.......................................................................
 *
 * Main program
 *
 */

int main(int argc, char *argv[])
{
  int i,j,k,m;             /* Looping variables */
  int no_error=1;          /* Flag set to 0 on error */
  int calc_offsets=0;      /* Set to 1 to calculate offsets */
  int nfiles;              /* Number of input catalogs */
  int firstfile=1;         /* argv index of first file */
  int fileindex;           /* Index of current file being read */
  int *nlines=NULL;        /* Number of lines in the input catalog */
  int ninit=0;             /* Number of lines in initial catalog */
  int nadd;                /* Number of additional lines when combining */
  int ncat;                /* Number of lines in the final catalog */
  int ncent;               /* Number of lines in cent pos file (should be 1) */
  int *format=NULL;        /* Format of input/output files */
  int mindex[10];          /* Catalog index of closest match */
  int **nmatch={NULL};      /* Number of matches within dmatch */
  int nmatch0;             /* Number of matches within dmatch - single source */
  int nmatchmax=0;         /* Maximum value of nmatch for the full comparison */
  int *iptr;               /* Pointer to navigate int arrays */
  int **id={NULL};         /* Output IDs */
  float **outmag={NULL};   /* Output magnitudes */
  float *fptr;             /* Pointer for navigating float arrays */
  double **dp={NULL};      /* Offsets between catalogs */
  double *dptr;            /* Pointer for navigating dp */
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
  Secat mindist[10];       /* Catalog entry for minimum offsets */
  Secat **incat={NULL};    /* Input catalogs */
  Secat *mastercat=NULL;   /* Data array from catalog, after purging */
  Secat *centpos=NULL;     /* Central position, if offsets are requested */
  Secat *sptr1,*sptr2;     /* Pointers to navigate catalogs */
  Secat *sptr3;            /* Pointers to navigate catalogs */

  /*
   * Check the command line invocation
   */

  if(argc < 5) {
    help_catcomb();
    return 1;
  }
  printf("\n");

  /*
   * Need a check on even number of passed arguments, too
   */

  if(((argc-1)/2.0 - floor((argc-1)/2))>0){
    fprintf(stderr,"\nError. Odd number of arguments!\n\n");
    help_catcomb();
    return 1;
  }

  /*
   * Check whether offset calculation is desired
   */

  if(strcmp(argv[1],"-c") == 0) {
    printf("Found -c flag.  Will calculate offsets\n");
    calc_offsets = 1;
    strcpy(posfile,argv[2]);
    nfiles = (argc - 3)/2;
    firstfile = 3;
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
  nmatch = (int **) malloc(sizeof(int *) * nfiles);
  if(!nmatch) {
    fprintf(stderr,"ERROR:  Insufficient memory for output id array.\n");
    return 1;
  }
  id = (int **) malloc(sizeof(int *) * nfiles);
  if(!id) {
    fprintf(stderr,"ERROR:  Insufficient memory for output id array.\n");
    return 1;
  }
  dp = (double **) malloc(sizeof(double *) * nfiles);
  if(!dp) {
    fprintf(stderr,"ERROR:  Insufficient memory for output dpos array.\n");
    return 1;
  }
  outmag = (float **) malloc(sizeof(float *) * nfiles);
  if(!outmag) {
    fprintf(stderr,"ERROR:  Insufficient memory for output outmag array.\n");
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
   * Allocate memory for master array and output arrays
   */

  if(no_error) {
    printf("\n%d total lines in all %d catalogs\n",ninit,nfiles);
    if(!(mastercat = new_secat(ninit))) {
      fprintf(stderr,"ERROR: Insufficient memory for combined catalog.\n");
      no_error = 0;
    }
  }

  i=0;
  while(no_error && i<nfiles) {
    if(!(id[i] = new_intarray(ninit,1)))
      no_error = 0;
    if(!(dp[i] = new_doubarray(ninit)))
      no_error = 0;
    if(!(outmag[i] = new_array(ninit,1)))
      no_error = 0;
    if(format[i] == 15 || format[i] == 16)
      for(j=0,fptr=outmag[i]; j<ninit; j++,fptr++)
	*fptr = -1.0;
    i++;
  }

  /*
   * Load master array with first catalog
   */

  if(no_error) {
    ncat = nlines[0];
    for(i=0,sptr1=mastercat,sptr2=incat[0],iptr=id[0],fptr=outmag[0]; 
	i<nlines[0]; i++,sptr1++,sptr2++,iptr++,fptr++) {
      *sptr1 = *sptr2;
      *iptr = sptr1->id;
      switch(format[0]) {
      case 6: case 7: case 9: case 10:
	*fptr = sptr1->mtot;
	break;
      case 15: case 16:
	*fptr = sptr1->zspec;
	break;
      default:
	*fptr = 0.0;
      }
    }
    printf("\nAfter catalog 1, there are %d master catalog members\n",ncat);
  }

  /*
   * Allocate memory for the skypos array
   */

  if(!(skypos = new_skypos(1)))
    no_error = 0;
  else
    skptr = skypos;

  /*
   * Loop over other catalogs, calculating distances for each object
   */

  i=1;
  while(no_error && i<nfiles) {

    nadd = 0;
    sptr3 = mastercat+ncat;
    iptr = id[i];

    /*
     * Get maximum distance for matching
     */

    printf("\nEnter maximum offset (in arcsec) for a match ");
    printf("between catalog %d and\n previous catalog(s): ",i+1);
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
       * Loop through the master catalog, and find all matches within
       * the specified match radius
       */

      if(no_error) {
	nmatch0 = 0;
	for(k=0,sptr1=mastercat; k<ncat; k++,sptr1++) {
	  *skptr = sptr1->skypos;
	  if(!(dpos = dspos2xy(currpos,skypos,1)))
	    no_error = 0;
	  if((tmpd = sqrt(dpos->x * dpos->x + dpos->y * dpos->y)) <
		  dmatch) {
	    mindist[nmatch0] = *sptr1;
	    mindist[nmatch0].dpos = tmpd;
	    mindist[nmatch0].matchid[0] = k;
	    switch(format[i]) {
	    case 6: case 7: case 9: case 10:
	      mindist[nmatch0].mtot = sptr2->mtot;
	      break;
	    case 15: case 16:
	      mindist[nmatch0].mtot = sptr2->zspec;
	      break;
	    }
	    nmatch0++;
	  }
	  dpos = del_pos(dpos);
	}
      }
      if(nmatch0>nmatchmax) {
	nmatchmax = nmatch0;
      }

      /*
       * Check to see if there are matches within the match radius.  
       * If not, add the current source to the end of the master catalog.
       */

      if(nmatch0>0) {
	/* Find the closest master catalog match */
	if(nmatch0>1) {
	  qsort(mindist,nmatch0,sizeof(mindist[0]),dposcmp);
	}

	/* Link the matched objects with the master catalog objext */
	for(k=0; k<nmatch0; k++) {
	  sptr1 = mastercat + mindist[k].matchid[0];
	  m = sptr1->nmatch;
	  sptr1->matchid[m] = sptr2->id;
	  sptr1->matchcat[m] = i;
	  sptr1->sep[m] = mindist[k].dpos;
	  sptr1->nmatch++;
	}
	iptr = id[i] + mindist[0].matchid[0];
	*iptr = sptr2->id;
	dptr = dp[i] + mindist[0].matchid[0];
	*dptr = mindist[0].dpos;
	fptr = outmag[i] + mindist[0].matchid[0];
	switch(format[i]) {
	case 6: case 7: case 9: case 10: case 15: case 16:
	  *fptr = mindist[0].mtot;
	  break;
	default:
	  *fptr = 0.0;
	}
      }

      /* Adding the current source to the master catalog */
      else {
	*sptr3 = *sptr2;
	sptr3->id = i * 10000 + sptr2->id;
	sptr3->matchid[0] = 0;
#if 0
	sptr1->dx = 0.0;
	sptr1->dy = 0.0;
	sptr1->dpos = 0.0;
#endif
	iptr = id[i] + ncat + nadd;
	*iptr = sptr2->id;
	fptr = outmag[i] + ncat + nadd;
	switch(format[i]) {
	case 6: case 7: case 9: case 10:
	  *fptr = sptr2->mtot;
	  break;
	case 15: case 16:
	  *fptr = sptr2->zspec;
	  break;
	default:
	  *fptr = 0.0;
	}
	nadd++;
	sptr3++;
      }
    }
    ncat += nadd;
    printf("After catalog %d, there are %d master catalog members\n",
	   i+1,ncat);
    i++;
    printf("%d\n",nmatchmax);
  }

  /*
   * If requested, (re)compute offsets between each object and the given
   *  central position.
   */

  if(calc_offsets && no_error) {
    printf("\nCalculating offsets\n");
    printf("--------------------------------------\n");
    if(!(centpos = read_distcalc(posfile,'#',&ncent,0))) {
      no_error = 0;
    }
    else {
      if(secat2offset(mastercat,ncat,DEGS,centpos[0],HMS) == ERROR)
	no_error = 0;
    }
  }

  /*
   * Write output file
   */

  if(no_error) {
    sprintf(outfile,"catcomb.out");
    printf("\n");
    if(write_matchcat(mastercat,ncat,id,dp,outmag,infiles,format,nfiles,outfile))
      no_error = 0;
  }

  /*
   * Clean up and exit
   */

  mastercat = del_secat(mastercat);
  skypos = del_skypos(skypos);
  dpos = del_pos(dpos);
  for(i=0; i<nfiles; i++) {
    incat[i] = del_secat(incat[i]);
    id[i] = del_intarray(id[i]);
    dp[i] = del_doubarray(dp[i]);
    outmag[i] = del_array(outmag[i]);
    infiles[i] = del_string(infiles[i]);
  }
  if(incat)
    free(incat);
  if(infiles)
    free(infiles);
  nlines = del_intarray(nlines);
  format = del_intarray(format);
  if(no_error) {
    printf("\nProgram catcomb finished.\n\n");
    return 0;
  }
  else {
    fprintf(stderr,"\nERROR.  Exiting catcomb.\n\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function write_matchcat
 *
 * Writes out combined catalog
 *
 * Inputs:
 *  Secat *matchcat            catalog of matches
 *  int ncat                   number of members of matchcat
 *  int **id                   id numbers of matched sources
 *  double **dp                offsets of matched sources
 *  float **outmag             magnitudes of matched sources
 *  char **infiles             names of the input files
 *  int *format                formats of the input files
 *  int nfiles                 number of input catalogs
 *  char *outfile              name of output file
 *
 * Output:
 *  int (0 or 1)
 *
 */

int write_matchcat(Secat *matchcat, int ncat, int **id, double **dp,
		   float **outmag, char **infiles, int *format, int nfiles, 
		   char *outfile)
{
  int i,j;         /* Looping variables */
  Secat *mptr;     /* Pointer to navicage matchcat */
  FILE *ofp=NULL;  /* Pointer to output file */

  if(!(ofp = open_writefile(outfile))) {
    fprintf(stderr,"ERROR: write_matchcat\n");
    return 1;
  }

  /*
   * Write output file header
   */

  fprintf(ofp,"# Result of running catcomb.c on the following data files:\n");
  for(i=0; i<nfiles; i++)
    fprintf(ofp,"#  %2d. %s\n",i+1,infiles[i]);
  fprintf(ofp,"#\n");
  fprintf(ofp,"# Column explanations:\n");
  fprintf(ofp,"#   1. Unique ID for fully combined catalogs\n");
  fprintf(ofp,"#   2. RA in decimal degrees for object (set by the lowest ");
  fprintf(ofp,"numbered catalog in which the object first appears)\n");
  fprintf(ofp,"#   3. Dec in decimal degrees (same setting algorithm)\n");
  fprintf(ofp,"#   4. RA  offset from fiducial object in arcsec\n");
  fprintf(ofp,"#   5. Dec offset from fiducial object in arcsec\n");
  for(i=0; i<nfiles; i++) {
    fprintf(ofp,"#  %2d. ID number in catalog %d (0 ==> no detection)\n",
	    3*i+6,i+1);
    if(i==0)
      fprintf(ofp,"#  %2d. Should be identically zero\n",3*i+7);
    else {
      fprintf(ofp,"#  %2d. Displacement between catalog %d position",3*i+7,
	      i+1);
      fprintf(ofp," and position in first catalog in which object appears\n");
    }
    if(format[i] == 15 || format[i] == 16)
      fprintf(ofp,"#  %2d. Spectroscopic redshift from catalog %d.\n",
	      3*i+8,i+1);
    else
      fprintf(ofp,"#  %2d. Catalog %d representative magnitude\n",3*i+8,i+1);
  }
  fprintf(ofp,"#\n");
  fprintf(ofp,"#Comb                         d_alpha  d_delta  ");
  for(i=0; i<nfiles; i++)
    fprintf(ofp,"        dpos%-2d        ",i+1);
  fprintf(ofp,"\n");
  fprintf(ofp,"# ID       RA        Dec      (arcsec) (arcsec) ");
  for(i=0; i<nfiles; i++) {
    if(format[i] == 15 || format[i] == 16)
      fprintf(ofp," ID%-2d (arcsec)  zspec ",i+1);
    else
      fprintf(ofp," ID%-2d (arcsec)  mag%-2d ",i+1,i+1);
  }
  fprintf(ofp,"\n");
  fprintf(ofp,"#---- ----------- ----------- -------- -------- ");
  for(i=0; i<nfiles; i++)
    if(format[i] == 15 || format[i] == 16)
      fprintf(ofp,"----- -------- ------- ");
    else
      fprintf(ofp,"----- -------- ------ ");
  fprintf(ofp,"\n");

  /*
   * Write main part of output file
   */

  for(i=0,mptr=matchcat; i<ncat; i++,mptr++) {
    fprintf(ofp,"%05d %11.7f %+11.7f %8.2f %8.2f ",
	    mptr->id,mptr->alpha,mptr->delta,mptr->dx,mptr->dy);
    for(j=0; j<nfiles; j++) {
      fprintf(ofp,"%5d %8.2f ",*(id[j]+i),*(dp[j]+i));
      switch(format[j]) {
      case 6: case 7: case 9: case 10:
	fprintf(ofp,"%6.2f ",*(outmag[j]+i));
	break;
      case 15: case 16:
	fprintf(ofp,"%7.4f ",*(outmag[j]+i));
	break;
      default:
	fprintf(ofp,"  0.00 ");
      }
    }
    fprintf(ofp,"\n");
  }

  /*
   * Clean up and exit
   */

  if(ofp)
    fclose(ofp);
  return 0;
}

/*.......................................................................
 *
 * Function help_catcomb
 *
 * Prints useful information for running catcomb
 *
 * Inputs: (none)
 *
 * Output: (none)
 */

void help_catcomb()
{
  char line[MAXC];  /* General string */

  fprintf(stderr,"\nUsage: \n");
  fprintf(stderr,"  catcomb [catfile1] [format1] ... [catfileN] ");
  fprintf(stderr,"[formatN]\n\n");
  fprintf(stderr," ** or **\n\n");
  fprintf(stderr,"  catcomb -c [posfile] [catfile1] [format1] ... [catfileN] ");
  fprintf(stderr,"[formatN]\n\n");
  fprintf(stderr," At least 2 input catalogs are required\n\n");
  fprintf(stderr,"OPTIONS:\n");
  fprintf(stderr," (no flag)   Just match catalogs on positions\n");
  fprintf(stderr,
	  " -c          Match catalogs on positions, but also calculate the\n");
  fprintf(stderr,
	  "              offset of each object from a central position\n");
  fprintf(stderr,"             The central position is given in posfile,\n");
  fprintf(stderr,
	  "              format: label ra_hr ra_min ra_sec dec_deg dec_amin");
  fprintf(stderr," dec_asec\n\n");
  fprintf(stderr,"The format flags indicate the formats of the input ");
  fprintf(stderr,"catalogs.\n");
  fprintf(stderr,"Hit return to see format options: ");
  fgets(line,MAXC,stdin);
  secat_format();

}
