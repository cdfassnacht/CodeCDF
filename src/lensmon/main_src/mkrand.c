/*
 * mkrand.c
 *
 * Usage: lcurve [1608 filename] [1634/1635 filename]
 *
 * Takes the observed light curves and randomizes the days while
 *  keeping the fluxes to create new light curves with the same flux
 *  distributions.  The main work in this program is done by a call
 *  to rand_curves, which is in the monte.c library.
 *
 * 22Jun98 CDF
 * v25Aug98 CDF, Flat-field 1608 curves before passing them to rand_curves.
 *               Modified the bad-day flagging to take account of new
 *                versions of function in lc_funcs.c.
 * v05Sep98 CDF, Inserted improved bad-day flagging.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "structdef.h"
#include "dataio.h"
#include "nr.h"
#include "lc_funcs.h"
#include "monte.h"

/*.......................................................................
 *
 * Main program
 *
 */

int main(int argc, char *argv[])
{
  int i;                      /* Looping variable */
  int no_error=1;             /* Flag set to 0 on error */
  int nlines;                 /* Number of data lines in lens input file */
  int nclines;                /* Number of data lines in comp input file */
  int nbad[N08];              /* Number of bad days in bad day file */
  int ncurves=10000;          /* Number of randomized curves to create */
  int flagbad=0;              /* Flag set to 1 for bad day flagging */
  float *flat=NULL;           /* Light curve "flat field" */
  float fracrms;              /* Fractional rms scatter in 1634/1635 ratio */
  char lensfile[MAXC];        /* Input file name for 1608 data */
  char compfile[MAXC];        /* Input file name for 1634 and 1635 data */
  char line[MAXC];            /* General string for reading input */
  Fluxrec *fl08[N08]={NULL};  /* 1608 light curves */
  Fluxrec *ff08[N08]={NULL};  /* Flat-fielded 1608 light curves */
  Fluxrec *fl34=NULL;         /* 1634 light curve */
  Fluxrec *fl35=NULL;         /* 1635 light curve */
  FILE *lfp=NULL;             /* 1608 data file pointer */
  FILE *cfp=NULL;             /* 1634/1635 data file pointer */

  /*
   * Check input line
   */

  if(argc != 3) {
    fprintf(stderr,"\nUsage: mkrand [1608 filename] [1634/1635 filename]\n\n");
    return 1;
  }

  /*
   * Check to see that the model files can be opened.
   */

  strcpy(lensfile,argv[1]);
  while((lfp = fopen(lensfile,"r")) == NULL) {
    printf("ERROR:  %s does not exist.\n",lensfile);
    printf("  Enter input file name again:  ");
    gets(lensfile);
  }

  strcpy(compfile,argv[2]);
  while((cfp = fopen(compfile,"r")) == NULL) {
    printf("ERROR:  %s does not exist.\n",compfile);
    printf("  Enter input file name again:  ");
    gets(compfile);
  }


  /*
   * First count the number of data input lines in order to set
   *  sizes of arrays.
   */

  nlines = n_lines(lfp,'#');
  nclines = n_lines(cfp,'#');

  if(nlines == 0 || nclines == 0) {
    fprintf(stderr,"ERROR. No valid data lines.\n");
    no_error = 0;
  }
  else if(nlines != nclines) {
    fprintf(stderr,
	    "ERROR. Number of lines in lens and comp files don't match.\n");
    no_error = 0;
  }
  else {
    printf("\n%d data lines in lens and comp input files.\n",nlines);
  }

  /*
   * Put pointers back to beginning of input files
   */

  if(no_error) {
    rewind(lfp);
    rewind(cfp);
  }

  /*
   * Allocate arrays
   */

  if(no_error) {
    for(i=0; i<N08; i++) {
      if(!(fl08[i] = new_fluxrec(nlines)))
	no_error = 0;
    }
  }

  if(no_error)
    if(!(fl34 = new_fluxrec(nlines)))
      no_error = 0;

  if(no_error)
    if(!(fl35 = new_fluxrec(nlines)))
      no_error = 0;

  /*
   * Read in data
   */

  if(no_error)
    if(read_data(fl34,fl35,fl08,nlines,lfp,cfp))
      no_error = 0;

  /*
   * See if bad-day flagging is requested
   */

  if(no_error) {
    printf("Flag bad days? (1/0) [%d] ",flagbad);
    gets(line);
    if(strcmp(line,"") != 0) {
      while(sscanf(line,"%d",&flagbad) != 1 || flagbad < 0) {
	fprintf(stderr,"ERROR: bad input.  Enter value again:  ");
	gets(line);
      }
    }
  }

  /*
   * Flag bad days, if requested
   */

  if(flagbad && no_error) {
    if(flag_bad(fl34,fl35,fl08,nlines,&flagbad,"badday.list",nbad)) {
      no_error = 0;
    }
    else {
      if(flagbad == 2)
	printf("Bad days flagged: %d %d %d %d\n",nbad[0],nbad[1],nbad[2],
	       nbad[3]);
      else
	printf("Bad days flagged: %d\n\n",nbad[0]);
    }
  }

  /*
   * Get number of curves to create
   */

  if(no_error) {
    printf("Enter number of randomized curves to create: [%d] ",ncurves);
    gets(line);
    if(strcmp(line,"") != 0) {
      while(sscanf(line,"%d",&ncurves) != 1 && ncurves < 0) {
	fprintf(stderr,"ERROR: Bad input value.  Enter number again:  ");
	gets(line);
      }
    }
  }

  /*
   * Make the "flat field"
   */

  if(no_error)
    if(!(flat = make_flat(fl34,fl35,nlines)))
      no_error = 0;

  /*
   * Calculate fractional error in 1634/1635 flux ratio
   */

  if(no_error)
    if(ratio_err(fl34,fl35,nlines,&fracrms))
      no_error = 0;

  /*
   * "Flat-field" the data
   */

  if(no_error) {
    for(i=0; i<N08; i++) {
      if(!(ff08[i] = flat_field(fl08[i],flat,nlines,fracrms,NULL,0)))
	no_error = 0;
    }
  }

  /*
   * Create random curves
   */

  if(no_error)
    rand_curves(ff08,nlines,nbad[0],ncurves);

  /*
   * Clean up
   */

  printf("\nCleaning up\n\n");

  fl34 = del_fluxrec(fl34);
  fl35 = del_fluxrec(fl35);
  for(i=0; i<4; i++) {
    fl08[i] = del_fluxrec(fl08[i]);
    ff08[i] = del_fluxrec(fl08[i]);
  }
  flat = del_array(flat);

  if(lfp)
    fclose(lfp);
  if(cfp)
    fclose(cfp);

  if(no_error)
    return 0;
  else {
    fprintf(stderr,"ERROR: Exiting program\n");
    return 1;
  }
}
