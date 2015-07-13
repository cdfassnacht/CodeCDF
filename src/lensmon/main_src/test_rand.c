/*
 * test_rand.c
 *
 * Usage: test_rand
 *
 * Generates a time series of zero-mean Gaussian-distributed random noise.
 *
 * 11Sep98 CDF
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "structdef.h"
#include "dataio.h"
#include "nr.h"
#include "lc_funcs.h"

/*.......................................................................
 *
 * Main program
 *
 */

int main(int argc, char *argv[])
{
  int i,j,k;                  /* Looping variables */
  int no_error=1;             /* Flag set to 0 on error */
  int nlines;                 /* Number of data lines in lens input file */
  long randseed=-97867564;    /* Initial seed for random number generator */
  char tmname[MAXC];          /* Name for template curve files */
  char outname[MAXC];         /* Output file name */
  char lensfile[MAXC];        /* Name of 1608 data file */
  Fluxrec *fl[N08]={NULL};    /* 1608 light curves */
  Fluxrec *aptr;              /* Pointer to navigate fl[0] */
  Fluxrec *bptr;              /* Pointer to navigate fl[1] */
  Fluxrec *cptr;              /* Pointer to navigate fl[2] */
  Fluxrec *dptr;              /* Pointer to navigate fl[3] */
  FILE *ofp;                  /* Output file pointer */
  FILE *lfp=NULL;             /* 1608 data file pointer */

  /*
   * Read in templates
   */
#if 0
  for(i=0; i<N08; i++) {
    sprintf(tmname,"smooth%d",i);
    if(!(fl[i] = read_fluxrec(tmname,'#',&nlines)))
      no_error = 0;
  }
#endif
  sprintf(lensfile,"mc_g.00001");
  while((lfp = fopen(lensfile,"r")) == NULL) {
    printf("ERROR:  %s does not exist.\n",lensfile);
    printf("  Enter input file name again:  ");
    gets(lensfile);
  }

  /*
   * First count the number of data input lines in order to set
   *  sizes of arrays.
   */

  nlines = n_lines(lfp,'#');

  if(nlines == 0) {
    fprintf(stderr,"ERROR. No valid data lines.\n");
    no_error = 0;
  }

  /*
   * Put pointers back to beginning of input files
   */

  if(no_error)
    rewind(lfp);

  /*
   * Allocate memory for random noise curve containers
   */

  for(i=0; i<N08; i++)
    if(!(fl[i] = new_fluxrec(nlines)))
      no_error = 0;

  j = 0;
  sprintf(outname,"randgauss%05d",j+1);

  /*
   * Read in data
   */

  if(no_error)
    if(read_1608(fl,nlines,lfp,lensfile))
      no_error = 0;

  for(j=0; j<1000; j++) {
    randseed -= (j+1);

  /*
   * Open output file
   */

    sprintf(outname,"randgauss%05d",j+1);
    if(!(ofp = fopen(outname,"w")))
      no_error = 0;

    /*
     * Fill random noise containers
     */

    if(no_error) {
      for(i=0; i<N08; i++) {
	if(gauss_noise(fl[i],nlines,0.02,&randseed))
	  no_error = 0;
      }
    }
#if 0
    /*
     * Do the cross-correlation
     */

    if(call_xcorr_time(fl,fl,nlines,nlines,NULL,ofp,"corrout.dat",1))
      no_error = 0;
#endif
    /*
     * Print output
     */

    if(no_error) {
      for(i=0,aptr=fl[0],bptr=fl[1],cptr=fl[2],dptr=fl[3];
	  i<nlines; i++,aptr++,bptr++,cptr++,dptr++)
	fprintf(ofp,"%6.1f %9f %9f %9f %9f %5f %5f %5f %5f\n",aptr->day,
		aptr->flux,bptr->flux,cptr->flux,dptr->flux,
		aptr->err,bptr->err,cptr->err,dptr->err);
    }

    if(ofp)
      fclose(ofp);
  }

  /*
   * Clean up
   */

  for(i=0; i<N08; i++) {
    fl[i] = del_fluxrec(fl[i]);
  }

  if(lfp)
    fclose(lfp);
  if(ofp)
    fclose(ofp);

  if(no_error)
    return 0;
  else {
    fprintf(stderr,"ERROR: Exiting program\n");
    return 1;
  }
}
