/*
 * tdelays.c
 *
 * Usage: tdelays ***NEEDS TO BE FINALIZED***  ([setup_file]),
 *   where each input file contains a light curve in the format:
 *    day  flux error
 *   The setup file is optional.
 *
 * Finds delays between light curves.  For now, the only method allowed
 *  is  the Pelt et al. dispersion method.
 *
 * 02Sep2002 CDF,   A modification of delays.c
 * v01Aug2005 CDF,  Slight tweaks to make more user-friendly
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "structdef.h"
#include "dataio.h"
#include "lc_setup.h"
#include "lc_funcs.h"
#include "noninterp_fns.h"

/*.......................................................................
 *
 * Main program
 *
 */

int main(int argc, char *argv[])
{
  int i;                        /* Looping variable */
  int junki;
  int no_error=1;               /* Flag set to 0 on error */
  int ncurves=2;                /* Number of input light curves */
  int *npoints=NULL;            /* Number of points in each light curve */
  int *index;                   /* Default index for dispersion method */
  char tmpname[MAXC];
  char setupfile[MAXC];         /* Name of setup file */
  Fluxrec **lc={NULL};          /* Array of light curves */
  Setup *setup=NULL;            /* Container for setup information */
  LCdisp bestdisp;              /* Container for best-fit results */
  FILE *ofp=NULL;               /* File pointer for optional output file */
 
  /*
   * Check input line
   */

  if(argc < 3) {
    fprintf(stderr,"\nUsage: tdelays ***TO BE DETERMINED*** ");
    fprintf(stderr,"([setup_file]),\n");
    fprintf(stderr,"Each of the input files contains one lightcurve in ");
    fprintf(stderr,"the following format:\n");
    fprintf(stderr,"  day  flux  flux_error\n\n");
    return 1;
  }

  /*
   * Allocate memory for arrays for number of points and index.
   */

  if(!(npoints = new_intarray(ncurves,1)))
    no_error = 0;
#if 0
  if(!(index = new_intarray(ncurves,1)))
    no_error = 0;
#endif

  /*
   * Load the light curves
   */

  lc = load_light_curves(argc,argv,ncurves,npoints);

  /*
   * Temporary check - print out input light curves
   */

  for(i=0; i<ncurves; i++) {
    sprintf(tmpname,"foo_%d.txt",i+1);
    junki = write_fluxrec(lc[i],npoints[i],tmpname,0,0.1);
  }
#if 0

  /*
   * Clean up and return
   */

  if(no_error)
    printf("\nCleaning up\n");

  for(i=0; i<ncurves; i++) {
    lc[i] = del_fluxrec(lc[i]);
  }
  if(lc)
    free(lc);
#endif

  if(no_error) {
    printf("\nFinished with tdelays.c\n");
    return 0;
  }
  else {
    fprintf(stderr,"ERROR: Exiting program tdelays.c\n");
    return 1;
  }
}
