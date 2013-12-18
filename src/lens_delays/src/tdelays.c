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
 * Function declarations
 *
 */

void tdelays_help();

/*.......................................................................
 *
 * Main program
 *
 */

int main(int argc, char *argv[])
{
  int i;                        /* Looping variable */
  int junki;
  int fresult;                  /* Result of running a function */
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
    tdelays_help();
    return 1;
  }

  /*
   * Check the command line arguments and, if the format is OK, create and 
   *  do the initial population of the setup container
   */

  if(!(setup = setup_from_command_line(argv,argc))) {
    tdelays_help();
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

  lc = load_light_curves(setup,npoints);

  /*
   * Temporary check - print out input light curves
   */

  for(i=0; i<ncurves; i++) {
    sprintf(tmpname,"goo_%d.txt",i+1);
    junki = write_fluxrec(lc[i],npoints[i],tmpname,0,0.1);
  }

  /*
   * Get the rest of the setup container parameters
   */

  fresult = get_setup_params(setup,lc);
  if(fresult == ERROR)
    no_error = 0;

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

  if(no_error) {
    printf("\nFinished with tdelays.c\n");
    return 0;
  }
  else {
    fprintf(stderr,"ERROR: Exiting program tdelays.c\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function tdelays_help
 *
 * Gives information on how to run this program
 *
 */

void tdelays_help()
{
  fprintf(stderr,"\n");
  fprintf(stderr,"tdelays.c - Calculates time delays between light curves\n\n");
  fprintf(stderr,"Usage: tdelays flags input_file1 (setup_file)\n");
  fprintf(stderr,"           --- or ----\n");
  fprintf(stderr,"Usage: tdelays flags input_file1 input_file2 (setup_file)\n");
  fprintf(stderr,"\nThe setup file is optional.  ");
  fprintf(stderr,"The input file format is:\n");
  fprintf(stderr,"  day flux1 err1 flux2 err2 (for 1 input file)\n");
  fprintf(stderr,"  day flux err              (for 2 input files)\n\n");
  fprintf(stderr,"Flag  Description\n");
  fprintf(stderr,"----  ----------------\n");
  fprintf(stderr," -1   One input file\n");
  fprintf(stderr," -2   Two input files\n\n");
}
