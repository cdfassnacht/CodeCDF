/*
 * tdelays.c
 *
 * Usage: tdelays [file_1] [file_2]  ([setup_file]),
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
  int no_error=1;               /* Flag set to 0 on error */
  int ncurves=2;                /* Number of input light curves */
  int *index;                   /* Default index for dispersion method */
  char infile[MAXC];            /* Container for input file names */
  char setupfile[MAXC];         /* Name of setup file */
  Fluxrec **lc={NULL};          /* Array of light curves */
  Setup *setup=NULL;            /* Container for setup information */
  LCdisp bestdisp;              /* Container for best-fit results */
  FILE *ofp=NULL;               /* File pointer for optional output file */
 
  /*
   * Check input line
   */

  if(argc < 3) {
    fprintf(stderr,"\nUsage: tdelays [file1] [file2] ");
    fprintf(stderr,"([setup_file]),\n");
    fprintf(stderr,"Each of the input files contains one lightcurve in ");
    fprintf(stderr,"the following format:\n");
    fprintf(stderr,"  day  flux  flux_error\n\n");
    return 1;
  }

  /*
   * Allocate first level of pointers to light curves and initialize
   *  light curve containers.
   */

  lc = (Fluxrec **) malloc(sizeof(Fluxrec *) * ncurves);
  if(!lc) {
    fprintf(stderr,"ERROR:  Insufficient memory for light curve array.\n");
    return 1;
  }
  for(i=0; i<ncurves; i++)
    lc[i] = NULL;

  /*
   * Allocate memory for arrays for number of points and index.
   */

  if(!(index = new_intarray(ncurves,1)))
    no_error = 0;

  /*
   * Initialize the Setup container
   */

  if(no_error) {
    if(!(setup = new_setup(1)))
      no_error = 0;
    else {
      setup->ncurves = ncurves;
      setup->infile[0] = argv[1];
      setup->infile[1] = argv[2];
    }
  }

  /*
   * Read input light curves and set up default index
   */

  for(i=0; i<ncurves; i++) {
    index[i] = i;
    strcpy(infile,argv[i+1]);
    if(!(lc[i] = read_fluxrec_1curve(infile,'#',&setup->npoints[i])))
      no_error = 0;
  }

  /*
   * Hard-wire for dispersion method only, unless override comes
   *  from optional setup file.
   */

  if(no_error) {
    setup->dochi = NO;
    setup->doxcorr = NO;
    setup->doacorr = NO;
    setup->dodisp = YES;
    setup->dodcf = NO;
    setup->docurvefit = NO;
  }

  /*
   * Put setup parameters into setup structure from setup file
   */

  if(argc == 4) {
    strcpy(setupfile,argv[3]);
    if(setup_file(setup,setupfile))
      no_error = 0;
  }

  /*
   * Fill in parts of the setup structure that weren't filled in
   *  from setup file.
   */

  if(no_error)
    if(setup_delays(setup))
      no_error = 0;

  /*
   * Summarize setup parameters
   */

  if(no_error)
    setup_delays_summary(setup);

  /*
   * Find the initial guesses for the flux ratios of the light
   *  curves and the delays between them.  Put these values into the 
   *  Setup container.
   */

  if(no_error) {
    set_tau_grid(lc,setup->npoints,index,setup);
    set_mu_grid(lc,setup->npoints,setup);
  }
#if 0
  /*
   * Calculate chi-squared with shifts
   */

  if(setup->dochi && no_error)
    if(do_chi(lc1,lc2,n1,n2,setup,chilogfp,1,SLOW) == 1)
      no_error = 0;

  /*
   * Calculate cross-correlation between curves
   */

  if(setup->doxcorr && no_error) {
    if(CORRFUNC == 3) {
      if(call_xcorr_time(lc1,lc1,n1,n1,setup,xclogfp,
			 "corrout.dat",1))
	no_error = 0;
    }
    else {
      if(call_xcorr_fft(lc1,n1,"corrout.dat",1,xclogfp,1))
	no_error = 0;
    }
  }

  /*
   * Calculate discrete correlation function curves
   */

  if(setup->dodcf && no_error)
    if(call_dcf(lc1,n1,"dcfout.dat",NULL,1))
      no_error = 0;
#endif
  /*
   * Dispersion spectrum calculation
   */

  if(setup->dodisp && no_error)
    if(disp_setup(lc,2,setup->npoints,index,setup,&bestdisp,"disp.out",1))
      no_error = 0;
    else if(setup->outfile) {
      printf("\n");
      ofp = open_appendfile(setup->outfile,1);
      fprintf(ofp,"%s %s %6.2f %6.4f %f\n",setup->infile[0],setup->infile[1],
	     bestdisp.tau,bestdisp.mu,bestdisp.disp);
      fclose(ofp);
    }

  /*
   * Clean up
   */

  if(no_error)
    printf("\nCleaning up\n");

  for(i=0; i<ncurves; i++) {
    lc[i] = del_fluxrec(lc[i]);
  }
  if(lc)
    free(lc);
  setup = del_setup(setup);
  index = del_intarray(index);

#if 0
  if(chilogfp)
    fclose(chilogfp);
  if(xclogfp)
    fclose(xclogfp);
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
