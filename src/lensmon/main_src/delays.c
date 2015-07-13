/*
 * delays.c
 *
 * Usage: delays [1608_raw file] [1608_interp_file]  ([setup_file]),
 *   where the 1608 interpolated lightcurve has already been interpolated 
 *   onto a regular grid (for the chisq and cross-correlation methods) and
 *   the setup file is optional.
 *
 * Finds delays between light curves that have already been "flat-fielded"
 *  AND smoothed/interpolated.  The delays can be found by numerous methods,
 *  including chisq-minimization, cross-correlation, and the Pelt et al.
 *  dispersion method.
 * NB: This program takes over the delay-finding functions of the former
 *  lcurve.c and lcurve_mc.c.  The interpolation functionality of those
 *  programs is now contained in interp.c.  For the methods needing
 *  an interpolated curve (e.g. chisq minimization, cross-correlation)
 *  it is necessary to run interp.c before running this program.
 *
 * 09Jul99 CDF,  A modification of lcurve_mc.c
 * v16Sep99 CDF, Slight modifications to deal with new version of lc_setup.c
 * v06Sep00 CDF, Changed to take two input light curves instead of one.
 *                The first can be a raw light curve while the second one
 *                must be interpolated.  This allows the proper running
 *                of the do_chi function.
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
#include "lc_chisq.h"
#include "lc_interp.h"
#include "correlate.h"
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
  int nraw=0;                   /* Number of points in raw curve */
  int ninterp=0;                /* Number of points in interpolated curve */
  int nbad[N08]={0,0,0,0};      /* Number of bad days in bad day file */
  char rawfile[MAXC];           /* Input file name for raw data */
  char interpfile[MAXC];        /* Input file name for interpolated data */
  char setupfile[MAXC];         /* Name of setup file */
  Fluxrec *raw[N08]={NULL};     /* Raw 1608 light curves */
  Fluxrec *interp[N08]={NULL};  /* Interpolated 1608 light curves */
  Setup *setup=NULL;            /* Container for setup information */
  FILE *rfp=NULL;               /* File pointer for raw data file */
  FILE *ifp=NULL;               /* File pointer for interpolated data file */
  FILE *chilogfp=NULL;          /* Chisq logfile pointer */
  FILE *xclogfp=NULL;           /* Cross-correlation logfile pointer */
 
  /*
   * Check input line
   */

  if(argc < 3) {
    fprintf(stderr,"\nUsage: delays [1608_raw_file] [1608_interp_file] ");
    fprintf(stderr,"([setup_file]),\n");
    fprintf(stderr," where the second 1608 curve has already been smoothed ");
    fprintf(stderr," and/or interpolated \n");
    fprintf(stderr," (for the chisq and xcorr methods) and the setup file is");
    fprintf(stderr," optional.\n\n");
    fprintf(stderr,"NB: The interpolation, if necessary, should be carried ");
    fprintf(stderr,"out by running interp.c\n");
    fprintf(stderr," BEFORE running this program.\n\n");
    return 1;
  }

  /*
   * Initialize the Setup container
   */

  if(!(setup = new_setup(1))) {
    fprintf(stderr,"ERROR. Exiting program.\n");
    return 1;
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
   * Open lens light curve files 
   */

  strcpy(rawfile,argv[1]);
  strcpy(interpfile,argv[2]);
  if(!(rfp = open_readfile(rawfile)))
    no_error = 0;
  if(!(ifp = open_readfile(interpfile)))
    no_error = 0;

  /*
   * Count the number of data input lines in order to set
   *  sizes of arrays.
   */

  if(no_error) {
    nraw = n_lines(rfp,'#');
    ninterp = n_lines(ifp,'#');
  }

  if(nraw == 0) {
    fprintf(stderr,"ERROR. No valid data lines in %s.\n",rawfile);
    no_error = 0;
  }
  else if(ninterp == 0) {
    fprintf(stderr,"ERROR. No valid data lines in %s.\n",interpfile);
    no_error = 0;
  }
  else {
    rewind(rfp);
    rewind(ifp);
  }

  /*
   * Allocate arrays
   */

  if(no_error) {
    for(i=0; i<N08; i++) {
      if(!(raw[i] = new_fluxrec(nraw)))
	no_error = 0;
      if(!(interp[i] = new_fluxrec(ninterp)))
	no_error = 0;
    }
  }

  /*
   * Read in data
   */

  if(no_error) {
    if(read_1608(raw,nraw,rfp,rawfile))
      no_error = 0;
    if(read_1608(interp,ninterp,ifp,interpfile))
      no_error = 0;
  }

  /*
   * Open log files if requested
   */

  if(no_error) {
    if(strcmp(setup->chilog,"stdout") != 0) {
      if(!(chilogfp = open_appendfile(setup->chilog,1))) {
	fprintf(stderr,"ERROR.  Cannot open chisq log file.  ");
	fprintf(stderr,"Output will be sent to stdout\n");
      }
    }
  }
  if(no_error) {
    if(strcmp(setup->xclog,"stdout") != 0) {
      if(!(xclogfp = open_appendfile(setup->xclog,1))) {
	fprintf(stderr,"ERROR.  Cannot open cross-correlation log file.  ");
	fprintf(stderr,"Output will be sent to stdout\n");
      }
    }
  }

  /*
   * Calculate the initial guesses for the flux ratios of the light
   *  curves and put them into the Setup container.
   */

  if(no_error)
    set_mu0(raw,nraw,setup);

  /*
   * Calculate chi-squared with shifts
   */

  if(setup->dochi && no_error)
    if(do_chi(raw,interp,nraw,ninterp,setup,chilogfp,1,SLOW) == 1)
      no_error = 0;

  /*
   * Calculate auto-correlation between curves to get normalization
   */
#if 0
  printf("\nCalculating auto-correlations....\n\n");
  if(no_error)
    if(call_acorr(raw,nraw,"acorr.dat",0))
      no_error = 0;
#endif
  /*
   * Calculate cross-correlation between curves
   */

  if(setup->doxcorr && no_error) {
    if(CORRFUNC == 3) {
      if(call_xcorr_time(raw,raw,nraw,nraw,setup,xclogfp,
			 "corrout.dat",1))
	no_error = 0;
    }
    else {
      if(call_xcorr_fft(raw,nraw,"corrout.dat",1,xclogfp,1))
	no_error = 0;
    }
  }

  /*
   * Calculate discrete correlation function curves
   */

  if(setup->dodcf && no_error)
    if(call_dcf(raw,nraw,"dcfout.dat",NULL,1))
      no_error = 0;

  /*
   * Dispersion spectrum calculation
   */

  if(setup->dodisp && no_error)
    if(call_disp(raw,nraw,setup,1))
      no_error = 0;

  /*
   * Simultaneous curve-fitting and chisq minimization
   */

  if(setup->docurvefit && no_error)
    if(fit_1608(raw,nraw,nbad,1.0))
      no_error = 0;

  /*
   * Clean up
   */

  if(no_error)
    printf("\nCleaning up\n\n");

  for(i=0; i<4; i++) {
    raw[i] = del_fluxrec(raw[i]);
    interp[i] = del_fluxrec(interp[i]);
  }
  setup = del_setup(setup);

  if(rfp)
    fclose(rfp);
  if(ifp)
    fclose(ifp);
  if(chilogfp)
    fclose(chilogfp);
  if(xclogfp)
    fclose(xclogfp);

  if(no_error)
    return 0;
  else {
    fprintf(stderr,"ERROR: Exiting program delays.c\n");
    return 1;
  }
}
