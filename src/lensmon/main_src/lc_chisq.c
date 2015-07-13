/*
 * lc_chisq.c
 *
 * A library of functions to perform chisq minimization on pairs of
 *  1608 light curves.
 *
 * 27Jul00 CDF,  Split off from lc_funcs.c
 * v30Jul00 CDF, Changed step_delays (*** BUT NOT STEP_DELAYS_2 ***) to
 *                search from a min to a max delay, rather than the whole
 *                range of delays.  This necessitates the creation of the
 *                new Prange (for parameter range) structure, to keep the
 *                number of passed variables down.  This also requires
 *                many changes to shift_chi.
 * v07Sep00 CDF, Changed shift_chi to search a pre-set range of 
 *                magnifications rather than just taking the default values.
 * v08Sep00 CDF, Added a new structure definition, LCchisq, to hold
 *                the data from the chisq search grid.  This required the
 *                addition of new functions new_lcchisq, del_lcchisq, and
 *                read_lcchisq.
 *               Fixed a bug in step_delays and modified shift_chi to 
 *                contain a flag for parabola fitting.
 * v11Sep00 CDF, Changed the structure type from Fluxrec to LCchisq in
 *                several function in which this new structure is much
 *                more suitable.  Functions affected are: do_chi, shift_chi,
 *                step_delays, find_min_chisq.
 *               *** Got rid of all the _2 functions since development on
 *                   them had stopped and it was becoming too hard to keep
 *                   things parallel. ***
 *               Moved chisq-slice printing into new function 
 *                print_chisq_slice.
 * v30Sep00 CDF, Added more output flexibility by moving output files for
 *                shift_chi from hard-wired values to variable names 
 *                contained in Setup structure.
 * v10Oct00 CDF, Fixed small bugs in print_chisq_slice.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "structdef.h"
#include "dataio.h"
#include "nrutil.h"
#include "nr.h"
#include "lc_setup.h"
#include "lc_funcs.h"
#include "lc_chisq.h"

/*.......................................................................
 *
 * Function new_lcchisq
 *
 * Allocates dynamic memory for a pointer array of LCchisq structures
 *
 * Input:  int size            size of the array
 *
 * Output: LCchisq *newstruct  pointer to the new array.  NULL if error
 *
 */

LCchisq *new_lcchisq(int size)
{
  int i;
  LCchisq *newstruct;
  LCchisq *lcptr;
 
  newstruct = (LCchisq *) malloc(sizeof(LCchisq) * size);
  if(!newstruct) {
    fprintf(stderr,"new_lcchisq:  Insufficient memory for array.\n");
    return NULL;
  }

  /*
   * Initialize values in the array
   */

  for(i=0,lcptr=newstruct; i<size; i++,lcptr++) {
    lcptr->tau = lcptr->mu = lcptr->chisq = lcptr->rchisq = 0.0;
    lcptr->ndof = lcptr->arraypos = 0;
  }

  return newstruct;
}

/*.......................................................................
 *
 * Function del_lcchisq
 *
 * Frees up memory allocated to Lcchisq array
 *
 * Input:  LCchisq *lcchisq    array to be freed
 *
 * Output: NULL
 */

LCchisq *del_lcchisq(LCchisq *lcchisq)
{
  if(lcchisq)
    free(lcchisq);
 
  return NULL;
}

/*.......................................................................
 *
 * Function read_lcchisq
 *
 * Reads a grid of LCchisq structures from one of the chi*dat files
 *  produced by delays.c.  The filename is passed as one of the
 *  arguments.
 *
 * Inputs: char *filename      input filename
 *         Axisinfo *tau       information about delay values
 *         Axisinfo *mu        information about magnification values
 *
 * Output: LCchisq *lcchisq    filled LCchisq structure.  NULL on error.
 *
 */

LCchisq *read_lcchisq(char *filename, Axisinfo *tau, Axisinfo *mu)
{
  int no_error=1;         /* Flag set to 0 on error */
  int nlines;             /* Number of data lines in input file */
  char line[MAXC];        /* General string variable for getting input */
  LCchisq *lcchisq=NULL;  /* LCchisq array to be filled */
  LCchisq *lcptr;         /* Pointer to navigate lcchisq */
  FILE *ifp=NULL;         /* Input file pointer */

  /*
   * Open input file
   */

  if(!(ifp = open_readfile(filename))) {
    fprintf(stderr,"ERROR: read_lcchisq\n");
    return NULL;
  }

  /*
   * Count number of lines in input file
   */

  if((nlines = n_lines(ifp,'#')) == 0) {
    fprintf(stderr,"ERROR: read_lcchisq. No data in input file %s.\n",
	    filename);
    no_error = 0;
  }
  else {
    printf("read_lcchisq: %s contains %d data lines.\n",filename,nlines);
    rewind(ifp);
  }

  /*
   * Allocate memory for LCchisq array
   */

  if(!(lcchisq = new_lcchisq(nlines)))
    no_error = 0;

  /*
   * Read in header information, which should be contained on the
   *  first line of the input file.
   */

  if(no_error) {
    fgets(line,MAXC,ifp);
    if(sscanf(line,"# Gridsize %d %d",&tau->nval,&mu->nval) != 2) {
      fprintf(stderr,"ERROR: read_lcchisq.  %s missing header information.\n",
	      filename);
      no_error = 0;
    }
    else
      printf("read_lcchisq: Input array from %s has dimensions %d x %d\n",
	     filename,tau->nval,mu->nval);
  }

  /*
   * Read in data
   */

  if(no_error) {
    tau->minval = mu->minval = 9999.9;
    tau->maxval = mu->maxval = -9999.9;
    lcptr = lcchisq;
    while(fgets(line,MAXC,ifp) != NULL && no_error) {
      if(line[0] != '#') {
	if(sscanf(line,"%f %f %f %d %f",&lcptr->tau,&lcptr->mu,
		  &lcptr->chisq,&lcptr->ndof,&lcptr->rchisq) != 5) {
	  fprintf(stderr,"ERROR: read_lcchisq. Bad input format.\n");
	  no_error = 0;
	}
	else {
	  if(lcptr->tau > tau->maxval)
	    tau->maxval = lcptr->tau;
	  if(lcptr->tau < tau->minval)
	    tau->minval = lcptr->tau;
	  if(lcptr->mu > mu->maxval)
	    mu->maxval = lcptr->mu;
	  if(lcptr->mu < mu->minval)
	    mu->minval = lcptr->mu;
	  lcptr++;
	}
      }
    }
  }

  /*
   * Print out axis info
   */

  if(no_error) {
    printf("read_lcchisq: tau_min = %7.2f, tau_max = %7.2f, ntau = %d\n",
	   tau->minval,tau->maxval,tau->nval);
    printf("read_lcchisq: mu_min  =  %6.4f, mu_max  =  %6.4f, nmu  = %d\n",
	   mu->minval,mu->maxval,mu->nval);
  }

  /*
   * Clean up and exit
   */

  if(ifp)
    fclose(ifp);

  if(no_error)
    return lcchisq;
  else {
    fprintf(stderr,"ERROR: read_lcchisq.\n");
    return del_lcchisq(lcchisq);
  }
}

/*.......................................................................
 *
 * Function do_chi
 *
 * A "meta-function" to call shift_chi a number of times.
 *
 * Inputs: Fluxrec *raw[]      control light curves
 *         Fluxrec *interp[]   comparsion light curves
 *         int nraw            number of values in control light curves
 *         int nint            number of values in comparison light curves
 *         Setup *setup        smoothing/interpolation info
 *         FILE *logfp         logfile pointer
 *         int doprint         flag set to 1 if output files are desired
 *         int speed           if the input curves are of the same size and
 *                              have the same spacing, this is set to FAST
 *                              and a quicker form of calculation can occur.
 *                              Otherwise, it will be set to SLOW.
  *
 * Output: int (0 or 1)        0 ==> success, 1 ==> failure
 *
 * v25Jun98 CDF, Made the interpolated curves into an array of Fluxrec
 *                arrays.
 * v09Jul98 CDF, Moved the interpolation of the light curves up into the
 *                calling function (e.g. lcurve.c).
 *               Changed output format to put all logfile output on same line.
 *               Added doprint flag as a passed parameter to control whether
 *                output gets printed.
 * v13Jul98 CDF, Better output format
 * v23Aug98 CDF, Cleaner output code.
 * v30Aug98 CDF, Allow two types of input curves: raw (sparsely sampled) and
 *                smoothed/interpolated.  Before both curves were interpolated.
 * v30Jan99 CDF, Added possibility of passing initial guess for flux density
 *                ratios to shift_chi.
 * v11Sep00 CDF, Changed the best* variables from Fluxrec to LCchisq
 *                structures.
 * v30Sep00 CDF, Changed shift_chi output files from hard-wired chiba.dat, etc.
 *                to values contained in the setup structure.
 */

int do_chi(Fluxrec *raw[], Fluxrec *interp[], int nraw, int nint, 
	   Setup *setup, FILE *logfp, int doprint, int speed)
{
  LCchisq bestba;            /* B-A lag with lowest chisq */
  LCchisq bestab;            /* A-B lag with lowest chisq */
  LCchisq bestbc;            /* B-C lag with lowest chisq */
  LCchisq bestcb;            /* C-B lag with lowest chisq */
  LCchisq bestbd;            /* B-C lag with lowest chisq */
  LCchisq bestdb;            /* D-B lag with lowest chisq */

  /*
   * Call shift_chi with input light curves, calculating chisq for curves
   *  A, C, and D with respect to curve B.  In each call to shift_chi
   *  the raw light curve will be the "control" curve.
   *
   * First call B and A
   */

  printf("\nCalculating chisq curves....\n\n");

  if(shift_chi(raw[1],interp[0],nraw,nint,setup->tau0[0],setup->mu0[0],setup,
	       speed,&bestba,setup->achifile,doprint)) {
    fprintf(stderr,"ERROR: do_chi\n");
    return 1;
  }

  /*
   * B and C
   */
  
  if(shift_chi(raw[1],interp[2],nraw,nint,setup->tau0[2],setup->mu0[2],setup,
	       speed,&bestbc,setup->cchifile,doprint)) {
    fprintf(stderr,"ERROR: do_chi\n");
    return 1;
  }

  /*
   * B and D
   */

  if(shift_chi(raw[1],interp[3],nraw,nint,setup->tau0[3],setup->mu0[3],setup,
	       speed,&bestbd,setup->dchifile,doprint)) {
    fprintf(stderr,"ERROR: do_chi\n");
    return 1;
  }

  printf(" shift_chi: ---------------------------------------------------");
  printf("---------------\n");

#if 0
  /*
   * Now adjust the day and flux values that need to be changed to
   *  produce a uniform output.
   */

  bestab.tau *= -1.0;
  bestcb.tau *= -1.0;
  bestdb.tau *= -1.0;
  bestab.mu = 1.0 / bestab.mu;
  bestcb.mu = 1.0 / bestcb.mu;
  bestdb.mu = 1.0 / bestdb.mu;
#endif
  /*
   * Print out results
   */

  printf("\n do_chi: Best-fit values:\n");
  printf("         Pair     lag     flux ratio    chisq\n");
  printf("         -----  -------   ----------   --------\n");
  printf("          B-A  %7.2f      %6.4f     %8.4f\n",bestba.tau,
	 bestba.mu,bestba.rchisq);
  printf("          B-C  %7.2f      %6.4f     %8.4f\n",bestbc.tau,
	 bestbc.mu,bestbc.rchisq);
  printf("          B-D  %7.2f      %6.4f     %8.4f\n",bestbd.tau,
	 bestbd.mu,bestbd.rchisq);

  if(logfp) {
#if 0
    fprintf(logfp,"# Chisq results\n");
    fprintf(logfp,"#         Pair     lag     flux ratio    chisq\n");
    fprintf(logfp,"#         -----  -------   ----------   --------\n");
    fprintf(logfp,"#          B-A  %7.2f      %6.4f     %8.4f\n",bestba.tau,
	    bestba.mu,bestba.rchisq);
    fprintf(logfp,"#          B-C  %7.2f      %6.4f     %8.4f\n",bestbc.tau,
	    bestbc.mu,bestbc.rchisq);
    fprintf(logfp,"#          B-D  %7.2f      %6.4f     %8.4f\n",bestbd.tau,
	    bestbd.mu,bestbd.rchisq);
#endif
    fprintf(logfp,"%7.2f %6.4f %8.4f ",bestba.tau,bestba.mu,bestba.rchisq);
    fprintf(logfp,"%7.2f %6.4f %8.4f ",bestbc.tau,bestbc.mu,bestbc.rchisq);
    fprintf(logfp,"%7.2f %6.4f %8.4f ",bestbd.tau,bestbd.mu,bestbd.rchisq);
    fprintf(logfp,"\n");
  }

  return 0;
}

/*.......................................................................
 *
 * Function shift_chi
 *
 * Computes the chi-squared obtained by comparing two light curves, the
 *  "control" one, with measured values, and the "comparison" one, with
 *  interpolated values.  The comparison curve is shifted by different
 *  lags and flux ratios, and the reduced chi-squared is computed at each 
 *  (lag, flux-ratio) pair.
 *
 * Inputs: Fluxrec *control    control light curve
 *         Fluxrec *comp       comparison light curve
 *         int nctrl           number of points in control curve
 *         int ncomp           number of points in comparison curve
 *         float tau0          initial guess for delay
 *         float mu0           initial guess at flux density ratio (calculated
 *                              by this function if mu0 == 0.0)
 *         Setup *setup        smoothing/interpolation info
 *         Fluxrec *bestlag    lag with best-fit chisq (set by this function)
 *         char *outname       name of output file
 *         int doprint         flag set to 1 if output files are desired
 *         int speed           if the input curves are of the same size and
 *                              have the same spacing, this is set to FAST
 *                              and a quicker form of calculation can occur.
 *                              Otherwise, it will be set to SLOW.
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 * v12Feb97 CDF, Moved comparison curve normalization to after the
 *                point where the interpolated curve has been shifted.
 *               Only include points in dummy array that match points
 *                in comparison curve.
 * v10Jun98 CDF, Transfer interpolation of light curves into this function
 *                from main.
 *               Put in steps in flux ratio rather than just normalizing
 *                the two curves.
 *               Moved the printing to the output file into this function.
 * v16Jun98 CDF, Added flag to choose interpolation or not
 * v19Jun98 CDF, Moved interpolation up one level to do_chi.
 *               Added the return of the lag with best-fitting chisq
 * v09Jul98 CDF, Added doprint flag as a passed parameter to control whether
 *                output gets printed.
 * v19Jul98 CDF, Changed steps in flux density ratio to be fractions of
 *                ratio0 rather than absolute values.
 * v23Aug98 CDF, Added a call to find_min_chisq which fits a parabola to the
 *                lag-chisq curve at the grid points with the lowest chisq
 *                to find the "true" minimum chisq value.
 * v30Aug98 CDF, Changed to allow control and comp to be of different sizes,
 *                e.g. if the control is a raw curve and the comparison curve
 *                has been smoothed.
 * v24Sep98 CDF, Change N_dof calculation to use the estimated number of 
 *                independent points in the overlap region rather than
 *                just the number of overlapping points (this only applies
 *                in the case in which both curves are interpolted, i.e.
 *                when speed == FAST).
 * v30Jan99 CDF, Moved delay loop into step_delays.
 *               Changed ratio0 from mean(control)/mean(comp) to
 *                mean(comp)/mean(control).
 *               Added possibility of getting initial guess for flux density
 *                ratios from do_chi rather than just taking the ratio of
 *                the mean values.
 * v22Feb99 CDF, Changed output from function.  Now, in addition to printing
 *                out full 2D chisq surface as a function of (ratio, delay)
 *                pairs, also print out a constant-ratio cut through the
 *                surface to get chisq as a function of delay only at the
 *                best-fitting value for the delay.
 * v30Jul00 CDF, Replaced the ratio0, dtau, nmu, and ntau variables
 *                with the two new Prange structure variables mu and tau.
 * v07Sep00 CDF, Forced setting of the mu variable to be in the set_mu0
 *                function in lc_funcs.c so that the value of mu0 passed to
 *                this function will always be the desired value.
 * v08Sep00 CDF, Added a flag to the find_min_chisq function to turn on or
 *                off the parabola fitting functionality.
 * V11Sep00 CDF, Changed the bestlag, bestchi, chisq, and chiptr variable2
 *                 from Fluxrec to LCchisq structures.
 *               Moved chisq-slice printing to new function print_chisq_slice.
 */

int shift_chi(Fluxrec *control, Fluxrec *comp, int nctrl, int ncomp, 
	      float tau0, float mu0, Setup *setup, int speed, 
	      LCchisq *bestlag, char *outname, int doprint)
{
  int i;                 /* Looping variable */
  int no_error=1;        /* Flag set to 0 on error */
  float ratio;           /* Flux ratio between two curves */
  char slicename[MAXC];  /* Name for the slice output file */
  LCchisq *bestchi;      /* Container for data associated w/min chisq */
  LCchisq *chisq=NULL;   /* Curve of (day,ratio,chisq) group */
  Prange mu;             /* Parameter search values for mu grid */
  Prange tau;            /* Parameter search values for tau grid */
  FILE *fullfp=NULL;     /* Output file pointer for full chisq surface */

  /*
   * Initialize bestlag->rchisq to something ridiculous
   */

  bestlag->rchisq = -999.0;

  /*
   * Transfer proper information into the mu Prange variable from mu0
   *  and the setup structure.
   */

  printf(" shift_chi: ---------------------------------------------------");
  printf("---------------\n");
  mu.val0 = mu0;
  mu.nval = setup->nmu;
  mu.dval = FLUXSTEP;
  printf(" shift_chi: mu0=%7.4f, mu_min=%7.4f, mu_max=%7.4f, ",mu.val0,
	 mu.val0*(1.0-mu.nval*mu.dval),mu.val0*(1.0+mu.nval*mu.dval));
  printf("dmu=%6.4f, nmu=%d\n",mu.dval,2*mu.nval+1);

  /*
   * Find the lag axis spacing, central delay, and number of delay
   *  steps.  The central delay and number of steps can be set in the Setup 
   *  structure and passed as non-zero values to shift_chi.  If they are not, 
   *  then set them to be 0 and ncomp/2, respectively.
   */

  if(no_error) {
    tau.dval = (comp+1)->day - comp->day;
    if(setup->ntau == 0) {
      tau.val0 = 0;
      tau.nval = ncomp/2;
      tau.minstep = -tau.nval;
      tau.maxstep = tau.nval;
    }
    else {
      printf(" shift_chi: Using tau0 value contained in Setup structure.\n");
      tau.val0 = tau0;
      tau.nval = setup->ntau;
      tau.minstep = (int) (tau.val0/tau.dval) - tau.nval;
      tau.maxstep = (int) (tau.val0/tau.dval) + tau.nval;
    }
  }

  /*
   * Print out tau-search information
   */

  printf(" shift_chi: tau0=%6.1f, tau_min=%6.1f, tau_max=%6.1f, ",tau.val0,
	 tau.minstep*tau.dval,tau.maxstep*tau.dval);
  printf("dtau=%5.2f, ntau=%d\n",tau.dval,2*tau.nval+1);

  /*
   * Allocate memory for chisq array.  It will contain 2*tau.nval+1 members.
   */

  if(!(chisq = new_lcchisq(2*tau.nval+1)))
    no_error = 0;

  /*
   * Open the output files and write header line
   */

  if(doprint && no_error) {
    if(!(fullfp = fopen(outname,"w")))
      no_error =0;
    if(no_error) {
      fprintf(fullfp,"# Gridsize %d %d\n",2*tau.nval+1,2*mu.nval+1);
      fprintf(fullfp,"#\n");
      fprintf(fullfp,"# tau     mu     chi^2   N_DOF red. chi^2\n");
      fprintf(fullfp,"#------ ------ --------- ----- ----------\n");
    }
  }

  /*
   * Start the calculations.  Make the outer loop on flux ratios to get 
   *  data into proper Fortran order for making contour plots.
   */
   
  i = -mu.nval;
  while(i<mu.nval+1 && no_error) {

    /*
     * Set value of flux ratio for this iteration
     */

    ratio = mu.val0 * (1  + i*mu.dval);

    /*
     * Step through the delays at this ratio, putting reduced chisq
     *  values into the chisq array.
     */
    
    if(step_delays(control,comp,nctrl,ncomp,ratio,tau,chisq, 
		   setup,speed,fullfp,doprint))
      no_error = 0;

    /*
     * Find the minimum value of chisq
     */

    bestchi = find_min_chisq(chisq,2*tau.nval+1,0);
    if(bestchi->rchisq < 0.0)
      no_error = 0;

    /*
     * Check to see if we have a better chisq than the previous
     *  best chisq or if this is the first time around 
     *  (bestlag->err = -999.0).  If either one is true, put the
     *  current value of bestchi into bestlag.
     * 
     */

    if(bestlag->rchisq < 0 || bestchi->rchisq < bestlag->rchisq)
      *bestlag = *bestchi;

    i++;
  }

  /*
   * Print out slice through 2D chisq surface at the best-fit value of
   *  the flux density ratio.  First, create the slice by calling step_delays
   *  again, with the ratio set to the proper ratio (contained in 
   *  bestlag->flux).
   */

  if(no_error)
    if(step_delays(control,comp,nctrl,ncomp,bestlag->mu,tau,
		   chisq,setup,speed,NULL,0))
      no_error = 0;

  /*
   * Now print out the values contained in the chisq array.
   */

  if(doprint && no_error) {
    sprintf(slicename,"%s_slice",outname);
    printf(" ");
    if(print_chisq_slice(chisq,2*tau.nval+1,slicename))
      no_error = 0;
  }

  /*
   * Clean up
   */

  if(fullfp)
    fclose(fullfp);

  chisq = del_lcchisq(chisq);

  if(no_error)
    return 0;
  else {
    fprintf(stderr,"ERROR: shift_chi\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function step_delays
 *
 * Given two input curves and a magnification ratio between them, steps
 *  through a number of delays between them, computing a reduced chisq 
 *  statistic for each delay.  This functionality used to be contained
 *  in shift_chi.
 *
 * Inputs: Fluxrec *control    control light curve
 *         Fluxrec *comp       comparison light curve
 *         int nctrl           number of points in control curve
 *         int ncomp           number of points in comparison curve
 *         float ratio         flux ratio between the two curves
 *         Prange tau          delay parameter search values
 *         LCchisq *chisq      array of (day,ratio,chisq) structure
 *         Setup *setup        smoothing/interpolation info
 *         int speed           if the input curves are of the same size and
 *                              have the same spacing, this is set to FAST
 *                              and a quicker form of calculation can occur.
 *                              Otherwise, it will be set to SLOW.
 *         FILE *ofp           output file pointer
 *         int doprint         flag set to 1 if output files are desired
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 * v30Jul00 CDF, Modified the delay steps to run from the passed tau.minstep
 *                to tau.maxstep instead of running the full range from
 *                -ntau to +ntau.  This required the creation of a new
 *                Prange structure to reduce the number of passed parameters.
 * v08Sep00 CDF, Fixed bug in chisq calculation in which only one of the
 *                error quantities was used.
 *               Changed output format slightly.
 * v11Sep00 CDF, Changed chisq from a Fluxrec structure to a LCchisq
 *                structure.
 */

int step_delays(Fluxrec *control, Fluxrec *comp, int nctrl, int ncomp, 
		float ratio, Prange tau, LCchisq *chisq, 
		Setup *setup, int speed, FILE *ofp, int doprint)
{
  int i,j,m;             /* Looping variables */
  int pos;               /* Array position */
  int nmatch;            /* Number of overlapping points in light curves */
  float chidiff;         /* Difference in values used to calculate chisq */
  float sig12,sig22;     /* Squared error values used to calculate chisq */
  float sum;             /* Sum of individual chisqs */
  Fluxrec *fptr1;        /* Pointer to navigate control curve */
  Fluxrec *fptr2;        /* Pointer to navigate comparison curve */
  LCchisq *chiptr;       /* Pointer to navigate chisq array */


  /*
   * Step through lags going from -ntau to ntau
   * At each lag i, step through control (incrementing j) and compute
   *  the chisq using control[j] and comp[m] if where m is set such that
   *  control[j].day - lag = comp[m].day.
   */

  for(i=tau.minstep,chiptr=chisq; i<tau.maxstep+1; i++,chiptr++) {

    /*
     * Set the value of the day and ratio members of the chisq array
     */

    chiptr->tau = i*tau.dval;
    chiptr->mu = ratio;

    /*
     * Initialize the running counter and sum variables.
     */

    nmatch = 0;
    sum = 0.0;

    /*
     * Now step through the control curve and compute the chisq contribution
     *  from the control and comparison (with the given lag), as long as
     *  the control point is a measured value (match = 1).  We can use the
     *  FAST version if the control and comparison curves have the same
     *  size and spacing.  Otherwise, we have to do the more tedious
     *  SLOW calculations.
     */

    if(speed == FAST) {

      /*
       * Fast case
       */
	
      for(j=0,fptr1=control,fptr2=comp; j<nctrl; j++,fptr1++,fptr2++) {
	pos = j + i;
	if(fptr1->match >= 0 && pos >= 0 && pos < ncomp) {
	  chidiff = fptr1->flux - (fptr2+i)->flux / ratio;
	  sig12 = fptr1->err * fptr1->err;
	  sig22 = (fptr2+i)->err * (fptr2+i)->err / (ratio * ratio);
	  sum += chidiff * chidiff / (sig12 + sig22);
	  nmatch++;
	}
      }
    }
    else {

      /*
       * Slow case (the default)
       */

      for(j=0,fptr1=control; j<nctrl; j++,fptr1++) {
	if(fptr1->match == 1) {
	  for(m=0,fptr2=comp; m<ncomp; m++,fptr2++) {
	    if(fabs(fptr1->day + chiptr->tau - fptr2->day) < 0.2*tau.dval) {
	      chidiff = fptr1->flux - fptr2->flux / ratio;
	      sig12 = fptr1->err * fptr1->err;
	      sig22 = fptr2->err * fptr2->err / (ratio * ratio);
	      sum += chidiff * chidiff / (sig12 + sig22);
	      nmatch++;
	      break;
	    }
	  }
	}
      }
    }

    /*
     * If both curves are interpolated (speed == FAST), compute the number
     *   of independent points in the overlap region by dividing the
     *   time of overlap (nmatch * tau.dval) by the average spacing of points
     *   in the raw curves.
     */

    if(speed == FAST) {
#if 0
      switch(setup->smtype) {
      case SMBOXCAR: case SMTRIANGLE:
	nmatch = floor(nmatch * tau.dval * WINFRAC / setup->smwidth);
	break;
      default:
      }
#endif
      nmatch = floor(nmatch * tau.dval / 3.7);
    }

    /*
     * Compute reduced chisq
     */

    if(nmatch <= NFIT) {
      fprintf(stderr,
	      "step_delays: Warning! Overlap between curves is too small.\n");
      chiptr->rchisq = -999.0;
      chiptr->ndof  = -99;
    }
    else {
      chiptr->chisq  = sum;
      chiptr->ndof   = nmatch - NFIT;
      chiptr->rchisq = sum/(nmatch-NFIT);
    }

    /*
     * Print out lag and reduced chisq
     */

    if(doprint) {
      fprintf(ofp,"%7.2f %6.4f %9.2f %5d  %8.3f\n",chiptr->tau,chiptr->mu,
	      chiptr->chisq,chiptr->ndof,chiptr->rchisq);
    }
  }

  return 0;
}

/*.......................................................................
 *
 * Function find_min_chisq
 *
 * Takes a Fluxrec array containing points of the nature
 *     fluxrec->day = day
 *     fluxrec->err = chisq
 *  and finds the array member with the lowest value of chisq.  Then
 *  fits a parabola (with the function fit_parab) to the points surrounding
 *  the minimum-chisq array member to get a "true" minimum chisq.
 *
 * Inputs: Fluxrec *chisq      chisq array
 *         int npoints         number of points in the array
 *         int fitparabola     flag set to 1 for parabola_fitting
 *
 * Output: Fluxrec *bestchi    day giving "true" minimum chisq
 *
 * v08Sep00 CDF, Added the fitparabola flag to toggle whether the
 *                fit_parab function is called.
 * v11Sep00 CDF, Changed output of function to be a LCchisq structure.
 *               Also changed chisq, bestchi, and chiptr.
 *
 */

LCchisq *find_min_chisq(LCchisq *chisq, int npoints, int fitparabola)
{
  int i;             /* Looping variable */
  int no_error=1;    /* Flag set to 0 on error */
  float a,b,c;       /* Parameters of the parabolic fit */
  LCchisq *bestchi;  /* "True" min chisq */
  LCchisq *chiptr;   /* Pointer to navigate chisq */

  /*
   * Check to make sure there are at least 3 grid points in chisq
   */

  if(npoints < 3) {
    fprintf(stderr,"ERROR: find_min_chisq.  ");
    fprintf(stderr,"There must be at least 3 grid points to find min chisq.\n");
    bestchi->rchisq = -999.0;
    return bestchi;
  }

  /*
   * Initialize bestchi
   */

  bestchi = chisq;

  /*
   * Find minimum chisq on grid (must be > 0)
   */

  for(i=0,chiptr=chisq; i<npoints; i++,chiptr++) {
    if(chiptr->rchisq < bestchi->rchisq && chiptr->rchisq >= 0.0)
      bestchi = chiptr;
  }

  /*
   * If the fitparabola flag is set to 1, call fit_parab to fit a parabola 
   *  to lowest grid point + 2 neighboring points.  There are three possible 
   *  cases:
   *   1. bestchi is at the first grid point ==> just take bestchi
   *   2. bestchi is at the last grid point ==> just take bestchi
   *   3. bestchi is any other grid point ==> take bestchi plus the points
   *      on either side of it.
   */

  if(fitparabola) {
    if(bestchi != chisq && bestchi != chisq+npoints-1) {
      if(fit_parab((bestchi-1)->tau,(bestchi-1)->rchisq,bestchi->tau,
		   bestchi->rchisq,(bestchi+1)->tau,(bestchi+1)->rchisq,
		   &a,&b,&c,0))
	no_error = 0;

      if(no_error) {
	bestchi->tau = -b / (2*a);
	bestchi->rchisq = a * bestchi->tau * bestchi->tau + 
	  b * bestchi->tau + c;
      }
      else {
	bestchi->rchisq = -999.0;
      }
    }
  }

  if(bestchi->rchisq < 0.0)
    fprintf(stderr,"ERROR: find_min_chisq.  Best chisq is less than 0!\n");
  return bestchi;
}

/*.......................................................................
 *
 * Function print_chisq_slice
 *
 * Prints out the reduced chisq as a function of tau.  This is a slice
 *  from the (tau, mu) plane, hence the name of the function.
 *
 * Inputs: LCchisq *chisq      structure containing chisq information
 *         int size            size of array
 *         char *outname       output file name
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 */

int print_chisq_slice(LCchisq *chisq, int size, char *outname)
{
  int i;             /* Looping variable */
  LCchisq *chiptr;   /* Pointer to navigate chisq */
  FILE *ofp=NULL;    /* Output file pointer */

  /*
   * Open output file
   */

  if(!(ofp = open_writefile(outname))) {
    fprintf(stderr,"ERROR: print_chisq_slice.\n");
    return 1;
  }

  /*
   * Print header
   */

  fprintf(ofp,"# tau     mu     chi^2   N_DOF red. chi^2\n");
  fprintf(ofp,"#------ ------ --------- ----- ----------\n");

  /*
   * Print data
   */

  for(i=0,chiptr=chisq; i<size; i++,chiptr++)
    fprintf(ofp,"%7.2f %6.4f %9.2f %5d  %8.3f\n",chiptr->tau,chiptr->mu,
	    chiptr->chisq,chiptr->ndof,chiptr->rchisq);

  /*
   * Clean up and exit
   */

  if(ofp)
    fclose(ofp);
  return 0;
}
