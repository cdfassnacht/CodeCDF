/*
 * noninterp_fns.c
 *
 * A library of functions to find time delays without interpolating
 *  the 1608 light curves.
 *
 * 30Nov98 CDF,  Created by moving some DCF functions from correlate.c
 * v02Dec98 CDF, Added call_disp, disp_setup, and disp_lovell to do
 *                dispersion analysis as described in the Lovell et al.
 *                1830 time delay paper.
 * v03Dec98 CDF, Added disp_d1 and disp_d2 to calculate D^2_1 and D^2_2
 *                functions as described in the Pelt et al. papers.
 * v04Dec98 CDF, Moved delta parameter for the D^2_2 analysis into the
 *                Setup structure.
 * v07Dec98 CDF, Rearranged things to make disp_d1 and disp_d2 more general.
 * v14Jan99 CDF, Moved a lot of the function of disp_setup into two_curve_disp,
 *                which gets called by disp_setup.  Then copied and modified
 *                two_curve_disp to become four_curve_disp to deal with the
 *                dispersion methods which treat all four light curves at once.
 * v09Oct00 CDF, Modified disp_setup and two_curve_disp to use Prange
 *                structures to define the search ranges rather than the
 *                previously used arrays.
 *               Added new_lcdisp, del_lcdisp, and read_lcdisp functions
 *                to handle the new LCdisp structure memory allocation and
 *                file input.  These are direct modifications of the similar
 *                function to handle the LCchisq structure in lc_chisq.c.
 * v10Oct00 CDF, Modified call_disp to pass a LCdisp structure array to
 *                disp_setup and hence to the dispersion-calculating 
 *                functions.  This allows the grid-points producing the
 *                lowest dispersion(s) to be passed back up to call_disp.
 *               Added printing out of each dispersion calculated.
 *               Added print_disp_slice to print out a (tau,disp) curve.
 * v19Apr01 CDF, Moved calculation of composite curve for D^2_1 and D^2_2
 *                methods from the disp_d1 and disp_d2 functions up into
 *                the calling functions (two_curve_disp or four_curve_disp).
 *                This makes the disp_d1 and disp_d2 functions more general.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "structdef.h"
#include "dataio.h"
#include "lc_funcs.h"
#include "noninterp_fns.h"

#define DAYSTEP 0.5

/*.......................................................................
 *
 * Function new_lcdisp
 *
 * Allocates dynamic memory for a pointer array of Lcdisp structures
 *
 * Input:  int size            size of the array
 *
 * Output: Lcdisp *newstruct   pointer to the new array.  NULL if error
 *
 */

LCdisp *new_lcdisp(int size)
{
  int i;
  LCdisp *newstruct;
  LCdisp *lcptr;
 
  newstruct = (LCdisp *) malloc(sizeof(LCdisp) * size);
  if(!newstruct) {
    fprintf(stderr,"new_lcdisp:  Insufficient memory for array.\n");
    return NULL;
  }

  /*
   * Initialize values in the array
   */

  for(i=0,lcptr=newstruct; i<size; i++,lcptr++) {
    lcptr->tau = lcptr->mu = lcptr->disp = 0.0;
    lcptr->arraypos = 0;
  }

  return newstruct;
}

/*.......................................................................
 *
 * Function del_lcdisp
 *
 * Frees up memory allocated to Lcdisp array
 *
 * Input:  Lcdisp *lcdisp      array to be freed
 *
 * Output: NULL
 */

LCdisp *del_lcdisp(LCdisp *lcdisp)
{
  if(lcdisp)
    free(lcdisp);
 
  return NULL;
}

/*.......................................................................
 *
 * Function read_lcdisp
 *
 * Reads a grid of Lcdisp structures from one of the disp*dat files
 *  produced by delays.c.  The filename is passed as one of the
 *  arguments.
 *
 * Inputs: char *filename      input filename
 *         Axisinfo *tau       information about delay values
 *         Axisinfo *mu        information about magnification values
 *
 * Output: Lcdisp *lcdisp      filled Lcdisp structure.  NULL on error.
 *
 */

LCdisp *read_lcdisp(char *filename, Axisinfo *tau, Axisinfo *mu)
{
  int no_error=1;         /* Flag set to 0 on error */
  int nlines;             /* Number of data lines in input file */
  char line[MAXC];        /* General string variable for getting input */
  LCdisp *lcdisp=NULL;  /* Lcdisp array to be filled */
  LCdisp *lcptr;         /* Pointer to navigate lcdisp */
  FILE *ifp=NULL;         /* Input file pointer */

  /*
   * Open input file
   */

  if(!(ifp = open_readfile(filename))) {
    fprintf(stderr,"ERROR: read_lcdisp\n");
    return NULL;
  }

  /*
   * Count number of lines in input file
   */

  if((nlines = n_lines(ifp,'#')) == 0) {
    fprintf(stderr,"ERROR: read_lcdisp. No data in input file %s.\n",
	    filename);
    no_error = 0;
  }
  else {
    printf("read_lcdisp: %s contains %d data lines.\n",filename,nlines);
    rewind(ifp);
  }

  /*
   * Allocate memory for Lcdisp array
   */

  if(!(lcdisp = new_lcdisp(nlines)))
    no_error = 0;

  /*
   * Read in header information, which should be contained on the
   *  first line of the input file.
   */

  if(no_error) {
    fgets(line,MAXC,ifp);
    if(sscanf(line,"# Gridsize %d %d",&tau->nval,&mu->nval) != 2) {
      fprintf(stderr,"ERROR: read_lcdisp.  %s missing header information.\n",
	      filename);
      no_error = 0;
    }
    else
      printf("read_lcdisp: Input array from %s has dimensions %d x %d\n",
	     filename,tau->nval,mu->nval);
  }

  /*
   * Read in data
   */

  if(no_error) {
    tau->minval = mu->minval = 9999.9;
    tau->maxval = mu->maxval = -9999.9;
    lcptr = lcdisp;
    while(fgets(line,MAXC,ifp) != NULL && no_error) {
      if(line[0] != '#') {
	if(sscanf(line,"%f %f %f",&lcptr->tau,&lcptr->mu,&lcptr->disp) != 3) {
	  fprintf(stderr,"ERROR: read_lcdisp. Bad input format.\n");
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
    printf("read_lcdisp: tau_min = %7.2f, tau_max = %7.2f, ntau = %d\n",
	   tau->minval,tau->maxval,tau->nval);
    printf("read_lcdisp: mu_min  =  %6.4f, mu_max  =  %6.4f, nmu  = %d\n",
	   mu->minval,mu->maxval,mu->nval);
  }

  /*
   * Clean up and exit
   */

  if(ifp)
    fclose(ifp);

  if(no_error)
    return lcdisp;
  else {
    fprintf(stderr,"ERROR: read_lcdisp.\n");
    return del_lcdisp(lcdisp);
  }
}

/*.......................................................................
 *
 * Function call_disp
 *
 * Calls one of the dispersion functions (via disp_setup) for all pairs of 
 *  light curves
 *
 * Inputs: Fluxrec *flux[]     container for the 1608 light curves
 *         int npoints         number of points in each light curve
 *         Setup *setup        container for dispersion method info
 *         int doprint         flag set to 1 to produce output files
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 * v10Oct00 CDF, Added bestdisp container for gridpoints producing
 *                lowest dispersions.
 *               Moved printing out of lowest-dispersion points from
 *                two_curve_disp to this function.
 * v31Aug02 CDF, Got rid of nbad array a passed parameter (no longer
 *                needed for make_compos).  Introduced the nptarr
 *                array for new version of make_compos.
 */

int call_disp(Fluxrec *flux[], int npoints, Setup *setup, int doprint)
{
  int i;                 /* Looping variable */
  int no_error=1;        /* Flag set to 0 on error */
  int nptarr[N08];       /* Array version of number of pts in each curve */
  int index[N08];        /* Array showing which curves are being compared */
  LCdisp *bestdisp=NULL; /* Container for gridpoint(s) giving lowest disp. */
  LCdisp *dptr;          /* Pointer for navigating bestdisp */
  FILE *bfp=NULL;        /* File pointer for lowest-dispersion results */

  /*
   * Initialize index
   */

  for(i=0; i<N08; i++)
    index[i] = 0;

  /*
   * Convert number of points into an array
   */

  for(i=0; i<N08; i++)
    nptarr[i] = npoints;

  /*
   * Allocate memory for bestdisp array.
   */

  if(!(bestdisp = new_lcdisp(N08-1))) {
    fprintf(stderr,"ERROR: call_disp.\n");
    return 1;
  }

  /*
   * Call the dispersion function
   */

  switch(setup->dispchoice) {

  /*
   * Dispersion methods for pairs of curves (D21, D22, DLOVELL)
   */

  case D21: case D22: case DLOVELL:

    /*
     * Print out informative headers
     */

    printf("\nCalculating dispersion spectra (pairs of curves)...\n\n");

    /*
     * B and A
     */

    index[0] = 1;
    index[1] = 0;
    if(disp_setup(flux,2,nptarr,index,setup,bestdisp,"dispba.dat",
		  doprint))
      no_error = 0;

    /*
     * B and C
     */

    index[1] = 2;
    if(no_error)
      if(disp_setup(flux,2,nptarr,index,setup,bestdisp+1,"dispbc.dat",
		    doprint))
	no_error = 0;

    /*
     * B and D
     */

    index[1] = 3;
    if(no_error)
      if(disp_setup(flux,2,nptarr,index,setup,bestdisp+2,"dispbd.dat",
		    doprint))
	no_error = 0;


    break;

  /*
   * Dispersion methods for all curves at once (D21M)
   */

  case D21M:

    /*
     * Print out informative headers
     */

    printf("\nCalculating dispersion spectra (all curves at once)...\n\n");

    /*
     * Set up indexing
     */

    index[0] = 1;
    index[1] = 0;
    index[2] = 2;
    index[3] = 3;

    /*
     * Call disp_setup
     */
    if(no_error)
      if(disp_setup(flux,4,nptarr,index,setup,bestdisp,"dispall.dat",0))
	no_error = 0;

    break;

  default:
    fprintf(stderr,"ERROR: call_disp.  Invalid dispersion method.\n");
    no_error = 0;
  }

  /*
   * Print out results
   */

  if(no_error) {
    printf(" disp_setup:----------------------------------------------------");
    printf("---------------\n\n");
    printf("call_disp: Best fit results:\n");
    printf("         Pair     tau         mu        Disp.\n");
    printf("         -----  -------   ----------   --------\n");
    dptr = bestdisp;
    printf("          B-A  %7.2f      %6.4f     %8.4f\n",dptr->tau,
	   dptr->mu,dptr->disp);
    dptr = bestdisp + 1;
    printf("          B-C  %7.2f      %6.4f     %8.4f\n",dptr->tau,
	   dptr->mu,dptr->disp);
    dptr = bestdisp + 2;
    printf("          B-D  %7.2f      %6.4f     %8.4f\n",dptr->tau,
	   dptr->mu,dptr->disp);

    if(!(bfp = open_appendfile("disp.log",0)))
      no_error = 0;
    else {
      for(i=0,dptr=bestdisp; i<N08-1; i++,dptr++)
	fprintf(bfp,"%7.2f %6.4f %8.4f ",dptr->tau,dptr->mu,dptr->disp);
      fprintf(bfp,"\n");
    }
  }

  /*
   * Clean up and exit
   */

  if(bfp)
    fclose(bfp);
  if(no_error)
    return 0;
  else {
    fprintf(stderr,"ERROR: call_disp\n");
    return 1;
  }
}


/*.......................................................................
 *
 * Function disp_setup
 *
 * Sets up the proper variables for calls to dispersion
 *
 * Inputs: Fluxrec *flux[]     input light curves
 *         int ncurves         number of curves being compared
 *         int *npoints        number of points in each light curve
 *         int *index          array showing which curves are being compared
 *         Setup *setup        container for dispersion method info
 *         char *outname       name of output file
 *         int doprint         flag set to 0 for no output file
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 * v07Dec98 CDF, Moved conversion of tau and mu to arrays into this
 *                function in order to make disp_d1 and disp_d2 more general.
 * v09Oct00 CDF, Changed mu0 and tau0 arrays into Prange structures, which
 *                contain more info in a more compact form.
 * v10Oct00 CDF, Added a passed variable *bestdisp to bring the grid
 *                point(s) giving the lowest dispersion back to call_disp.
 * v31Aug02 CDF, Changed npoints to an array.  Got rid of nbad array
 *                passed parameter (no longer needed for make_compos).
 */

int disp_setup(Fluxrec *flux[], int ncurves, int *npoints, 
	       int *index, Setup *setup, LCdisp *bestdisp, char *outname, 
	       int doprint)
{
  int i;                  /* Looping variables */
  int no_error=1;         /* Flag set to 0 on error */
  Prange *mu0=NULL;       /* Parameter search values for mu grid(s) */
  Prange *tau0=NULL;      /* Parameter search values for tau grid(s) */
  Prange *pptr;           /* Pointer to navigate Prange structures */

  /*
   * Allocate memory for Prange arrays.  The size will be one less than
   *  the number of curves.
   */

  if(!(mu0 = new_prange(ncurves-1)))
    no_error = 0;
  else
    if(!(tau0 = new_prange(ncurves-1)))
      no_error = 0;

  /*
   * Put the proper mu information into the Prange container.  This
   *  information is stored in the Setup container, and the assignment
   *  is controlled by the index array.  Note that index[0] will always
   *  contain the index of the UNSHIFTED/UNMAGNIFIED curve.  This means
   *  that we only load the values for index[1] -> index[ncurves-1].
   */

  if(no_error) {
    printf(" disp_setup:----------------------------------------------------");
    printf("---------------\n");
    for(i=0,pptr=mu0; i<ncurves-1; i++,pptr++) {
      pptr->val0 = setup->mu0[index[i+1]];
      pptr->nval = setup->nmu;
      pptr->dval = FLUXSTEP;
      printf(" disp_setup: mu0=%7.4f, mu_min=%7.4f, mu_max=%7.4f, ",
	     pptr->val0,pptr->val0*(1.0 - pptr->nval * pptr->dval),
	     pptr->val0*(1.0 + pptr->nval * pptr->dval));
      printf("dmu=%6.4f, nmu=%d\n",pptr->dval * pptr->val0,2 * pptr->nval + 1);

    }
  }

  /*
   * Calculate the number of time delays to consider if this value has
   *  not already been set in the Setup function.
   * For a spacing of one day, this is just half the number of days spanned
   *  by the observations, since we will consider delays from
   *  -ndays/2 to ndays/2
   */

  if(no_error) {
    if(setup->ntau == 0)
      setup->ntau = ((int) (flux[index[0]]+npoints[index[0]]-1)->day - 
		     flux[index[0]]->day);
  }

  /*
   * Set tau grid spacing to default value if hasn't been already set.
   */

  if(no_error && setup->dtau == 0.0)
    setup->dtau = DAYSTEP;

  /*
   * Put the tau information into the Prange container.  This
   *  information is stored in the Setup container, and the assignment
   *  is controlled by the index array.  Note that index[0] will always
   *  contain the index of the UNSHIFTED/UNMAGNIFIED curve.  This means
   *  that we only load the values for index[1] -> index[ncurves-1].
   */

  if(no_error) {
    for(i=0,pptr=tau0; i<ncurves-1; i++,pptr++) {
      pptr->val0 = setup->tau0[index[i+1]];
      pptr->nval = setup->ntau;
      pptr->dval = setup->dtau;
      pptr->minstep = ((int) (pptr->val0/pptr->dval)) - pptr->nval;
      pptr->maxstep = ((int) (pptr->val0/pptr->dval)) + pptr->nval;
      printf(" disp_setup: tau0=%6.1f, tau_min=%6.1f, tau_max=%6.1f, ",
	     pptr->val0,pptr->minstep * pptr->dval,
	     pptr->maxstep * pptr->dval);
      printf("dtau=%5.2f, ntau=%d\n",pptr->dval,2 * pptr->nval + 1);
    }
  }

  /*
   * Call the two-curve or four-curve function, depending on the
   *  value of setup->dispchoice
   */

  if(no_error) {
    switch(setup->dispchoice) {
    case D21: case D22: case DLOVELL:
      if(two_curve_disp(flux,npoints,index,tau0,mu0,setup,bestdisp,
			outname,doprint))
	no_error = 0;
      break;
    case D21M:
      printf(" disp_setup:-----------------------------------------------");
      printf("--------------------\n");
      if(four_curve_disp(flux,npoints,index,tau0,mu0,setup,bestdisp,
			 outname,0))
	no_error = 0;
      break;
    default:
      fprintf(stderr,"ERROR: disp_setup: Unsupported dispersion method.\n");
      no_error = 0;
    }
  }

  if(no_error)
    printf(" disp_setup: bestdisp(%d): %6.1f %8.5f %f\n",index[1],
	   bestdisp->tau,bestdisp->mu,bestdisp->disp);

  /*
   * Clean up and exit
   */

  mu0 = del_prange(mu0);
  tau0 = del_prange(tau0);

  if(no_error)
    return 0;
  else {
    fprintf(stderr,"ERROR: disp_setup\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function two_curve_disp
 *
 * Sets things up for calling one of the dispersion routines for calculating
 *  dispersions between two curves only (as opposed to multiple curves, 
 *  which are treated in four_curve_disp).
 *
 * Inputs: Fluxrec *flux[]     input light curves
 *         int *npoints        number of points in each light curve
 *         int *index          array showing which curves are being compared
 *         int ntau            maximum time delay (relative to initial guess) 
 *                              considered
 *         Prange *tau0        parameters for tau grid search
 *         Prange *mu0         parameters for mu grid search
 *         Setup *setup        container for dispersion method info
 *         char *outname       name of output file
 *         int doprint         flag set to 0 for no output file
 *
 * Output: LCdisp *d2min       (mu,tau) pair(s) with lowest dispersion(s)
 *
 * v11Feb99 CDF, Moved saving of best-fit spectrum to outside of mu loop
 *                to avoid unnecessary repetition.
 * v09Oct00 CDF, Modified to take advantage of passing Prange structures
 *                for mu0 and tau0 rather than float arrays.
 * v10Oct00 CDF, Added a passed variable *bestdisp to bring the grid
 *                point(s) giving the lowest dispersion back to disp_setup.
 *               Moved printing of lowest-dispersion point into call_disp
 *                function.
 *               Moved printing of slice through plane to new 
 *                print_disp_slice function.
 * v19Apr01 CDF, Moved creation of composite curve up to this function
 *                from the disp_d1 and disp_d2 functions.
 *               Also moved printing of output info up from disp_d1 and
 *                disp_d2 into this function.
 * v31Aug02 CDF, Changed npoints to an array.  Got rid of nbad array
 *                a passed parameter (no longer needed for make_compos).
 */

int two_curve_disp(Fluxrec *flux[], int *npoints, int *index,
		   Prange *tau0, Prange *mu0, Setup *setup, 
		   LCdisp *bestdisp, char *outname, int doprint)
{
  int i,j;                  /* Looping variables */
  int no_error=1;           /* Flag set to 0 on error */
  int ncurves=2;            /* Number of curves being combined */
  int ncompos=0;            /* Number of points in the composite curve */
  float tau[N08];           /* Time delays between the curves */
  float mu[N08];            /* Flux density ratios between the curves */
  char slicename[MAXC];     /* Name of dispersion spectrum file */
  Fluxrec *compos=NULL;     /* Composite curve */
  LCdisp d2minj;            /* Min value of D^2 in one loop */
  LCdisp *d2min=NULL;       /* Absolute min value of D^2 */
  LCdisp *d2arr=NULL;       /* Array containing dispersion values */
  LCdisp *dptr;             /* Pointer to navigate d2arr */
  FILE *ofp=NULL;           /* Output file pointer */

  /*
   * Initialize mu and tau arrays.
   */

  for(i=0; i<N08; i++) {
      tau[i] = 0.0;
      mu[i] = 1.0;
  }

  /*
   * Allocate memory for minimum dispersion container.
   */

  if(!(d2min = new_lcdisp(ncurves-1))) {
    fprintf(stderr,"ERROR: two_curve_disp.\n");
    return 1;
  }

  /*
   * Open output file.
   */

  if(doprint && no_error) {
    printf(" ");
    if(!(ofp = open_writefile(outname))) {
      no_error = 0;
      fprintf(stderr,"ERROR: two_curve_disp\n");
    }
    else {
      fprintf(ofp,"# Gridsize %d %d\n",2 * tau0->nval + 1,2 * mu0->nval + 1);
      fprintf(ofp,"#\n");
      fprintf(ofp,"# tau     mu     disp   \n");
      fprintf(ofp,"#------ ------ ---------\n");
    }
  }

  /*
   * Allocate memory for a slice through the dispersion surface at a
   *  constant value of mu
   */

  if(!(d2arr = new_lcdisp(2 * tau0->nval + 1)))
    no_error = 0;

  /*
   * Outer loop on flux density ratio
   */

  if(no_error) {
    for(i=-mu0->nval; i < mu0->nval+1; i++) {

      /*
       * Set mu for this iteration
       */

      mu[index[1]] = mu0->val0 * (1.0 + i*mu0->dval);

      /*
       * Inner loop on tau
       */

      for(j=tau0->minstep,dptr=d2arr; j < tau0->maxstep+1; j++,dptr++) {
	tau[index[1]] = j * tau0->dval;
	dptr->tau = tau[index[1]];
	dptr->mu = mu[index[1]];

	/*
	 * For the D21 and D22 methods, create the composite curve (through 
	 *  a call to make_compos) before calling the function that calculates
	 *  the dispersion.
	 */

	switch(setup->dispchoice) {
	case D21:
	  if(!(compos = make_compos(flux,ncurves,npoints,index,tau,mu,
				    &ncompos,0)))
	    no_error = 0;
	  dptr->disp = disp_d1(compos,ncompos,setup->d2delta);
	  compos = del_fluxrec(compos);
	  break;
	case D22:
	  if(!(compos = make_compos(flux,ncurves,npoints,index,tau,mu,
				    &ncompos,0)))
	    no_error = 0;
	  dptr->disp = disp_d2(compos,ncompos,setup->d2delta);
	  compos = del_fluxrec(compos);
	  break;
	case DLOVELL:
	  dptr->disp = disp_lovell(flux[index[0]],flux[index[1]],npoints[0],
				   tau[index[1]],mu[index[1]],setup->d2delta);
	  break;
	default:
	  fprintf(stderr,"ERROR: two_curve_disp. Invalid dispersion method.\n");
	  fprintf(stderr," Using D^2_1 method.\n");
	  if(!(compos = make_compos(flux,ncurves,npoints,index,tau,mu,
				    &ncompos,0)))
	    no_error = 0;
	  dptr->disp = disp_d1(compos,ncompos,setup->d2delta);
	  compos = del_fluxrec(compos);
	}

	/*
	 * Print output if desired
	 */

	if(doprint)
	  fprintf(ofp,"%7.2f %6.4f %7.4f\n",tau[index[1]],mu[index[1]],
		  dptr->disp);

	/*
	 * Hold value giving minimum dispersion FOR THIS LOOP
	 */

	if(i == -mu0->nval && j == tau0->minstep)
	  d2minj = *dptr;
	else if(dptr->disp < d2minj.disp)
	  d2minj = *dptr;
      }

      /*
       * Save values giving ABSOLUTE minimum dispersion
       */

      if(i == -mu0->nval) {
	*d2min = d2minj;
      }
      else if(d2minj.disp < d2min->disp) {
	*d2min = d2minj;
      }
    }
  }

  /*
   * Make another call to the dispersion-calculating function with
   *  mu set to its best-fit value.  This will produce the dispersion
   *  spectrum with the best-fitting parameters.  Print out the results
   *  to the output file.
   */

  if(no_error && doprint) {

    /*
     * Set output filename
     */

    sprintf(slicename,"%s_slice",outname);

    /*
     * Now loop through the delays with mu set to d2min->mu, calculating
     *  the dispersion at each step.
     */


    for(j=tau0->minstep,dptr=d2arr; j < tau0->maxstep + 1; j++,dptr++) {
      tau[index[1]] = j*tau0->dval;
      mu[index[1]] = d2min->mu;

      /*
       * For the D21 and D22 methods, create the composite curve (through 
       *  a call to make_compos) before calling the function that calculates
       *  the dispersion.
       */

      switch(setup->dispchoice) {
      case D21:
	if(!(compos = make_compos(flux,ncurves,npoints,index,tau,mu,
				  &ncompos,0)))
	  no_error = 0;
  	dptr->disp = disp_d1(compos,ncompos,setup->d2delta);
	compos = del_fluxrec(compos);
	break;
      case D22:
	if(!(compos = make_compos(flux,ncurves,npoints,index,tau,mu,
				  &ncompos,0)))
	  no_error = 0;
	dptr->disp = disp_d2(compos,ncompos,setup->d2delta);
	compos = del_fluxrec(compos);
	break;
      case DLOVELL:
	dptr->disp = disp_lovell(flux[index[0]],flux[index[1]],npoints[0],
				 tau[index[1]],mu[index[1]],setup->d2delta);
	break;
      default:
	fprintf(stderr,"ERROR: two_curve_disp. Invalid dispersion method.\n");
	fprintf(stderr," Using D^2_1 method.\n");
	if(!(compos = make_compos(flux,ncurves,npoints,index,tau,mu,
				  &ncompos,0)))
	  no_error = 0;
  	dptr->disp = disp_d1(compos,ncompos,setup->d2delta);
	compos = del_fluxrec(compos);
      }
    }

    /*
     * Print values to output file
     */

    if(print_disp_slice(d2arr,2 * tau0->nval + 1,slicename))
      no_error = 0;
  }


  /*
   * Transfer lowest-dispersion information to bestdisp.
   */

  if(no_error) {
    for(i=0,dptr=d2min; i<ncurves-1; i++,dptr++)
      *(bestdisp+i) = *dptr;
  }

  /*
   * Clean up and exit
   */

  d2min = del_lcdisp(d2min);
  d2arr = del_lcdisp(d2arr);
  if(ofp)
    fclose(ofp);

  if(no_error)
    return 0;
  else {
    fprintf(stderr,"ERROR: two_curve_disp\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function four_curve_disp
 *
 * Sets things up for calling one of the dispersion routines for calculating
 *  dispersions between all four curves at once.
 *
 * Inputs: Fluxrec *flux[]     input light curves
 *         int *npoints        number of points in each light curve
 *         int *index          array showing which curves are being compared
 *         int ntau            maximum time delay (relative to initial guess) 
 *                              considered
 *         Prange *tau0        parameters for tau grid search
 *         Prange *mu0         parameters for mu grid search
 *         Setup *setup        container for dispersion method info
 *         char *outname       name of output file
 *         int doprint         flag set to 0 for no output file
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 * v09Oct00 CDF, Modified to take advantage of passing Prange structures
 *                for mu0 and tau0 rather than float arrays.
 * v19Apr01 CDF, Moved creation of composite curve up to this function
 *                from the disp_d1.
 * v31Aug02 CDF, Changed npoints to an array.  Got rid of nbad array
 *                a passed parameter (no longer needed for make_compos).
 */

int four_curve_disp(Fluxrec *flux[], int *npoints, int *index,
		    Prange *tau0, Prange *mu0, Setup *setup, 
		    LCdisp *bestdisp, char *outname, int doprint)
{
  int i,j,k,m,n,p;             /* Looping variables */
  int no_error=1;              /* Flag set to 0 on error */
  int count=0;                 /* Number of curve combinations tried */
  int ncurves=4;               /* Number of curves being combined */
  int ncompos=0;            /* Number of points in the composite curve */
  float tau[N08];              /* Time delays between the curves */
  float mu[N08];               /* Flux density ratios between the curves */
  Fluxrec *compos=NULL;     /* Composite curve */
  LCdisp *d2min=NULL;          /* Absolute min value of D^2 */
  LCdisp *d2arr[N08]={NULL};   /* Array containing dispersion values */
  LCdisp *minarr[N08]={NULL};  /* Dispersion spectrum containing min D^2 */
  LCdisp *dptr1;               /* Pointer to navigate d2arr[1] */
  LCdisp *dptr2;               /* Pointer to navigate d2arr[2] */
  LCdisp *dptr3;               /* Pointer to navigate d2arr[3] */
  LCdisp *mptr1;               /* Pointer to navigate minarr[1] */
  LCdisp *mptr2;               /* Pointer to navigate minarr[2] */
  LCdisp *mptr3;               /* Pointer to navigate minarr[3] */
  Prange *pmp1,*pmp2;          /* Pointers for mu0 array */
  Prange *ptp1,*ptp2;          /* Pointers for tau0 array */
  FILE *ofp=NULL;              /* Output file pointer */
  FILE *tmpfp=NULL;            /* Container for temporary output */

  /*
   * Initialize mu and tau arrays.
   */

  for(i=0; i<N08; i++) {
      tau[i] = 0.0;
      mu[i] = 1.0;
  }

  /*
   * Allocate memory for minimum dispersion container.
   */

  if(!(d2min = new_lcdisp(ncurves-1))) {
    fprintf(stderr,"ERROR: four_curve_disp.\n");
    return 1;
  }

  /*
   * Allocate memory for two-dimensional slices through the dispersion 
   *  hypersurface at constant values of mu
   */

  for(i=0; i<N08; i++) {
    if(!(d2arr[i] = new_lcdisp(2 * setup->ntau + 1)))
      no_error = 0;
    if(!(minarr[i] = new_lcdisp(2 * setup->ntau + 1)))
      no_error = 0;
  }

  /*
   * Start nested loops.  For each of the three comparison curves, have 
   *  an outer loop on flux ratio and an inner loop on time step.
   */

  
  for(i=-mu0->nval; i < mu0->nval + 1; i++) {

    mu[index[1]] = mu0->val0 * (1.0 + i * mu0->dval);

    /*
     * Loop on tau for first curve
     */

    for(j=tau0->minstep,dptr1=d2arr[1]; j < tau0->maxstep +1; j++,dptr1++) {
      tau[index[1]] = j*tau0->dval;
      dptr1->tau = tau[index[1]];
      dptr1->mu = mu[index[1]];

      /*
       * Loop on mu for second curve
       */

      pmp1 = mu0 + 1;
      for(k=-pmp1->nval; k < pmp1->nval + 1; k++) {
	
	mu[index[2]] = pmp1->val0 * (1.0 + k * pmp1->dval);

	/*
	 * Loop on tau for second curve
	 */

	ptp1 = tau0 + 1;
	for(m=ptp1->minstep,dptr2=d2arr[2]; m < ptp1->maxstep + 1; 
	    m++,dptr2++) {
	  tau[index[2]] = m * ptp1->dval;
	  dptr2->tau = tau[index[2]];
	  dptr2->mu = mu[index[2]];

	  /*
	   * Loop on mu for third curve
	   */

	  pmp2 = mu0 + 2;
	  for(n=-pmp2->nval; n < pmp2->nval + 1; n++) {

	    mu[index[3]] = pmp2->val0 * (1.0 + n * pmp2->dval);

	    /*
	     * Loop on tau for third curve
	     */

	    ptp2 = tau0 + 2;
	    for(p=ptp2->minstep,dptr3=d2arr[3]; p < ptp2->maxstep + 1; 
		p++,dptr3++) {
	      tau[index[3]] = p * ptp2->dval;
	      dptr3->tau = tau[index[3]];
	      dptr3->mu = mu[index[3]];

	      /*
	       * Create the composite curve.
	       */

	      if(!(compos = make_compos(flux,ncurves,npoints,index,tau,mu,
					&ncompos,0)))
		no_error = 0;
	      switch(setup->dispchoice) {
	      case D21M:
		dptr3->disp = disp_d1(compos,ncompos,setup->d2delta);
		break;
	      default:
		fprintf(stderr,"ERROR: four_curve_disp. ");
		fprintf(stderr,"Invalid dispersion method.\n");
		fprintf(stderr," Using D^2_1 method.\n");
		dptr3->disp = disp_d1(compos,ncompos,setup->d2delta);
	      }
	      compos = del_fluxrec(compos);

	      /*
	       * Hold value giving minimum dispersion FOR THIS LOOP
	       */

	      count++;
	      if(count == 1) {

		/*
		 * Initialize d2min
		 */

		*d2min = *dptr3;
		*(d2min+1) = *dptr3;
		*(d2min+2) = *dptr3;

		/*
		 * Open temporary output file
		 */
		
		if(!(tmpfp = fopen("four_curve.update","a")))
		  no_error = 0;

		printf("\nStart at:                   %6.1f %6.4f %6.1f ",
		       tau[index[1]],mu[index[1]],tau[index[2]]);
		printf("%6.4f %6.1f %6.4f %8.4f\n",mu[index[2]],
		       tau[index[3]],mu[index[3]],d2min->disp);
		fprintf(tmpfp,
			"\nStart at %6.1f %8.5f %6.1f %8.5f %6.1f %8.5f %f\n",
			tau[index[1]],mu[index[1]],tau[index[2]],mu[index[2]],
			tau[index[3]],mu[index[3]],d2min->disp);
		if(tmpfp)
		  fclose(tmpfp);
	      }
	      else if(dptr3->disp < d2min->disp) {
		*d2min = *dptr1;
		*(d2min+1) = *dptr2;
		*(d2min+2) = *dptr3;
		d2min->disp = (d2min+1)->disp = (d2min+2)->disp = dptr3->disp;

		/*
		 * Open temporary output file
		 */
		
		if(!(tmpfp = fopen("four_curve.update","a")))
		  no_error = 0;

		printf("(%4d,%4d,%4d,%4d,%4d): %6.1f %6.4f %6.1f %6.4f ",
		       i,j,k,m,n,
		       tau[index[1]],mu[index[1]],tau[index[2]],mu[index[2]]);
		printf("%6.1f %6.4f %8.4f\n",
		       tau[index[3]],mu[index[3]],d2min->disp);
		fprintf(tmpfp,
			"Loop %7d: %6.1f %8.5f %6.1f %8.5f %6.1f %8.5f %f\n",
			count,
			tau[index[1]],mu[index[1]],tau[index[2]],mu[index[2]],
			tau[index[3]],mu[index[3]],d2min->disp);
		if(tmpfp)
		  fclose(tmpfp);
	      }
	    }
	  }
	}
      }
    }
#if 0
    /*
     * Hold two-dimensional slices of dispersion surface containing minimum D^2
     */

    if(i == -NFLUXSTEP) {
      d2min = d2minj;
      for(j=-setup->ntau,dptr=d2arr,mptr1=minarr[1]; j < setup->ntau + 1; 
	  j++,dptr++,mptr1++) {
	*mptr = *dptr;
      }
    }
    else if(d2minj.disp < d2min.disp) {
      d2min = d2minj;
      for(j=-setup->ntau,dptr=d2arr,mptr1=minarr[1]; j < setup->ntau + 1; 
	  j++,dptr++,mptr1++) {
	*mptr1 = *dptr;
      }
    }
#endif
  }

  /*
   * Print out dispersion spectrum with minimum D^2 and best-fitting
   *  parameters
   */

  if(no_error) {
    if(doprint) {
      if(!(ofp = open_writefile(outname))) {
	no_error = 0;
      }
      else {
	for(i=-setup->ntau,mptr1=minarr[1]; i < setup->ntau + 1; i++,mptr1++)
	  fprintf(ofp,"%8.5f %6.1f %f\n",mptr1->mu,mptr1->tau,mptr1->disp);
      }
    }
    printf("    Best-fit Values        \n");
    printf("    ---------------        \n\n");
    printf("  tau       mu       D^2   \n");
    printf(" -------- -------- --------\n");

    for(i=1,dptr1=d2min; i<N08; i++,dptr1++)
      printf("%6.1f %8.5f %f\n",dptr1->tau,dptr1->mu,dptr1->disp);
  }

  /*
   * Clean up and exit
   */

  for(i=0; i<N08; i++) {
    d2arr[i] = del_lcdisp(d2arr[i]);
    minarr[i] = del_lcdisp(minarr[i]);
  }
  d2min = del_lcdisp(d2min);

  if(ofp)
    fclose(ofp);
  if(tmpfp)
    fclose(tmpfp);

  if(no_error)
    return 0;
  else {
    fprintf(stderr,"ERROR: four_curve_disp\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function disp_d1
 *
 * Computes the dispersion of a combined light curve based on two or more
 *  input light curve which have been shifted in time and flux
 *  density.  The dispersion is calculated using the method of
 *  Pelt et al (1994, A&A, 286, 775; 1996, A&A, 305, 97),
 *  with their D^2_1 dispersion.
 *
 * In this method, one curve is shifted (by tau) and scaled (by 1/mu)
 *  and then all of the curves are combined into the composite curve 
 *  C_i(tau,mu).  This C_i(tau,mu) curve is the one that is passed to
 *  this function and called the "compos" curve.  The dispersion, D^2_1, 
 *  is then
 *
 *                     sum(i)_0^(N-1) W_i,i+1 G_i,i+1 (C_i+1 - C_i)^2
 *  D^2_1 (tau) = min  ----------------------------------------------
 *                 mu         2 sum(i)_0^(N-1) W_i,i+1 G_i,i+1 
 *
 *  where
 *
 *            { 1   if i and i+1 come from different original light curves
 *  G_i,i+1 = {
 *            { 0   if i and i+1 come from the same original light curve
 *
 *  and
 *
 *          W_i W_j
 *  W_ij = ---------  (W_i = 1/sigma_i^2)
 *         W_i + W_j
 *
 * Inputs: Fluxrec *compos     composite light curve
 *         int ncompos         number of points in composite light curve.
 *         float delta         timescale for "near pairs" (dummy variable
 *                              for this function, but needed for disp_d2)
 *
 * Output: float d2            D^2, the dispersion
 *
 * v10Oct00 CDF, Added output file printing.
 * v19Apr01 CDF, Moved calculation of composite curve up to the calling
 *                function (e.g., two_curve_disp).  This simplifies the
 *                disp_d1 function call considerably, and allows it to
 *                be much more general.
 *               Moved output file printing up to calling function.
 */

float disp_d1(Fluxrec *compos, int ncompos, float delta)
{
  int i;                    /* Looping variables */
  int npairs=0;             /* Number of valid pairs used in calculation */
  float fdiff;              /* C_i+1 - C_i */
  float wi,wj,wij;          /* Statistical weights */
  float sum=0.0;            /* Weighted sum of (C_i+1 - C_i)^2 */
  float wtsum=0.0;          /* Sum of weights */
  float d2;                 /* Dispersion */
  Fluxrec *cptr;            /* Pointer to navigate compos */

  /*
   * Loop through the composite curve
   */

  for(i=0,cptr=compos; i<ncompos-1; i++,cptr++) {

    /*
     * Check to see if the points cptr and cptr+1 are from the
     *  same original light curve.  If they are NOT, then make the
     *  computations.  Otherwise, go on to the next iteration of this
     *  for loop.
     */

    if(cptr->match != (cptr+1)->match) {
      npairs++;

      /*
       * Calculate difference in flux density
       */

      fdiff = (cptr+1)->flux - cptr->flux;

      /*
       * Calculate weights
       */

      wi = 1.0 / (cptr->err * cptr->err);
      wj = 1.0 / ((cptr+1)->err * (cptr+1)->err);
      wij = wi * wj / (wi + wj);

      /*
       * Add appropriate values to the sums
       */

      sum += (wij * fdiff * fdiff);
      wtsum += wij;
    }    
  }

  /*
   * Calculate the dispersion
   */

  if(npairs == 0) {
    fprintf(stderr,"ERROR: disp_d1.  No valid pairs in curve.\n");
    d2 = 0.0;
  }
  else {
    d2 = sum / (2.0 * wtsum);
  }

  /*
   * Clean up and return
   */

  return d2;
}

/*.......................................................................
 *
 * Function disp_d1
 *
 * Computes the dispersion of a combined light curve based on two or more
 *  input light curve which have been shifted in time and flux
 *  density.  The dispersion is calculated using the method of
 *  Pelt et al (1994, A&A, 286, 775; 1996, A&A, 305, 97),
 *  with their D^2_1 dispersion.
 *
 * In this method, one curve is shifted (by tau) and scaled (by 1/mu)
 *  and then all of the curves are combined into the composite curve 
 *  C_i(tau,mu).  This C_i(tau,mu) curve is the one that is passed to
 *  this function and called the "compos" curve.  The dispersion, D^2_1, 
 *  is then
 *
 *                     sum(i)_0^(N-1) W_i,i+1 G_i,i+1 (C_i+1 - C_i)^2
 *  D^2_1 (tau) = min  ----------------------------------------------
 *                 mu         2 sum(i)_0^(N-1) W_i,i+1 G_i,i+1 
 *
 *  where
 *
 *            { 1   if i and i+1 come from different original light curves
 *  G_i,i+1 = {
 *            { 0   if i and i+1 come from the same original light curve
 *
 *  and
 *
 *          W_i W_j
 *  W_ij = ---------  (W_i = 1/sigma_i^2)
 *         W_i + W_j
 *
 * Inputs: Fluxrec *compos     composite light curve
 *         int ncompos         number of points in composite light curve.
 *         float delta         timescale for "near pairs" (dummy variable
 *                              for this function, but needed for disp_d2)
 *
 * Output: float d2            D^2, the dispersion
 *
 * v10Oct00 CDF, Added output file printing.
 * v19Apr01 CDF, Moved calculation of composite curve up to the calling
 *                function (e.g., two_curve_disp).  This simplifies the
 *                disp_d1 function call considerably, and allows it to
 *                be much more general.
 *               Moved output file printing up to calling function.
 */

float disp_d1b(Fluxrec *compos, int ncompos, float delta)
{
  int i;                    /* Looping variables */
  int npairs=0;             /* Number of valid pairs used in calculation */
  float fdiff;              /* C_i+1 - C_i */
  float wi,wj,wij;          /* Statistical weights */
  float sum=0.0;            /* Weighted sum of (C_i+1 - C_i)^2 */
  float wtsum=0.0;          /* Sum of weights */
  float d2;                 /* Dispersion */
  Fluxrec *cptr;            /* Pointer to navigate compos */

  /*
   * Loop through the composite curve
   */

  for(i=0,cptr=compos; i<ncompos-1; i++,cptr++) {

    /*
     * Check to see if the points cptr and cptr+1 are from the
     *  same original light curve.  If they are NOT, then make the
     *  computations.  Otherwise, go on to the next iteration of this
     *  for loop.
     */

    if(cptr->match != (cptr+1)->match) {
      npairs++;

      /*
       * Calculate difference in flux density
       */

      fdiff = (cptr+1)->flux - cptr->flux;

      /*
       * Calculate weights
       */

      wi = 1.0 / (cptr->err * cptr->err);
      wj = 1.0 / ((cptr+1)->err * (cptr+1)->err);
      wij = wi * wj / (wi + wj);

      /*
       * Add appropriate values to the sums
       */

      sum += (wij * fdiff * fdiff);
      wtsum += wij;
    }    
  }

  /*
   * Calculate the dispersion
   */

  if(npairs == 0) {
    fprintf(stderr,"ERROR: disp_d1.  No valid pairs in curve.\n");
    d2 = 0.0;
  }
  else {
    d2 = sum / (2.0 * wtsum);
  }

  /*
   * Clean up and return
   */

  return d2;
}

/*.......................................................................
 *
 * Function disp_d2
 *
 * Computes the dispersion of a combined light curve based on two
 *  input light curve which have been shifted in time and flux
 *  density.  The dispersion is calculated using the method of
 *  Pelt et al (1994, A&A, 286, 775; 1996, A&A, 305, 97),
 *  with their D^2_2 dispersion.
 *
 * In this method, one curve is shifted (by tau) and scaled (by 1/mu)
 *  and then the two curves are combined into the composite curve 
 *  C_i(tau,mu).  This C_i(tau,mu) curve is the one that is passed to
 *  this function and called the "compos" curve.  The dispersion, D^2_2, 
 *  is then
 *
 *                     sum_n sum_m W_n,m S_n,m G_n,m (C_n - C_m)^2
 *  D^2_2 (tau) = min  ----------------------------------------------
 *                 mu         2 sum_n sum_m W_n,m S_n,m G_n,m 
 *
 *  where the sum on n goes from 0 to N-1 and the sum on m goes from
 *  n+1 to N.  Additionally
 *
 *          { 1   if n and m come from different original light curves
 *  G_n,m = {
 *          { 0   if n and m come from the same original light curve,
 *
 *
 *           W_n W_m
 *  W_n,m = ---------  (W_i = 1/sigma_i^2),
 *          W_n + W_m
 *
 *  and S_n,m is the "nearness" weighting coeffiecient:
 *
 *          { (1 - |t_n - t_m|)/delta    if |t_n - t_m| <= delta
 *  S_n,m = {
 *          {           0                otherwise
 *
 *
 * Inputs: Fluxrec *compos     composite light curve
 *         int ncompos         number of points in composite light curve.
 *         float delta         timescale for "near pairs"
 *
 * Output: float d2            D^2, the dispersion
 *
 * v10Oct00 CDF, Added output file printing.
 * v19Apr01 CDF, Moved calculation of composite curve up to the calling
 *                function (e.g., two_curve_disp).  This simplifies the
 *                disp_d2 function call considerably, and allows it to
 *                be much more general.
 *               Moved output file printing up to calling function.
 */

float disp_d2(Fluxrec *compos, int ncompos, float delta)
{
  int n,m;                  /* Looping variables */
  int npairs=0;             /* Number of valid pairs used in calculation */
  float tdiff;              /* |t_n - t_m| */
  float fdiff;              /* C_n - C_m */
  float wn,wm,wnm;          /* Statistical weights */
  float snm;                /* "Nearness" weight */
  float sum=0.0;            /* Weighted sum of (C_n - C_m)^2 */
  float wtsum=0.0;          /* Sum of weights */
  float d2;                 /* Dispersion */
  Fluxrec *cptr1,*cptr2;    /* Pointers to navigate compos */

  /*
   * Loop through the composite curve
   */

  for(n=0,cptr1=compos; n<ncompos-1; n++,cptr1++) {
    for(m=n+1,cptr2=compos+m; m<ncompos; m++,cptr2++) {

      /*
       * Calculate tdiff
       */

      tdiff = fabs(cptr1->day - cptr2->day);

      /*
       * Check to see if the points cptr1 and cptr2 are from the
       *  same original light curve.  If they are not AND if 
       *  tdiff < delta, then make the computations.  Otherwise, go on to 
       *  the next iteration of this for loop.
       */

      if((cptr1->match != cptr2->match) && tdiff <= delta) {
	npairs++;

	/*
	 * Calculate difference in flux density
	 */

	fdiff = cptr1->flux - cptr2->flux;

	/*
	 * Calculate statistical weights
	 */

	wn = 1.0 / (cptr1->err * cptr1->err);
	wm = 1.0 / (cptr2->err * cptr2->err);
	wnm = wn * wm / (wn + wm);

	/*
	 * Calculate "nearness" weight
	 */

	snm = 1 - (tdiff/delta);

	/*
	 * Add appropriate values to the sums
	 */

	sum += (wnm * snm * fdiff * fdiff);
	wtsum += (wnm * snm);
      }
    }   
  }

  /*
   * Calculate the dispersion
   */

  d2 = sum / (2.0 * wtsum);

  /*
   * Clean up and return
   */

  return d2;
}

/*.......................................................................
 *
 * Function disp_lovell
 *
 * Computes the dispersion of a combined light curve based on two
 *  input light curve which have been shifted in time and flux
 *  density.  The dispersion is calculated using the method of
 *  Pelt et al (1994, A&A, 286, 775; 1996, A&A, 305, 97) as
 *  modified by Lovell et al. in the 1830 time delay paper.
 *
 * In this method, the dispersion, D^2, is
 *
 *                 sum(i,j) W_ij V'_ij (a_i - b_j)^2
 *  D^2 (tau,mu) = ----------------------------------
 *                   2 sum(i,j) W_ij V'_ij
 *
 *  where
 *
 *          { 1                                   if |t_i - t_j| < delta
 *  v'_ij = {
 *          {
 *          { (1 + ((delta - |t_i - t_j|)/(delta/2))^2)^{-1}  otherwise
 *
 *  and
 *
 *          W_i W_j
 *  W_ij = ---------  (W_i = 1/sigma_i^2)
 *         W_i + W_j
 *
 *  and
 *
 *  a_i is the first light curve and b_j is the second light curve, shifted
 *   in time (tau) and multiplied by 1/mu.
 *
 * Inputs: Fluxrec *a          first light curve
 *         Fluxrec *b          second light curve
 *         int npoints         number of points in light curves
 *         float tau           time delay
 *         float mu            flux density ratio
 *         float delta         timescale for "near pairs"
 *
 * Output: float d2            D^2, the dispersion
 *
 */

float disp_lovell(Fluxrec *a, Fluxrec *b, int npoints, float tau, float mu,
		  float delta)
{
  int i,j;              /* Looping variables */
  float tdiff;          /* |t_i - t_j| */
  float wi,wj,wij;      /* Statistical weights */
  float vij;            /* "Nearness weight */
  float abdiff;         /* a_i - b_j */
  float sum=0.0;        /* Weighted sum of (a_i - b_j)^2 */
  float wtsum=0.0;      /* Sum of weights */
  float d2;             /* Dispersion */
  Fluxrec *aptr,*bptr;  /* Pointers for navigating a and b */

  for(i=0,aptr=a; i<npoints; i++,aptr++) {
    for(j=0,bptr=b; j<npoints; j++,bptr++){

      /*
       * Calculate differences
       */

      tdiff = fabs(aptr->day - bptr->day + tau);
      abdiff = aptr->flux - bptr->flux / mu;

      /*
       * Calculate weights
       */

      wi = 1.0 / (aptr->err * aptr->err);
      wj = (mu * mu) / (bptr->err * bptr->err);
      if(aptr->match == -1  || bptr->match == -1)
	wij = 0.0;
      else
	wij = wi * wj / (wi + wj);
      if(tdiff < delta)
	vij = 1.0;
      else
	vij = 1.0 / (1.0 + 4.0*(delta-tdiff)*(delta-tdiff)/(delta*delta));

      /*
       * Add appropriate values to the sums
       */

      sum += (wij * vij * abdiff * abdiff);
      wtsum += (wij * vij);
    }    
  }

  /*
   * Calculate the dispersion
   */

  d2 = sum / (2.0 * wtsum);

  return d2;
}

/*.......................................................................
 *
 * Function print_disp_slice
 *
 * Prints out the dispersion as a function of tau.  This is a slice
 *  from the (tau, mu) plane, hence the name of the function.
 *
 * Inputs: LCdisp *disp        structure containing dispersion information
 *         int size            size of array
 *         char *outname       output file name
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 */

int print_disp_slice(LCdisp *disp, int size, char *outname)
{
  int i;           /* Looping variable */
  LCdisp *dptr;    /* Pointer to navigate disp */
  FILE *ofp=NULL;  /* Output file pointer */

  /*
   * Open output file
   */

  printf(" ");
  if(!(ofp = open_writefile(outname))) {
    fprintf(stderr,"ERROR: print_disp_slice.\n");
    return 1;
  }

  /*
   * Print header
   */

  fprintf(ofp,"# tau     mu     disp   \n");
  fprintf(ofp,"#------ ------ ---------\n");

  /*
   * Print data
   */

  for(i=0,dptr=disp; i<size; i++,dptr++)
    fprintf(ofp,"%7.2f %6.4f %7.4f\n",dptr->tau,dptr->mu,
	    dptr->disp);

  /*
   * Clean up and exit
   */

  if(ofp)
    fclose(ofp);
  return 0;
}

/*.......................................................................
 *
 * Function call_dcf
 *
 * Finds discrete correlation functions (DCFs) between pairs of light curves, 
 *  by calling the function discrete_corr for each pair of curves.  Puts the
 *  output into a file designated by filename.
 *
 * Inputs: Fluxrec *flux[]     light curves
 *         int size            number of points in each curve
 *         char *filename      name of output file
 *         FILE *logfp         logfile pointer
 *         int doprint         flag set to 1 if output files are desired
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> failure
 *
 */

int call_dcf(Fluxrec *flux[], int size, char *filename, FILE *logfp, 
	     int doprint)
{
  int i;                       /* Looping variable */
  int no_error=1;              /* Flag set to 0 on error */
  int corsize;                 /* Length of correlation curves */
  float maxlag=200.0;          /* Maximum lag in discrete routines */
  float dlag=4.0;             /* Bin width for steps in DCF */
  Fluxrec best_ba;             /* Highest discrete correlation for B-A */
  Fluxrec best_bc;             /* Highest discrete correlation for B-C */
  Fluxrec best_bd;             /* Highest discrete correlation for B-D */
  Fluxrec *zmean[N08]={NULL};  /* Zero-mean light curves */
  Fluxrec *corrbc=NULL;        /* B-C discrete correlation */
  Fluxrec *corrba=NULL;        /* B-A discrete correlation */
  Fluxrec *corrbd=NULL;        /* B-D discrete correlation */
  Fluxrec *cp1,*cp2,*cp3;      /* Pointers to navigate correlation arrays */
  FILE *ofp=NULL;              /* Output file */

  /*
   * Initialize the best lag containers to something ridiculous.  In this
   *  function the lag will be contained in the day member of the Fluxrec
   *  structure and the discrete correlation coefficient will be contained 
   *  in the flux member.
   */

  best_ba.flux = -999.0;
  best_bc.flux = -999.0;
  best_bd.flux = -999.0;

  /*
   * Normalize the curves since we are looking for correlations in
   *  fractional variations and then subtract 1 to create a zero-mean
   *  light curve.  This helps by getting rid of the DC component in
   *  the discrete correlation.
   */

  for(i=0; i<N08; i++) {
    if(no_error)
      if(!(zmean[i] = norm_zero_mean(flux[i],size)))
	no_error = 0;
  }

  /*
   * Set the size of the discrete correlation function arrays
   */

  corsize = ((int) maxlag/dlag) * 2 + 1;
  printf("\ncall_dcf: Max lag = %6.2f, d_lag = %7.4f, ndcf = %d\n",
	 maxlag,dlag,corsize);

  /*
   * Call the discrte correlation routines for each pair of
   *  zero-mean light curves.
   */

  if(no_error)
    if(!(corrbc = discrete_corr(zmean[1],zmean[2],size,dlag,maxlag,corsize)))
      no_error = 0;

  if(no_error)
    if(!(corrba = discrete_corr(zmean[1],zmean[0],size,dlag,maxlag,corsize)))
      no_error = 0;

  if(no_error)
    if(!(corrbd = discrete_corr(zmean[1],zmean[3],size,dlag,maxlag,corsize)))
      no_error = 0;

  /*
   * Open output file
   */

  if(doprint && no_error)
    if((ofp = fopen(filename,"w")) == NULL) {
      fprintf(stderr,"ERROR: call_dcf.  Cannot open %s\n",filename);
      no_error = 0;
    }

  /*
   * Print discrete correlation curves to output files and also find
   *  lags with highest discrete correlation values.
   */

  if(no_error) {
    if(doprint) {
      fprintf(ofp,"#  Lag    Corr_ba   Corr_bc   Corr_bd\n");
      fprintf(ofp,"#------- --------- --------- ---------\n");
    }
    for(i=0,cp1=corrba,cp2=corrbc,cp3=corrbd; i<corsize; 
	i++,cp1++,cp2++,cp3++) {
      if(doprint)
	fprintf(ofp,"%8.3f %9f %9f %9f\n",
	      -cp1->day,cp1->flux,cp2->flux,cp3->flux);
      if(cp1->flux > best_ba.flux)
	best_ba = *cp1;
      if(cp2->flux > best_bc.flux)
	best_bc = *cp2;
      if(cp3->flux > best_bd.flux)
	best_bd = *cp3;
    }
  }

  /*
   * Print best lags to logfile, if it exists.  Take the negative of
   *  best_b?.day because if B leads the other light curve, the
   *  discrete correlation will return a negative number for the best-fit
   *  lags.  Since we know/hope that B leads all other values and we
   *  are looking for the amount that the other curves lag B, just
   *  take the negative of best_b?.day.
   */

  if(logfp && no_error) {
    fprintf(logfp,"# Discrete correlation results\n");
    fprintf(logfp,"%7.2f %f ",-best_ba.day,best_ba.flux);
    fprintf(logfp,"%7.2f %f ",-best_bc.day,best_bc.flux);
    fprintf(logfp,"%7.2f %f\n",-best_bd.day,best_bd.flux);
  }

  /*
   * Also print output to screen
   */

  if(no_error) {
    printf("\n call_dcf: B-A best lag = %7.2f d (%f)\n",
	   -best_ba.day,best_ba.flux);
    printf(" call_dcf: B-C best lag = %7.2f d (%f)\n",
	   -best_bc.day,best_bc.flux);
    printf(" call_dcf: B-D best lag = %7.2f d (%f)\n",
	   -best_bd.day,best_bd.flux);
  }

  /*
   * Clean up
   */

  for(i=0; i<N08; i++)
    zmean[i] = del_fluxrec(zmean[i]);
  corrbc = del_fluxrec(corrbc);
  corrba = del_fluxrec(corrba);
  corrbd = del_fluxrec(corrbd);
  if(ofp)
    fclose(ofp);

  if(no_error)
    return 0;
  else {
    fprintf(stderr,"ERROR: call_dcf\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function discrete_corr
 *
 * Calculates the discrete correlation function (Edelson & Krolik 1988,
 *  ApJ, 333, 646) of two unevenly sampled time series.  This method
 *  does not require any interpolation of the light curves.  For the
 *  two time series a_i and b_j, define the unbinned discrete correlations
 *
 *                      (a_i - <a>)(b_j - <b>)
 *    UDCF_ij = -------------------------------------------
 *              sqrt((sigma_a^2 - e_a^2)(sigma_b^2 - e_b^2)
 *
 *  the binned discrete correlation then becomes
 *
 *    DCF(tau) = sum(UCDF_ij) / npair
 *
 *  where the sum is over the npair pairs for which delta t_ij = t_j - t_i
 *  is in the range (tau - (delta tau)/2) < t_ij < (tau + (delta tau)/2).
 *
 * NB: delete points from the DCF which have t_ij = 0.
 *
 * Inputs: Fluxrec *flux1      first light curve
 *         Fluxrec *flux2      second light curve
 *         int npoints         number of points in light curves
 *         float binsize       width of bin in lag space
 *         float maxlag        maximum lag
 *         int ndcf            number of points in the dcf curve
 *
 * Output: Fluxrec *dcf        discrete correlation function
 *
 */

Fluxrec *discrete_corr(Fluxrec *flux1, Fluxrec *flux2, int npoints, 
		       float binsize, float maxlag, int ndcf)
{
  int i,j;                 /* Looping variables */
  int no_error=1;          /* Flag set to 0 on error */
  int npair;               /* Number of pairs contained in any bin */
  Fluxrec *fptr1,*fptr2;   /* Pointers for navigating flux1 and flux2 */
  Fluxrec *udcf=NULL;      /* Container for unbinned discrete correlations */
  Fluxrec *uptr;           /* Pointer for navigating udcf */
  Fluxrec *dcf=NULL;       /* Container for binned discrete correlations */
  Fluxrec *dptr;           /* Pointer for navigating dcf */

  /*
   * Allocate memory for the containers
   */

  if(!(udcf = new_fluxrec(npoints*npoints))) {
    fprintf(stderr,"ERROR: discrete_corr\n");
    return NULL;
  }

  if(!(dcf = new_fluxrec(ndcf)))
    no_error = 0;

  /*
   * Calculate the unbinned discrete correlations.
   * NB: Since we have made the data sets zero-mean, we don't have to
   *      subtract the means in this calculation.
   * NB: Ignore the normalization factor for now.  Instead normalize by
   *      1/(err1 * err2)
   */

  if(no_error) {
    uptr = udcf;
    for(i=0,fptr1=flux1; i<npoints; i++,fptr1++) {
      for(j=0,fptr2=flux2; j<npoints; j++,fptr2++) {
	uptr->day = fptr2->day - fptr1->day;
	uptr->flux = fptr1->flux * fptr2->flux / (fptr1->err * fptr2->err);
	uptr++;
      }
    }
  }

  /*
   * Bin the UDCF results
   */

  if(no_error) {
    for(i=0,dptr=dcf; i<ndcf; i++,dptr++) {
      dptr->day = -maxlag + i*binsize;
      npair=0;
      dptr->flux = dptr->err = 0.0;
      for(j=0,uptr=udcf; j<npoints*npoints; j++,uptr++) {
	if(uptr->day - dptr->day != 0 && 
	   fabs(uptr->day - dptr->day) < binsize/2.0) {
	  dptr->flux += uptr->flux;
	  npair++;
	}
      }
      if(npair > 0)
	dptr->flux /= npair;
    }
  }

  /*
   * Clean up and exit
   */

  udcf = del_fluxrec(udcf);

  if(no_error)
    return dcf;
  else {
    fprintf(stderr,"ERROR: discrete_corr\n");
    return NULL;
  }
}

/*.......................................................................
 *
 * Function fit_delays
 *
 * This function is called by fit_1608.  It loops over all of the possible
 *  delays and creates a composite curve for each combination of delays.
 *  The composite curve is fit with fit_poly and the reduced chisq of the
 *  fit is compared to the best chisq found.  The values of the parameters
 *  associated with each improved fit are recorded.
 *
 * Inputs: Fluxrec *flux[]     input light curves
 *         int npoints         number of points in input light curves
 *         float mean[]        means of input light curves
 *         float d0[]          initial guesses for delays
 *         float dtmax         maximum delay allowed relative to initial guess
 *         float dt            step size for incrementing delays
 *         int ndeg            degree of function to fit
 *         Fluxrec *compos     composite curve container
 *         int ncompos         number of points in composite curve
 *         int *first          flag set to 1 for first iteration, 0 otherwise
 *         float *bestchi      current minimum reduced chisq
 *         Fluxrec bestpars[]  container for parameters (delays and fluxes) 
 *                              associated with the lowest reduced chisq
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 */
#if 0
int fit_delays(Fluxrec *flux[], int npoints, float mean[], float d0[], 
	       float dtmax, float dt, int ndeg, Fluxrec *compos, int ncompos, 
	       int *first, float *bestchi, Fluxrec bestpars[])
{
  int i,j;           /* Looping variables */
  int no_error=1;    /* Flag set to 0 on error */
  float chisq;       /* Reduced chisq for fit */
  float delays[N08]; /* Current value of delays */
  float *a=NULL;     /* Container for coefficients of fitting fn. */
  Fluxrec *compptr;  /* Pointer to navigate compos */
  Fluxrec *fptr;     /* Pointer to navigate input light curves */

  /*
   * Set the B delay to 0
   */

  delays[1] = d0[1];

  /*
   * Step through all possible combinations of delays given dtmax and dt.
   * Make the outermost loop on BA delay.
   */

  delays[0] = d0[0] - dtmax;
  while(delays[0] < d0[0] + dtmax+0.5*dt && no_error) {

    /*
     * Initialize BC delay for this iteration
     */

    delays[2] = d0[2] - dtmax;

    /*
     * Loop on BC delay
     */

    while(delays[2] < d0[2]+dtmax+0.5*dt && no_error) {

      /*
       * Initialize BD delay for this iteration
       */

      delays[3] = d0[3] - dtmax;

      /*
       * Loop on BD delay
       */

      while(delays[3] < d0[3]+dtmax+0.5*dt && no_error) {

	/*
	 * Fill the composite curve for this iteration
	 */

	compptr = compos;
	for(i=0; i<N08; i++){
	  for(j=0,fptr=flux[i]; j<npoints; j++,fptr++) {
	    if(fptr->match > -1) {
	      *compptr = *fptr;
	      compptr->day -= delays[i];
	      compptr->flux /= mean[i];
	      compptr->err /= mean[i];
	      compptr++;
	    }
	  }
	}

	/*
	 * Call the fitting function
	 */

	if(!(a = fit_poly(compos,ncompos,11,0)))
	  no_error = 0;

	/*
	 * Now calculate the reduced chisq
	 */

	if(no_error)
	  if((chisq = chisq_fit(compos,ncompos,a,ndeg+1,7,0)) < 0.0)
	    no_error = 0;

	if((chisq < *bestchi || *first) && no_error) {
	  *first = 0;
	  *bestchi = chisq;
	  printf("Reduced chisq = %f ",*bestchi);
	  for(i=0; i<N08; i++) {
	    bestpars[i].day = delays[i];
	    printf("%5.1f ",delays[i]);
	  }
	  for(i=0; i<N08; i++) {
	    bestpars[i].flux = mean[i];
	    printf("%6.3f ",mean[i]);
	  }
	  printf("\n");
	}
	delays[3] += dt;
      }
      delays[2] += dt;
    }
    delays[0] += dt;
  }

  printf("fit_delays: finished loop.\n");

  if(no_error)
    return 0;
  else{
    fprintf(stderr,"ERROR: fit_delays\n");
    return 1;
  }
}
#endif
