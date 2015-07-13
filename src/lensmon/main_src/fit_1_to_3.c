/*
 * fit_1_to_3.c
 *
 * Usage: fit_1_to_3 lcurve [1608 filename] [1634/1635 filename]
 *
 * Creates a composite light curve from three of the four 1608 lightcurves
 *  and interactively fits a function to them.  Then shifts the fourth curve 
 *  in steps, calculating the goodness of fit of the free curve with respect 
 *  to the fitted function at each step.
 *
 * 30Sep98 CDF
 * v01Oct98 CDF, Added sorting of compos to compos_3
 * v02Dec98 CDF, Moved actual creation of composite curve into make_compos
 *                function in lc_funcs.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "structdef.h"
#include "dataio.h"
#include "cpgplot.h"
#include "plotfuncs.h"
#include "lc_funcs.h"

/*.......................................................................
 *
 * Function declarations
 *
 */

Fluxrec *compos_3(Fluxrec *flux[], int npoints, int nbad[], int *ncompos,
		  int *shiftind);
float *fit_compos(Fluxrec *compos, int ncompos, int *ndeg);
int plot_fit_curve(float *a, int ncoeff, float tmin, float tmax);
int plot_resid(Fluxrec *flux, int npoints, float *a, int ncoeff,
	       float *meanrms);
int find_best_shift(Fluxrec *flux, int npoints, float *a, int ncoeff, 
		    float fitrms, float stepsize, float tmin, float tmax,
		    float *bestlag);
int plot_best_shift(Fluxrec *flux, int npoints, float mean, float *a, 
		    int ncoeff, float tmin, float tmax, float bestlag);

/*.......................................................................
 *
 * Main program
 *
 */

int main(int argc, char *argv[])
{
  int i;                      /* Looping variable */
  int no_error=1;             /* Flag set to 0 on error */
  int nlines=0;               /* Number of data lines in lens input file */
  int ncompos;                /* Number of points in composite curve */
  int flagbad=1;              /* Flag set to 1 for bad day flagging */
  int errtype=0;              /* 1 ==> include fracrms in errors */
  int ndeg;                   /* Degree of fitting function */
  int nbad[N08];              /* Number of bad days in bad day file */
  int shiftind;               /* Curve excluded from composite */
  float fracrms;              /* Fractional rms scatter in 1634/1635 ratio */
  float tmin,tmax;            /* Min and max times used for fit */
  float meanrms;              /* RMS about mean for residual curve */
  float bestlag;              /* Lag giving lowest reduced chisq */
  float mean[N08];            /* Means of the 4 light curves */
  float *flat=NULL;           /* "Flat field" light curve */
  float *a=NULL;              /* Coefficients of fitting function */
  char line[MAXC];            /* General string for reading input */
  char lensfile[MAXC];        /* Input file name for 1608 data */
  char calfile[MAXC];         /* Input file name for 1634 and 1635 data */
  Fluxrec *fl08[N08]={NULL};  /* Raw 1608 light curves */
  Fluxrec *ff08[N08]={NULL};  /* Flat-fielded 1608 light curves */
  Fluxrec *fl34=NULL;         /* 1634 light curve */
  Fluxrec *fl35=NULL;         /* 1635 light curve */
  Fluxrec *compos=NULL;       /* Composite curve created from 3 light curves */
  FILE *lfp=NULL;             /* 1608 data file pointer */
  FILE *cfp=NULL;             /* 1634/1635 data file pointer */

  /*
   * Check input line
   */
  if(argc < 3) {
    fprintf(stderr,"\nUsage: fit_1_to_3 [1608 filename] [1634/1635 filename]");
    fprintf(stderr,"\n\n");
    return 1;
  }

  /*
   * Open 1608 and secondary calibrator files.
   */

  strcpy(lensfile,argv[1]);
  strcpy(calfile,argv[2]);

  if(!(lfp = open_readfile(lensfile)))
    no_error = 0;
  if(!(cfp = open_readfile(calfile)))
    no_error = 0;

  /*
   * First count the number of data input lines in order to set
   *  sizes of arrays.
   */

  if(no_error)
    nlines = n_lines(lfp,'#');

  if(nlines == 0)
    no_error = 0;
  else if(nlines != n_lines(cfp,'#')) {
    fprintf(stderr,
	    "ERROR. Number of lines in lens and cal files don't match.\n");
    no_error = 0;
  }
  else {
    printf("\n%d data lines in lens and cal input files.\n\n",nlines);
    rewind(lfp);
    rewind(cfp);
  }

  /*
   * Allocate arrays
   */

  if(no_error)
    for(i=0; i<N08; i++) {
      if(!(fl08[i] = new_fluxrec(nlines)))
	no_error = 0;
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
    fgets(line,MAXC,stdin);
    if(line[0] != '\n') {
      while(sscanf(line,"%d",&flagbad) != 1 || flagbad < 0) {
	fprintf(stderr,"ERROR: bad input.  Enter value again:  ");
	fgets(line,MAXC,stdin);
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
	printf("Bad days flagged: %d\n",nbad[0]);
    }
  }

  /*
   * Make the "flat field" curve
   */

  if(no_error)
    if(!(flat = make_flat(fl34,fl35,nlines)))
      no_error = 0;

  /*
   * Set error determination
   */

  if(no_error) {
    printf("\nEnter choice of error determination\n");
    printf("  0. RMS noise in map only\n");
    printf("  1. Include fractional RMS from 1634/1635 ratio\n");
    printf("Choice? [%d] ",errtype);
    fgets(line,MAXC,stdin);
    if(line[0] != '\n') {
      while(sscanf(line,"%d",&errtype) != 1 || errtype < 0) {
	fprintf(stderr,"ERROR: Bad input.  Enter new value:  ");
	fgets(line,MAXC,stdin);
      }
    }
  }

  /*
   * Calculate fractional error in 1634/1635 flux ratio
   */

  if(no_error) {
    if(errtype == 0)
      fracrms = 0.0;
    else
      if(ratio_err(fl34,fl35,nlines,&fracrms))
	no_error = 0;
  }

  /*
   * Flat field the 1608 light curves
   */

  if(no_error) {
    for(i=0; i<N08; i++) {
      if(!(ff08[i] = flat_field(fl08[i],flat,nlines,fracrms,NULL,0)))
	no_error = 0;
    }
  }

  /*
   * Set the curve means to their default values
   */

  mean[0] = AMEAN;
  mean[1] = BMEAN;
  mean[2] = CMEAN;
  mean[3] = DMEAN;

  /*
   * Create the three-curve composite
   */

  if(no_error)
    if(!(compos = compos_3(ff08,nlines,nbad,&ncompos,&shiftind)))
      no_error = 0;

  /*
   * Get tmin and tmax from compos, which is easy since it is sorted
   */

  if(no_error) {
    tmin = compos[0].day;
    tmax = compos[ncompos-1].day;
    printf("tmin = %6.2f, tmax = %6.2f\n",tmin,tmax);
  }

  /*
   * Open the plot
   */

  if(no_error)
    plot_open(1.5);

  /*
   * Fit a function to the composite curve
   */

  if(no_error)
    if(!(a = fit_compos(compos,ncompos,&ndeg)))
      no_error = 0;

  /*
   * Summarize results
   */

  printf("\nUsing a fit of degree %d\n",ndeg);

  /*
   * Find and plot residuals
   */

  if(no_error)
    if(plot_resid(compos,ncompos,a,ndeg+1,&meanrms))
      no_error = 0;

  /*
   * Now shift the free curve in steps along the composite curve and
   *  create a new composite curve from the free + old composite.
   * Calculate a reduced chisq at each step by comparing the new composite
   *  with the fitted function.
   */

  if(no_error)
    if(find_best_shift(ff08[shiftind],nlines,a,ndeg+1,meanrms,1.0,tmin,tmax,
		       &bestlag))
      no_error = 1;

  /*
   * Plot out best fit
   */

  if(no_error)
    if(plot_best_shift(ff08[shiftind],nlines,mean[shiftind],a,ndeg+1,
		       tmin,tmax,bestlag))
      no_error = 0;

  /*
   * Clean up
   */

  if(no_error)
    printf("\nCleaning up\n\n");

  plot_close();
  a = del_array(a);
  fl34 = del_fluxrec(fl34);
  fl35 = del_fluxrec(fl35);
  compos = del_fluxrec(compos);
  for(i=0; i<4; i++) {
    fl08[i] = del_fluxrec(fl08[i]);
    ff08[i] = del_fluxrec(ff08[i]);
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

/*.......................................................................
 *
 * Function compos_3
 *
 * Makes a composite curve using 3 out of the 4 light curves
 *
 * Inputs: Fluxrec *flux[]     1608 light curves
 *         int npoints         number of points in light curves
 *         int nbad[]          number of bad points in each curve
 *         int *ncompos        number of points in composit curve (set by 
 *                              this function)
 *         int *shiftind       index of the fourth light curve, which will
 *                              shift w/r.t. the composite (set by this  
 *                              function).  0 = A, 1 = B, etc.
 */

Fluxrec *compos_3(Fluxrec *flux[], int npoints, int nbad[], int *ncompos,
		  int *shiftind)
{
  int i,j;              /* Looping variables */
  int index[3];         /* Array of indices */
  int nptarr[3];        /* Array version of npoints */
  float lag[4];         /* Lags for the 4 light curves */
  float mean[4];        /* Mean values for the 4 light curves */
  char line[MAXC];      /* General string for reading input */
  Fluxrec *compos=NULL; /* Composite light curve */
  Fluxrec *fptr;        /* Pointers to navigate individual light curves */
  Fluxrec *compptr;     /* Pointer to navigate compos */

  /*
   * Fill lag and mean arrays
   */

  lag[0] = ALAG;
  lag[1] = 0.0;
  lag[2] = CLAG;
  lag[3] = DLAG;
#if 0
  printf("\n");
  for(i=0; i<4; i++) {
    if(calc_swmean(flux[i],npoints,&mean[i])) {
      fprintf(stderr,"ERROR: fit_1608\n");
      return NULL;
    }
    printf("compos_3: Mean = %8.5f\n",mean[i]);
  }
#else
  mean[0] = AMEAN;
  mean[1] = BMEAN;
  mean[2] = CMEAN;
  mean[3] = DMEAN;
#endif

  /*
   * Find which light curve will be excluded from the composite
   */

  printf("\nThis function creates a composite from three of the four ");
  printf("1608 light curves.\n");
  printf("  0. A\n");
  printf("  1. B\n");
  printf("  2. C\n");
  printf("  3. D\n");
  printf("Which do you want to exclude?  ");
  fgets(line,MAXC,stdin);
  while(sscanf(line,"%d",shiftind) != 1 || *shiftind < 0 || *shiftind > 3) {
    fprintf(stderr,"ERROR: Bad input.  Enter new value:  ");
    fgets(line,MAXC,stdin);
  }

  /*
   * Set indices to exclude the chosen curve
   */

  switch(*shiftind) {
  case 0:
    index[0] = 1;
    index[1] = 2;
    index[2] = 3;
    break;
  case 1:
    index[0] = 0;
    index[1] = 2;
    index[2] = 3;
    break;
  case 2:
    index[0] = 0;
    index[1] = 1;
    index[2] = 3;
    break;
  case 3:
    index[0] = 0;
    index[1] = 1;
    index[2] = 2;
    break;
  default:
    fprintf(stderr,"ERROR: compos_3.  Not a valid option\n");
    return NULL;
  }

  /*
   * Actually create the composite curve by calling make_compos
   */

  for(i=0; i<3; i++)
    nptarr[i] = npoints;

  if(!(compos = make_compos(flux,3,nptarr,index,lag,mean,ncompos,1))) {
    fprintf(stderr,"ERROR: compos_3\n");
    return NULL;
  }
  
  /*
   * Clean up and return
   */

  return compos;

}

/*.......................................................................
 *
 * Function fit_compos
 *
 * Interactively fits a function to the composite curve.
 *
 * Inputs: Fluxrec *compos     composite curve
 *         int ncompos         number of points in composite curve
 *         int *ndeg           degree of fit (set by this function)
 *
 * Output: float *a            array containing coefficients of fit.
 *                             NULL on error.
 *
 */

float *fit_compos(Fluxrec *compos, int ncompos, int *ndeg)
{
  int no_error=1;   /* Flag set to 0 on error */
  int contin=1;     /* Flag set to 0 to break out of loop */
  float chisq;      /* Reduced chisq calculated in curve fitting */
  float *a=NULL;    /* Container for coefficients of fit */
  char line[MAXC];  /* General string for reading input */

  while(contin && no_error){

    /*
     * Plot points
     */

    if(no_error)
      if(plot_lcurve(compos,ncompos,"Days","Flux Density",""))
	no_error = 0;

    /*
     * Get the "degree" of the fitting function
     */

    if(no_error) {
      printf("\nEnter the desired degree of the fitting function: ");
      fgets(line,MAXC,stdin);
      while(sscanf(line,"%d",ndeg) != 1 || *ndeg < 1) {
	fprintf(stderr,"ERROR: Bad input.  Enter value again: ");
	fgets(line,MAXC,stdin);
      }
    }

    /*
     * Call the fitting function
     */

    if(no_error)
      if(!(a = fit_poly(compos,ncompos,*ndeg+1,1)))
	no_error = 0;

    /*
     * Plot fitted polynomial
     */

    if(no_error)
      if(plot_fit_curve(a,*ndeg+1,compos[0].day,compos[ncompos-1].day))
	no_error = 0;

    /*
     * Report reduced chisq
     */

    if(no_error)
      if((chisq = chisq_fit(compos,ncompos,a,*ndeg+1,0,1)) < 0.0)
	no_error = 0;

    /*
     * Continue?
     */

    printf("Another fit? (y/n) [y] ");
    fgets(line,MAX,stdin);
    if(strcmp(line,"") != 0) {
      if(line[0] == 'N' || line[0] == 'n')
	contin = 0;
    }
    else
      a = del_array(a);
  }

  /*
   * Return
   */

  if(no_error)
    return a;
  else {
    fprintf(stderr,"ERROR: fit_compos\n");
    return del_array(a);
  }
}

/*.......................................................................
 *
 * Function plot_fit_curve
 *
 * Plots the polynomial that has been fit to the light curve
 *
 * Inputs: float *a            polynomial coefficients
 *         int ncoeff          number of coefficients
 *         float tmin          earliest time for which the fit is valid
 *         float tmax          latest time for which the fit is valid
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int plot_fit_curve(float *a, int ncoeff, float tmin, float tmax)
{
  int i,j;            /* Looping variables */
  int nstep=500;      /* Number of steps used to compute curve */
  float xtmp,ytmp;    /* Values of x,y used in constructing the curve */
  float *polyx=NULL;  /* Polynomial in x */

  cpgsci(2);
  cpgslw(5);

  /*
   * Allocate memory for container for polynomial values
   */

  if(!(polyx = new_array(ncoeff,1))) {
    fprintf(stderr,"ERROR: plot_fit_curve\n");
    return 1;
  }

  /*
   * Loop through range of x values in step sizes determined by nstep.
   * At each value of x, compute the value of y by first calling poly
   *  to compute the values of the various powers of x at that step, and
   *  then multiplying by the coefficients contained in a.
   */

  for(i=0; i<nstep; i++) {

    /*
     * Compute the values of x and y at this step
     */

    xtmp = tmin + i*(tmax - tmin)/(1.0 * nstep);
    sinpoly(xtmp,polyx-1,ncoeff);
    ytmp = 0;
    for(j=0; j<ncoeff; j++)
      ytmp += a[j] * polyx[j];

    /*
     * Now connect this point to the previous one with a cpgdraw call
     */

    if(i == 0)
      cpgmove(xtmp,ytmp);
    else
      cpgdraw(xtmp,ytmp);
  }

  /*
   * Clean up and exit
   */

  cpgslw(1);
  cpgsci(1);
  polyx = del_array(polyx);

  return 0;
}

/*.......................................................................
 *
 * Function plot_resid
 *
 * Plots the residuals between a light curve and the function that has been 
 *  fit to the curve.  Also calculates the RMS uncertainty in the mean value
 *  (which should be very close to 0.0).
 *
 * Inputs: Fluxrec *flux       input light curve
 *         int npoints         number of points in light curve
 *         float *a            coefficients of the fitted function
 *         int ncoeff          number of coefficients
 *         float *meanrms      rms uncertainty in mean (set by this function)
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int plot_resid(Fluxrec *flux, int npoints, float *a, int ncoeff,
	       float *meanrms)
{
  int i,j;              /* Looping variables */
  int no_error=1;       /* Flag set to 0 on error */
  float ytmp;           /* Value of fitted function at a given day */
  float mean;           /* Mean of residual curve */
  float *polyx=NULL;    /* Polynomial in x */
  Fluxrec *resid=NULL;  /* Residual curve */
  Fluxrec *rptr;        /* Pointer to navigate resid */
  Fluxrec *fptr;        /* Pointer to navigate flux */

  /*
   * Allocate memory
   */

  if(!(polyx = new_array(ncoeff,1))) {
    fprintf(stderr,"ERROR: plot_resid\n");
    return 1;
  }

  if(!(resid = new_fluxrec(npoints)))
    no_error = 0;

  /*
   * Loop through flux, at each step compute the value of the function
   *  and then subtract it from the flux.
   */

  if(no_error) {
    for(i=0,fptr=flux,rptr=resid; i<npoints; i++,fptr++,rptr++) {

    /*
     * Compute the function value at this step
     */

    sinpoly(fptr->day,polyx-1,ncoeff);
    ytmp = 0;
    for(j=0; j<ncoeff; j++)
      ytmp += a[j] * polyx[j];

    /*
     * Now set resid->flux to be the difference between flux->flux and the
     *  fitted function
     */

    *rptr = *fptr;
    rptr->flux -= ytmp;
    }
  }

  /*
   * Now plot the points
   */

  if(no_error)
    if(plot_lcurve(resid,npoints,"Days","Residual Flux Density",""))
      no_error = 0;

  /*
   * Draw a line at 0.0 residual
   */

  if(no_error) {
    cpgsci(2);
    cpgslw(5);
    cpgmove(-999,0.0);
    cpgdraw(999,0.0);
    cpgslw(1);
    cpgsci(1);
  }

  /*
   * Calculate mean and rms
   */

  if(no_error) {
    if(calc_mean(resid,npoints,&mean,meanrms))
      no_error = 0;
    else {
      *meanrms /= sqrt(npoints);
      printf("\nplot_resid: Mean of residuals = %8.4f. ",mean);
      printf("RMS uncertainty in mean = %6.4f\n",*meanrms);
    }
  }

  /*
   * Clean up and exit
   */

  cpgslw(1);
  cpgsci(1);
  polyx = del_array(polyx);
  resid = del_fluxrec(resid);

  if(no_error)
    return 0;
  else {
    fprintf(stderr,"ERROR: plot_resid\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function find_best_shift
 *
 * A curve and a function are input, along with the min and max times
 *  for which the function is valid.  The curve is shifted in time and
 *  at each step, compared to the function in the region in which they
 *  overlap.  The reduced chisq is calculated at each step, and the
 *  step giving the best fit is recorded.
 *
 * NB: It is assumed that the curve is sorted in time order.
 *
 * Inputs: Fluxrec *flux       light curve to be shifted
 *         int npoints         number of points in light curve
 *         float *a            array containing coefficients of fit function
 *         int ncoeff          number of coefficients
 *         float fitrms        typical RMS error associated with the fitting
 *         float stepsize      size of steps when shifting first curve
 *         float tmin          earliest time for which fit function is valid
 *         float tmax          latest time for which fit function is valid
 *         float *bestlag      lag with lowest chisq (set by this function)
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int find_best_shift(Fluxrec *flux, int npoints, float *a, int ncoeff, 
		    float fitrms, float stepsize, float tmin, float tmax,
		    float *bestlag)
{
  int i;                    /* Looping variable */
  int no_error=1;           /* Flag set to 0 on error */
  int noverlap;             /* Number of points in overlap region */
  int first=1;              /* Flag set to indicate first time through loop */
  float lag;                /* Current value of lag */
  float lagday;             /* flux1->day - lag */
  float chisq;              /* Reduced chisq */
  float bestchi;            /* Lowest reduced chisq */
  float maxshift;           /* Maximum shift for flux */
  float mean;               /* Mean value for flux */
  Fluxrec *shift=NULL;      /* Shifted version of flux */
  Fluxrec *sptr;            /* Pointer to navigate shift */
  Fluxrec *fptr;            /* Pointer to navigate flux */
  FILE *ofp=NULL;           /* Output file pointer */

  /*
   * Find the max shift for flux
   */

  maxshift = ((flux+npoints-1)->day - flux->day)/2.0;

  /*
   * Calculate the mean value of flux
   */

  if(calc_swmean(flux,npoints,&mean)) {
    fprintf(stderr,"ERROR: find_best_shift\n");
    return 1;
  }

  /*
   * Open up the output file for writing
   */

  if(!(ofp = fopen("fit_1_to_3.out","w"))) {
    fprintf(stderr,"ERROR: find_best_shift\n");
    return 1;
  }

  /*
   * Allocate memory for shift
   */

  if(!(shift = new_fluxrec(npoints)))
    no_error = 0;


  /*
   * Loop over possible shifts for flux
   */

  lag = *bestlag = -maxshift;
  while(lag < maxshift && no_error) {

    /*
     * Initialize
     */

    noverlap = 0;
    sptr = shift;

    /*
     * Loop through flux, adding to shift if the point is good
     *  (match > -1) and if the value with the lag subtracted is between
     *  tmin and tmax.  For the errors, first take the error associated
     *  with the point and divide by its mean value, then add error
     *  associated with the fitting function in quadrature.
     */

    for(i=0,fptr=flux; i<npoints; i++,fptr++) {
      lagday = fptr->day - lag;
      if(fptr->match > -1 && lagday >= tmin && lagday <=tmax) {
	*sptr = *fptr;
	sptr->day = lagday;
	sptr->flux /= mean;
	sptr->err /= mean;
	sptr->err = sqrt((sptr->err * sptr->err) + fitrms*fitrms);
	sptr++;
	noverlap++;
      }
    }

    /*
     * Calculate the reduced chisq by calling chisq_fit
     */

    if((chisq = chisq_fit(shift,noverlap,a,ncoeff,1,0)) < 0) {
      no_error = 0;
      break;
    }
    else
      fprintf(ofp,"%6.2f %f %d\n",lag,chisq,noverlap);

    /*
     * Compare to previous best value of chisq
     */

    if(first) {
      bestchi = chisq;
      first = 0;
    }
    else if(chisq < bestchi) {
      bestchi = chisq;
      *bestlag = lag;
    }

    /*
     * Increment lag
     */

    lag += stepsize;
  }

  /*
   * Print out results
   */

  printf("\nfind_best_shift.  Best reduced chisq of %f at a lag of %6.2f\n",
	 bestchi,*bestlag);

  /*
   * Clean up and exit
   */

  shift = del_fluxrec(shift);
  if(ofp)
    fclose(ofp);

  if(no_error)
    return 0;
  else {
    fprintf(stderr,"ERROR: find_best_shift.\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function plot_best_shift
 *
 * Plots free curve at its best-fit lag, along with the fitted function.
 *
 * Inputs: Fluxrec *flux       free light curve
 *         int npoints         number of points in curve
 *         float mean          mean value of free light curve
 *         float *a            coefficients of fitted function
 *         int ncoeff          number of coefficients
 *         float tmin          earliest time for which the fit is valid
 *         float tmax          latest time for which the fit is valid
 *         float bestlag       best-fit lag
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int plot_best_shift(Fluxrec *flux, int npoints, float mean, float *a, 
		    int ncoeff, float tmin, float tmax, float bestlag)
{
  int i;                  /* Looping variable */
  int no_error=1;         /* Flag set to 0 on error */
  int ngood=0;            /* Number of good points in flux */
  float lagday;           /* Shifted version of flux->day */
  Fluxrec *shift=NULL;    /* Shifted version of flux */
  Fluxrec *sptr;          /* Pointer to navigate shift */
  Fluxrec *fptr;          /* Pointer to navigate flux */

  /*
   * Allocate memory for shift
   */

  if(!(shift = new_fluxrec(npoints))) {
    fprintf(stderr,"ERROR: plot_best_shift\n");
    return 1;
  }

  /*
   * Create shifted curve
   */

  sptr = shift;
  for(i=0,fptr=flux; i<npoints; i++,fptr++) {
    lagday = fptr->day - bestlag;
    if(fptr->match > -1) {
      *sptr = *fptr;
      sptr->day = lagday;
      sptr->flux /= mean;
      sptr->err /= mean;
      sptr++;
      ngood++;
    }
  }

  /*
   * Plot shifted curve
   */

  if(no_error)
    if(plot_lcurve(shift,ngood,"Days","Flux Density",""))
      no_error = 0;

  /*
   * Plot fitted function
   */
  
  if(no_error)
    if(plot_fit_curve(a,ncoeff,tmin,tmax))
      no_error = 0;

  /*
   * Clean up and return
   */

  shift = del_fluxrec(shift);
  return 0;
}


/*.......................................................................
 *
 * Weird function to make pgplot work
 *
 */

int MAIN_(void)
{
  
  return 0;
}
