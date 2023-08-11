/*
 * correlate.c
 *
 * A library of functions to calculate various types of correlations
 *  between light curves.
 *
 * 25Sep98 CDF,  Moved functions over from lc_funcs.c.  For history prior
 *                to this date, see lc_funcs.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "structdef.h"
#include "nr.h"
#include "nrutil.h"
#include "lc_funcs.h"
#include "lc_interp.h"
#include "correlate.h"

/*.......................................................................
 *
 * Function call_xcorr_fft
 *
 * Does a series of cross-correlations between pairs of light curves, by
 *  calling the function do_corr for each pair of curves.  Puts the
 *  output into a file designated by filename.  This function manipulates
 *  the input curves in a way appropriate for the FFT based 
 *  cross-correlation functions (cross_corr_fft and cross_corr_nr).
 *  The related function call_xcorr_time sets things up for the time-domain
 *  based cross-correlation (cross_corr_time).
 *
 * Inputs: Fluxrec *flux[]     light curves
 *         int size            number of points in each curve
 *         char *filename      name of output file
 *         int nbad            number of bad data points
 *         FILE *logfp         logfile pointer
 *         int doprint         flag set to 1 if output files are desired
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> failure
 *
 * v04Jul98 CDF, Modified logfile output to have all lags on one line
 *                rather than having separate output lines for the
 *                BA, BC, and BD lags.
 * v09Jul98 CDF, Added doprint flag as a passed parameter to control whether
 *                output gets printed.
 * v21Jul98 CDF, Changed to call_xcorr_fft to set things up for FFT-based
 *                cross-correlation functions.  The separate setup for
 *                time-domain-based cross-correlations has been split off
 *                into call_xcorr_time.
 *
 */

int call_xcorr_fft(Fluxrec *flux[], int size, char *filename, int nbad, 
		   FILE *logfp, int doprint)
{
  int i;                       /* Looping variable */
  int no_error=1;              /* Flag set to 0 on error */
  int corsize;                 /* Length of correlation curves */
  float meanba,meanbc,meanbd;  /* Mean values of cross-correlation curves */
  float rmsba,rmsbc,rmsbd;     /* RMS values of cross-correlation curves */
  float norm=1.0;              /* Normalization for curves */
  float normday;               /* Day from which normalization value comes */
  Fluxrec best_ba;             /* Highest cross-correlation for B-A */
  Fluxrec best_bc;             /* Highest cross-correlation for B-C */
  Fluxrec best_bd;             /* Highest cross-correlation for B-D */
  Fluxrec *zmean[N08]={NULL};  /* Zero-mean light curves */
  Fluxrec *acorrb=NULL;        /* B autocorrelation */
  Fluxrec *corrbc=NULL;        /* B-C cross-correlation */
  Fluxrec *corrba=NULL;        /* B-A cross-correlation */
  Fluxrec *corrbd=NULL;        /* B-D cross-correlation */
  Fluxrec *cp1,*cp2,*cp3;      /* Pointers to navigate cross-corr arrays */
  FILE *ofp=NULL;              /* Output file */

  /*
   * Initialize the best lag containers to something ridiculous.  In this
   *  function the lag will be contained in the day member of the Fluxrec
   *  structure and the cross-correlation coefficient will be contained 
   *  in the flux member.
   */

  best_ba.flux = -999.0;
  best_bc.flux = -999.0;
  best_bd.flux = -999.0;

  /*
   * Normalize the curves since we are looking for correlations in
   *  fractional variations and then subtract 1 to create a zero-mean
   *  light curve.  This helps by getting rid of the DC component in
   *  the cross-correlation.
   */

  printf("\nCalculating cross-correlations....\n\n");
  for(i=0; i<N08; i++) {
    if(no_error)
      if(!(zmean[i] = norm_zero_mean(flux[i],size)))
	no_error = 0;
  }

  /*
   * First do an autocorrelation on curve B and use the peak
   *  value (at zero lag) to normalize the cross-correlation
   *  curves.  
   */

  if(no_error)
    if(!(acorrb = do_corr(zmean[1],zmean[1],size,&corsize,nbad)))
      no_error = 0;

  if(no_error) {
    norm = -9999.0;
    for(i=0,cp1=acorrb; i<corsize; i++,cp1++) {
      if(cp1->flux > norm) {
	norm = cp1->flux;
	normday = cp1->day;
      }
    }
    printf("\ncall_xcorr_fft: Normalization = %f from day %5.1f\n\n",
	   norm,normday);
  }

  /*
   * Call the cross-correlation routines for each pair of
   *  zero-mean light curves.
   */

  if(no_error)
    if(!(corrba = do_corr(zmean[1],zmean[0],size,&corsize,nbad)))
      no_error = 0;

  if(no_error)
    if(!(corrbc = do_corr(zmean[1],zmean[2],size,&corsize,nbad)))
      no_error = 0;

  if(no_error)
    if(!(corrbd = do_corr(zmean[1],zmean[3],size,&corsize,nbad)))
      no_error = 0;

  /*
   * Open output file
   */

  if(doprint && no_error)
    if((ofp = fopen(filename,"w")) == NULL) {
      fprintf(stderr,"ERROR: call_xcorr_fft.  Cannot open %s\n",filename);
      no_error = 0;
    }

  /*
   * Print cross-correlation curves to output files and also find
   *  lags with highest cross-correlation values.
   */

  if(no_error) {
    if(doprint) {
      fprintf(ofp,"#  Lag    Corr_ba   Corr_bc   Corr_bd\n");
      fprintf(ofp,"#------- --------- --------- ---------\n");
    }
    for(i=0,cp1=corrba,cp2=corrbc,cp3=corrbd; i<corsize; 
	i++,cp1++,cp2++,cp3++) {
      cp1->flux /= norm;
      cp2->flux /= norm;
      cp3->flux /= norm;
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
   * Find rms (and mean) of cross-correlation curves.  This can be used
   *  to assess the significance of the highest peaks in the curves.
   */

  if(no_error) {
    calc_mean(corrba,corsize,&meanba,&rmsba);
    calc_mean(corrbc,corsize,&meanbc,&rmsbc);
    calc_mean(corrbd,corsize,&meanbd,&rmsbd);
#if 0
    printf(" call_xcorr_fft: B-A xcorr curve mean = %f, rms = %f\n",
	   meanba,rmsba);
    printf(" call_xcorr_fft: B-C xcorr curve mean = %f, rms = %f\n",
	   meanbc,rmsbc);
    printf(" call_xcorr_fft: B-D xcorr curve mean = %f, rms = %f\n\n",
	   meanbd,rmsbd);
#endif
  }

  /*
   * Print best lags to logfile, if it exists.  Take the negative of
   *  best_b?.day because if B leads the other light curve, the
   *  cross-correlation will return a negative number for the best-fit
   *  lags.  Since we know/hope that B leads all other values and we
   *  are looking for the amount that the other curves lag B, just
   *  take the negative of best_b?.day.
   */

  if(logfp && no_error) {
    fprintf(logfp,"# Cross correlation results\n");
    fprintf(logfp,"%7.2f %f %f ",-best_ba.day,best_ba.flux,rmsba);
    fprintf(logfp,"%7.2f %f %f ",-best_bc.day,best_bc.flux,rmsbc);
    fprintf(logfp,"%7.2f %f %f\n",-best_bd.day,best_bd.flux,rmsbd);
  }

  /*
   * Also print output to screen
   */

  if(no_error) {
    printf(" call_xcorr_fff: B-A best lag = %7.2f d (%f = %6.2f sigma)\n",
	   -best_ba.day,best_ba.flux,(best_ba.flux / rmsba));
    printf(" call_xcorr_fft: B-C best lag = %7.2f d (%f = %6.2f sigma)\n",
	   -best_bc.day,best_bc.flux,(best_bc.flux / rmsbc));
    printf(" call_xcorr_fft: B-D best lag = %7.2f d (%f = %6.2f sigma)\n",
	   -best_bd.day,best_bd.flux,(best_bd.flux / rmsbd));
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
    fprintf(stderr,"ERROR: call_xcorr_fft\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function do_corr
 *
 * Computes a correlation using an FFT.  If the data is not an integral
 *  power of 2, the data is gridded onto a grid with the nearest power of
 *  2 pts.  If the data does have an integral power of 2 number of points,
 *  no gridding is done.  This function calls the designated
 *  cross-correlation function.
 *
 * Inputs: Fluxrec *flux1      first light curve
 *         Fluxrec *flux2      second light curve
 *         int size            size of Fluxrec arrays
 *         int *corsize        size of correlation array (set by function)
 *         int nbad            number of bad data points
 *
 * Output: Fluxrec *correl     cross-correlation of the arrays
 *
 * v24Jun98 CDF, Moved normalization and creation of zero-mean curve
 *                up one level to call_xcorr.
 */

Fluxrec *do_corr(Fluxrec *flux1, Fluxrec *flux2, int size, int *corsize,
		 int nbad)
{
  int no_error=1;         /* Flag set to 0 on error */
  int nx;                 /* Size of gridded arrays */
  float xmin,xmax;        /* Min and max day numbers */
  float dx;               /* Grid step-size between days */
  float test;             /* Tests size to see if it's a power of 2 */
  float *zpad1=NULL;      /* Zero padded light curve */
  float *zpad2=NULL;      /* Zero padded light curve */
  Fluxrec *grflux1=NULL;  /* Gridded version of flux1 */
  Fluxrec *grflux2=NULL;  /* Gridded version of flux2 */
  Fluxrec *correl=NULL;   /* Cross-correlation of flux1 and flux2 */

  /*
   * Note that the day part of the flux array will be sorted by 
   *  construction, so the min and max values will be easy to find.
   */

  xmin = flux1->day;
  xmax = (flux1+size-1)->day;

  /*
   * Check if size is a power of 2 because the number of points in the
   *  curve MUST be a power of two if the cross-correlation function
   *  is cross_corr_fft or cross_corr_nr.
   */

  test = log((float)size)/log(2.0);

  /*
   * If it is or if CORRFUNC == 3 , then set the values for nx and dx, 
   *  as simple functions of size, xmin, and xmax, and then zero
   *  pad the input data
   */

  if(test-(int)test == 0) {
    nx = size;
    dx = (xmax-xmin)/(nx-1);

    if(!(zpad1 = zero_pad(flux1,nx)))
      no_error = 0;
    if(!(zpad2 = zero_pad(flux2,nx)))
      no_error = 0;
  }    


  /*
   * If size is not a power of 2 and CORRFUNC != 3, interpolate the data onto 
   *  an appropriately spaced grid with a size that is a power of 2.
   */

  else {

    /*
     * Compute new values for nx and dx
     */

    nx = pow(2.0f,floor((log((double)size)/log(2.0)+0.5)));
    dx = (xmax-xmin)/(nx-1);
    printf("do_corr: Using grid of %d points for cross-correlation.\n",nx);
    if(!(grflux1 = lin_interp(flux1,size,nx,dx,nbad)))
      no_error = 0;
      
    if(no_error)
      if(!(grflux2 = lin_interp(flux2,size,nx,dx,nbad)))
	no_error = 0;
      
    /*
     * Zero-pad the gridded data
     */
      
    if(no_error)
      if(!(zpad1 = zero_pad(grflux1,nx)))
	no_error = 0;
      
    if(no_error)
      if(!(zpad2 = zero_pad(grflux2,nx)))
	no_error = 0;
  }

  /*
   * Do the cross-correlation, via one of three styles:
   *  1. FFT method (Numerical Recipes FFTs, then mulitply and transform back
   *  2. Numerical Recipes correl function
   *  3. Time domain method with no FFTs
   */

  if(no_error) {
    nx *= ZEROFAC;
    switch(CORRFUNC) {
    case 1:
      if(!(correl = cross_corr_fft(zpad1,zpad2,nx,dx)))
	no_error = 0;
      break;
    case 2:
      if(!(correl = cross_corr_nr(zpad1,zpad2,nx,dx)))
	no_error = 0;
      break;
    default:
      fprintf(stderr,"ERROR: do_corr.  Bad choice of correlation function\n");
      no_error = 0;
    }
  }

  /*
   * Clean up
   */

  zpad1 = del_array(zpad1);
  zpad2 = del_array(zpad2);
  grflux1 = del_fluxrec(grflux1);
  grflux2 = del_fluxrec(grflux2);

  if(no_error) {
    *corsize = nx;
    return correl;
  }
  else {
    *corsize = 0;
    return NULL;
    fprintf(stderr,"ERROR: do_corr\n");
  }
}

/*.......................................................................
 *
 * Function cross_corr_fft
 *
 * Does the cross correlation of two Fluxrec arrays.  The arrays will
 *  have been already gridded by do_corr
 *
 * Inputs: float *flux1        first flux curve
 *         float *flux2        second flux curve
 *         int size            size of arrays
 *         float dx            step size for lag axis
 *
 * Output: Fluxrec *correl     correlation of the two arrays
 *
 * v05Jul97 CDF, Added zero padding
 * v22Jun98 CDF, Added normalization
 * v24Jun98 CDF, Moved zero-padding to a separate function and put the
 *                function call in do_corr.  Moved normalization up
 *                to call_xcorr.
 */

Fluxrec *cross_corr_fft(float *flux1, float *flux2, int size, float dx)
{
  int i;                  /* Looping variable */
  int no_error=1;         /* Flag set to 0 on error */
  float re1,im1;          /* Real and imag. parts of FFT of flux1 */
  float re2,im2;          /* Real and imag. parts of FFT of flux2 */
  float nfact;            /* Normalization factor for FFTs */
  float *y1=NULL;         /* Fluxes from first array -- complex format */
  float *y2=NULL;         /* Fluxes from second array -- complex format */
  float *yptr1,*yptr2;    /* Pointers for y arrays */
  float *fptr1,*fptr2;    /* Pointers to navigate flux1 and flux2 */
  Fluxrec *correl=NULL;   /* Cross-correlation of flux1 and flux2 */
  Fluxrec *cptr;          /* Pointer to navigate correl */

  /*
   * Allocate memory for float arrays
   */

  if(!(y1 = new_array(2*size,1))) {
    return NULL;
  }

  if(!(y2 = new_array(2*size,1)))
    no_error = 0;

  /*
   * Allocate memory for the cross-correlation array
   */

  if(no_error) {
    if(!(correl = new_fluxrec(size))) {
      no_error = 0;
    }
  }

  /*
   * Put information from Fluxrec arrays into float arrays, which
   *  are in the Numerical Recipes "complex" format which consists.
   *  of alternating real and complex values.  Since our arrays
   *  are real, set all the imaginary members to zero.  Note that
   *  new_array automatically sets all of its members to zero, so
   *  all we have to do is fill in the real members of the arrays with
   *  the members of flux1 and flux2.
   */

  yptr1 = y1;
  yptr2 = y2;

  if(no_error) {
#if 0
    for(i=0,fptr1=flux1,fptr2=flux2; i<size; i++,fptr1++,fptr2++) {
      *yptr1 = *fptr1;
      *yptr2 = *fptr2;
      yptr1 += 2;
      yptr2 += 2;
    }
#endif
    for(i=0; i<size; i++) {
      y1[2*i] = flux1[i];
      y2[2*i] = flux2[i];
    }
  }

  /*
   * Do the cross-correlation
   */

  if(no_error) {

    /*
     * Fourier transform the two arrays
     */

    four1(y1-1,size,1);
    four1(y2-1,size,1);

    /*
     * Now we have the transform -> compute the correlation.
     */

    nfact = 1.0;

    for(i=0; i<size; i++) {
      re1 = y1[2*i];
      im1 = y1[2*i+1];
      re2 = y2[2*i];
      im2 = y2[2*i+1];

      /*
       * Store the correlation in the first complex array.
       */

      y1[2*i] = (re1*re2 + im1*im2)/nfact;
      y1[2*i+1] = (im1*re2 - im2*re1)/nfact;
    }

    /*
     * Transform back.
     */

    four1(y1-1,size,-1);
  }

  /*
   * Put the results into the cross correlation array
   */

  if(no_error) {
#if 0
    yptr1 = y1;
    for(i=0,cptr=correl; i<size/2; i++,cptr++) {
      yptr2 = yptr1 + size;
      cptr->day = (i-size/2)*dx;
      (cptr+size/2)->day = i*dx;
      cptr->flux = *yptr2;
      (cptr+size/2)->flux = *yptr1;
      yptr1 += 2;
    }
#endif
    for(i=0; i<size/2; i++) {
      correl[i].day = (i-size/2)*dx;
      correl[i+size/2].day = i*dx;
      correl[i].flux = y1[size+2*i];
      correl[i+size/2].flux = y1[2*i];
    }
  }

  /*
   * Clean up
   */

  y1 = del_array(y1);
  y2 = del_array(y2);

  if(no_error)
    return correl;
  else {
    fprintf(stderr,"ERROR: cross_corr_fft\n");
    return NULL;
  }
}

/*.......................................................................
 *
 * Function cross_corr_nr
 *
 * Cross correlates two light curves using the Numerical Recipes
 *  correlation routine.
 *
 * Inputs: float *flux1        first light curve
 *         float *flux2        second light curve
 *         int npoints         number of points in light curves
 *         float dx            step size for lag axis
 *
 * Output: Fluxrec *xcorr      cross-correlation curve
 *
 */

Fluxrec *cross_corr_nr(float *flux1, float *flux2, int npoints, float dx)
{
  int i;                    /* Looping variable */
  int no_error=1;           /* Flag set to 0 on error */
  float *y1=NULL;           /* Light curve 1 fluxes */
  float *y2=NULL;           /* Light curve 2 fluxes */
  float *corrout=NULL;      /* Float output from correlation */
  float *yptr1,*yptr2;      /* Pointers to navigate y1, y2, and corrout */
  float *fptr1,*fptr2;      /* Pointers to navigate flux1 and flux2 */
  Fluxrec *xcorr=NULL;      /* Cross-correlation curve */
  Fluxrec *xptr;            /* Pointer to navigate xcorr */

  /*
   * Allocate memory for the float versions of the light curves
   */

  if(!(y1 = new_array(npoints,1))) {
    fprintf(stderr,"ERROR: cross_corr_nr\n");
    return NULL;
  }

  if(!(y2 = new_array(npoints,1)))
    no_error = 0;

  /*
   * Transfer data into y1 and y2.  We need to do this because the
   *  correl routine modifies the arrays passed to i.
   */

  if(no_error) {
    for(i=0,fptr1=flux1,fptr2=flux2,yptr1=y1,yptr2=y2; i<npoints;
	i++,fptr1++,fptr2++,yptr1++,yptr2++) {
      *yptr1 = *fptr1;
      *yptr2 = *fptr2;
    }
  }

  /*
   * Allocate memory for the output from the cross-correlation
   */

  if(no_error)
    if(!(corrout = new_array(npoints,2)))
      no_error = 0;

  /*
   * Call the correlation function.
   * NB: The -1 after y1, y2, and corrout is necessary because of the
   *  Numerical Recipes fortran convention.
   */

  if(no_error)
    correl(y1-1,y2-1,(unsigned long) npoints,corrout-1);

  /*
   * Allocate memory for output Fluxrec array
   */

  if(no_error)
    if(!(xcorr = new_fluxrec(npoints)))
      no_error = 0;

  /*
   * Fill array.  
   * *** NB: The output from correl is the cross-correlation curve with
   *         the positive values going from 0 to (npoints/2-1) and
   *         the negative values going from npoints to npoints/2.
   *         We want the array to be in the "correct" order with the most
   *         negative value at position 0 and monotonically increasing from
   *         there.  Therefore fill the array accordingly.
   */

  if(no_error) {
    for(i=0,yptr1=corrout,xptr=xcorr; i<npoints/2; i++,yptr1++,xptr++) {
      yptr2 = yptr1 + npoints/2;
      xptr->day = (i-npoints/2)*dx;
      xptr->flux = *yptr2;
      (xptr+npoints/2)->day = i*dx;
      (xptr+npoints/2)->flux = *yptr1;
    }
  }

  /*
   * Clean up and exit
   */

  y1 = del_array(y1);
  y2 = del_array(y2);
  corrout = del_array(corrout);

  if(no_error)
    return xcorr;
  else {
    fprintf(stderr,"ERROR: cross_corr_nr.\n");
    return NULL;
  }
}

/*.......................................................................
 *
 * Function call_xcorr_time
 *
 * Does a series of cross-correlations between pairs of light curves, by
 *  calling the function do_corr for each pair of curves.  Puts the
 *  output into a file designated by filename.  This function does the
 *  setup for time-domain-based cross-correlations.  It previously was
 *  part of the more general call_xcorr (now called call_xcorr_fft).
 *
 * Inputs: Fluxrec *raw[]      raw light curves
 *         Fluxrec *interp[]   interpolated light curves
 *         int nraw            number of points in raw curves
 *         int ninterp         number of points in interpolated curves
 *         Setup *setup        container for smoothing/interpolation info
 *         FILE *logfp         logfile pointer
 *         char *filename      name of output file
 *         int doprint         flag set to 1 if output files are desired
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> failure
 *
 * v22Aug98 CDF, Split off from call_xcorr (which now is called call_xcorr_fft
 *                and sets up for FFT-based cross-correlations).
 * v07Sep98 CDF, Now print out number of overlapping points, too.
 * v12Sep98 CDF, Now pass Setup container to this function, so that the
 *                type of smoothing/interpolation is known.
 * v24Sep98 CDF, Add containers for lowest Bevington probability, if
 *                calculated, because the lowest probability lag may not
 *                match the highest correlation coefficient.
 *
 */

int call_xcorr_time(Fluxrec *raw[], Fluxrec *interp[], int nraw, int ninterp, 
		    Setup *setup, FILE *logfp, char *filename, int doprint)
{
  int i;                       /* Looping variable */
  int no_error=1;              /* Flag set to 0 on error */
  float dx;                    /* Step size in interpolated curves */
  Fluxrec best_ba;             /* Highest cross-correlation for B-A */
  Fluxrec best_bc;             /* Highest cross-correlation for B-C */
  Fluxrec best_bd;             /* Highest cross-correlation for B-D */
  Fluxrec bestbev_ba;          /* Highest cross-correlation for B-A */
  Fluxrec bestbev_bc;          /* Highest cross-correlation for B-C */
  Fluxrec bestbev_bd;          /* Highest cross-correlation for B-D */
  Fluxrec *zmraw[N08]={NULL};  /* Zero-mean raw light curves */
  Fluxrec *zmint[N08]={NULL};  /* Zero-mean interpolated light curves */
  Fluxrec *corrbc=NULL;        /* B-C cross-correlation */
  Fluxrec *corrba=NULL;        /* B-A cross-correlation */
  Fluxrec *corrbd=NULL;        /* B-D cross-correlation */
  Fluxrec *cp1,*cp2,*cp3;      /* Pointers to navigate cross-corr arrays */
  FILE *ofp=NULL;              /* Output file */

  /*
   * Initialize the best lag containers to something ridiculous.  In this
   *  function the lag will be contained in the day member of the Fluxrec
   *  structure, the cross-correlation coefficient will be contained 
   *  in the flux member, and the Bevington probability, if computed,
   *  will be contained in the err member.
   */

  best_ba.flux = best_bc.flux = best_bd.flux = -999.0;
  bestbev_ba.err = bestbev_bc.err = bestbev_bd.err = 1.0;

  /*
   * Normalize the curves since we are looking for correlations in
   *  fractional variations and then subtract 1 to create a zero-mean
   *  light curve.  This helps by getting rid of the DC component in
   *  the cross-correlation.
   */

  for(i=0; i<N08; i++) {
    if(!(zmraw[i] = norm_zero_mean(raw[i],nraw)))
	no_error = 0;
    if(!(zmint[i] = norm_zero_mean(interp[i],ninterp)))
	no_error = 0;
  }

  /*
   * Calculate dx
   */

  if(no_error) {
    cp1 = interp[0];
    cp1++;
    dx = cp1->day - interp[0]->day;
    printf(" call_xcorr_time: dx = %7.3f\n",dx);
  }

  /*
   * Call the cross-correlation routines for each pair of
   *  zero-mean light curves.
   */

  if(no_error)
    if(!(corrba = cross_corr_time(zmraw[1],zmint[0],nraw,ninterp,dx,setup)))
      no_error = 0;

  if(no_error)
    if(!(corrbc =  cross_corr_time(zmraw[1],zmint[2],nraw,ninterp,dx,setup)))
      no_error = 0;

  if(no_error)
    if(!(corrbd =  cross_corr_time(zmraw[1],zmint[3],nraw,ninterp,dx,setup)))
      no_error = 0;

  /*
   * Open output file
   */

  if(doprint && no_error)
    if((ofp = fopen(filename,"a")) == NULL) {
      fprintf(stderr,"ERROR: call_xcorr_time.  Cannot open %s\n",filename);
      no_error = 0;
    }

  /*
   * Print cross-correlation curves to output files and also find
   *  lags with highest cross-correlation values.
   */

  if(no_error) {
    if(doprint) {
      fprintf(ofp,"#  Lag    Corr_ba   Corr_bc   Corr_bd  ");
      fprintf(ofp,"  P(ba)     P(bc)     P(bd)   N_overlap\n");
      fprintf(ofp,"#------- --------- --------- --------- ");
      fprintf(ofp,"--------- --------- --------- ---------\n");
    }
    for(i=0,cp1=corrba,cp2=corrbc,cp3=corrbd; i<ninterp; 
	i++,cp1++,cp2++,cp3++) {
      if(doprint)
	fprintf(ofp,"%8.3f %9f %9f %9f %6.3e %6.3e %6.3e %6d\n",
		-cp1->day,cp1->flux,cp2->flux,cp3->flux,
		cp1->err,cp2->err,cp3->err,cp1->match);
      if(cp1->flux > best_ba.flux)
	best_ba = *cp1;
      if(cp2->flux > best_bc.flux)
	best_bc = *cp2;
      if(cp3->flux > best_bd.flux)
	best_bd = *cp3;
      if(cp1->err < bestbev_ba.err)
	bestbev_ba = *cp1;
      if(cp2->err < bestbev_bc.err)
	bestbev_bc = *cp2;
      if(cp3->err < bestbev_bd.err)
	bestbev_bd = *cp3;
    }
  }

  /*
   * Print best lags to logfile, if it exists.  Take the negative of
   *  best_b?.day because if B leads the other light curve, the
   *  cross-correlation will return a negative number for the best-fit
   *  lags.  Since we know/hope that B leads all other values and we
   *  are looking for the amount that the other curves lag B, just
   *  take the negative of best_b?.day.
   */

  if(logfp && no_error) {
    fprintf(logfp,"# Cross correlation results\n");
    fprintf(logfp,"%7.2f %f 0.0 ",-best_ba.day,best_ba.flux);
    fprintf(logfp,"%7.2f %f 0.0 ",-best_bc.day,best_bc.flux);
    fprintf(logfp,"%7.2f %f 0.0 ",-best_bd.day,best_bd.flux);
    fprintf(logfp,"%7.2f %f 0.0 ",-bestbev_ba.day,log10(bestbev_ba.err));
    fprintf(logfp,"%7.2f %f 0.0 ",-bestbev_bc.day,log10(bestbev_bc.err));
    fprintf(logfp,"%7.2f %f 0.0\n",-bestbev_bd.day,log10(bestbev_bd.err));
  }

  /*
   * Also print output to screen
   */

  if(no_error) {
    printf(" call_xcorr_time: B-A best lag = %7.2f d (%f)\n",
	   -best_ba.day,best_ba.flux);
    printf(" call_xcorr_time: B-C best lag = %7.2f d (%f)\n",
	   -best_bc.day,best_bc.flux);
    printf(" call_xcorr_time: B-D best lag = %7.2f d (%f)\n",
	   -best_bd.day,best_bd.flux);
    printf(" call_xcorr_time: B-A best Bevington lag = %7.2f d (%f)\n",
	   -bestbev_ba.day,bestbev_ba.err);
    printf(" call_xcorr_time: B-C best Bevington lag = %7.2f d (%f)\n",
	   -bestbev_bc.day,bestbev_bc.err);
    printf(" call_xcorr_time: B-D best Bevington lag = %7.2f d (%f)\n",
	   -bestbev_bd.day,bestbev_bd.err);
  }

  /*
   * Clean up
   */

  for(i=0; i<N08; i++) {
    zmraw[i] = del_fluxrec(zmraw[i]);
    zmint[i] = del_fluxrec(zmint[i]);
  }
  corrbc = del_fluxrec(corrbc);
  corrba = del_fluxrec(corrba);
  corrbd = del_fluxrec(corrbd);
  if(ofp)
    fclose(ofp);

  if(no_error)
    return 0;
  else {
    fprintf(stderr,"ERROR: call_xcorr_time\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function cross_corr_time
 *
 * Cross correlates two light curves in the time domain -- without using any 
 *  FFTs.  it is possible that one of the light curves is unevenly sampled 
 *  while the other curve is interpolated onto a regular grid.  In this
 *  function, flux1 will be the raw curve (unevenly sampled) and flux2
 *  will be the interpolated curve.
 *
 * Inputs: Fluxrec *flux1      first light curve (raw)
 *         Fluxrec *flux2      second light curve (interpolated)
 *         int n1              number of points in first light curve
 *         int n2              number of points in second light curve
 *         float dx            step size for lag axis
 *         Setup *setup        container for smoothing info.
 *
 * Output: Fluxrec *xcorr      cross-correlation curve
 *
 * v22Aug98 CDF, Completely rewritten to use White and Peterson method
 *                for calculating CCF.
 * v07Sep98 CDF, Moved calculation of correlation coefficient into
 *                a separate function to give the opportunity of using
 *                either Bevington (corr_bevington) or White and Peterson
 *                (corr_white_peterson) methods.
 *               Record number of overlapping points in match member of xcorr.
 * v12Sep98 CDF, Setup container now passed to function so it can be used
 *                in corr_bevington.
 */

Fluxrec *cross_corr_time(Fluxrec *flux1, Fluxrec *flux2, int n1, int n2, 
			 float dx, Setup *setup)
{
  int i,j,k;                /* Looping variables */
  int no_error=1;           /* Flag set to 0 on error */
  int nprod;                /* Number of products computed in each sum */
  int nmeas;                /* Number of products using measured points */
  int nind;                 /* Number of independent points */
  int minlag,maxlag;        /* Limits on lags used for xcorr calculations */
  Fluxrec *fptr1,*fptr2;    /* Pointers to navigate flux1 and flux2*/
  Fluxrec *xcorr=NULL;      /* Cross-correlation curve */
  Fluxrec *xptr;            /* Pointer to navigate xcorr */
  Fluxrec *tmp1=NULL;       /* Container for overlapping portion of flux1 */
  Fluxrec *tmp2=NULL;       /* Container for overlapping portion of flux1 */
  Fluxrec *tp1,*tp2;        /* Pointers to navigate tmp1 and tmp2 */

  /*
   * Allocate memory for output Fluxrec array
   */

  if(!(xcorr = new_fluxrec(n2)))
    no_error = 0;

  /*
   * Allocate memory for tmp1 and tmp2.  Since the maximum overlap will
   *  be n2, just give them a size of n2 and then use nprod
   *  to inform corr_sum how big the arrays are for a given lag.
   */

  if(no_error)
    if(!(tmp1 = new_fluxrec(n2)))
      no_error = 0;

  if(no_error)
    if(!(tmp2 = new_fluxrec(n2)))
      no_error = 0;

  /*
   * Set min and max lags
   */

  if(1.0*n2/2.0 == (int) (1.0*n2/2.0)) {
    minlag = -n2/2;
    maxlag = n2/2;
  }
  else {
    minlag = -(n2 - 1)/2;
    maxlag = 1 + (n2 - 1)/2;
  }

  /*
   * Step through lags going from -(n2/2)*dx to [(n2/2)-1]*dx.
   * At each lag i, step through flux1 (incrementing j) and compute
   *  sums of flux1[j] * flux2[k] if where k is set such that
   *  flux1[j].day - lag = flux2[k].day.
   */

  if(no_error) {
    for(i=minlag,xptr=xcorr; i<maxlag; i++,xptr++) {

      /*
       * Initialize the values of xcorr, nprod, tp1, and tp2  at this lag
       */

      xptr->day = i*dx;
      xptr->flux = 0.0;
      nprod = 0;
      nmeas = 0;
      tp1 = tmp1;
      tp2 = tmp2;

      /*
       * Now step through flux1, putting the overlapping points of flux1
       *  and flux2 (with the given lag) into tmp1 and tmp2, as long as
       *  the flux1 point is good (match > -1)
       */

      for(j=0,fptr1=flux1; j<n1; j++,fptr1++) {
	if(fptr1->match >= 0) {
	  for(k=0,fptr2=flux2; k<n2; k++,fptr2++) {
	    if(fabs(fptr1->day - xptr->day - fptr2->day) < 0.2*dx) {
	      *tp1 = *fptr1;
	      *tp2 = *fptr2;
	      tp1++;
	      tp2++;
	      nprod++;
	      if(fptr1->match == 1)
		nmeas++;
	      break;
	    }
	  }
	}
      }

      /*
       * Put the number of overlapping points into the match member of
       *  xcorr.
       */

      xptr->match = nprod;

      /*
       * Estimate the number of independent points using the
       *  information in the Setup container.
       */

      if(setup) {
	switch(setup->smtype) {
	case SMBOXCAR: case SMTRIANGLE:
	  nind = floor(WINFRAC * ((tmp2+nprod-1)->day - tmp2->day) / 
		       setup->smwidth);
	  break;
	case SMVARBOX: case SMVARTRI:
	  nind = floor(nmeas / setup->nvar);
	  break;
	case SMGAUSS:
	  nind = floor(WINFRAC * ((tmp2+nprod-1)->day - tmp2->day) / 
		       (2.0 * setup->smwidth));
	  break;
	default:
	  nind = nprod;
	}
      }
      else {
	nind = nprod;
      }

      /*
       * Now call one of the functions that calculates the linear
       *  correlation coefficient for the overlapping points of the
       *  two curves.
       */
#if 1
      if(corr_bevington(tmp1,tmp2,nprod,xptr,nind))
	no_error = 0;
#else
      if(corr_white_peterson(tmp1,tmp2,nprod,&xptr->flux))
	no_error = 0;
#endif
    }
  }

  /*
   * Clean up and exit
   */

  tmp1 = del_fluxrec(tmp1);
  tmp2 = del_fluxrec(tmp2);

  if(no_error)
    return xcorr;
  else {
    fprintf(stderr,"ERROR: cross_corr_time.\n");
    return del_fluxrec(xcorr);
  }
}

/*.......................................................................
 *
 * Function corr_bevington
 *
 * Calculates a correlation coefficient according to method in Bevington
 *  where
 *            s_jk^2
 *    r_jk = --------
 *           s_j * s_k
 *
 *  and the covariance, variance and mean are defined in Book 6, p. 75.
 *
 * Inputs: Fluxrec *flux1      first flux (x_j)
 *         Fluxrec *flux2      second flux (x_k)
 *         int npoints         number of points in each curve (N)
 *         Fluxrec *r          linear correlation coefficient (set by this
 *                              function)
 *         int nind            estimated number of independent points
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 * v10Sep98 CDF, Make the correlation coefficient a Fluxrec rather than float
 *               Put Bevington probability (from bevington_prob) into
 *                the err member of r.
 */

int corr_bevington(Fluxrec *flux1, Fluxrec *flux2, int npoints, 
		   Fluxrec *r, int nind)
{
  int i;              /* Looping variable */
  float mean1,mean2;  /* Statistically weighted means of each curve */
  float diff1,diff2;  /* Difference between a point and the curve mean */
  float sig12,sig22;  /* Square of error for each point */
  float wsum1,wsum2;  /* Sum of weights of the two curves */
  float wcovar;       /* Weighting for each point in covariance calculation */
  float wsum12;       /* Sum of wcovar */
  float s1,s2;        /* Variances of each curve */
  float covar;        /* Covariance of the two curves */
  Fluxrec *fp1,*fp2;  /* Pointers to navigate flux1 and flux2 */

  /*
   * Calculate the weighted means of the two curves
   */

  if(calc_swmean(flux1,npoints,&mean1)) {
    fprintf(stderr,"ERROR: corr_bevington\n");
    return 1;
  }

  if(calc_swmean(flux2,npoints,&mean2)) {
    fprintf(stderr,"ERROR: corr_bevington\n");
    return 1;
  }

  /*
   * Initialize all the running sums
   */

  wsum1 = wsum2 = wsum12 = 0.0;
  s1 = s2 = covar = 0.0;

  /*
   * Now start the loop, calculating running sums as we go.
   */

  for(i=0,fp1=flux1,fp2=flux2; i<npoints; i++,fp1++,fp2++) {

    diff1 = fp1->flux - mean1;
    diff2 = fp2->flux - mean2;
    sig12 = fp1->err * fp1->err;
    sig22 = fp2->err * fp2->err;

    /*
     * Calculations for s1
     */

    s1 += diff1 * diff1 / sig12;
    wsum1 += 1.0 / sig12;

    /*
     * Calculations for s2
     */

    s2 += diff2 * diff2 / sig22;
    wsum2 += 1.0 / sig22;

    /*
     * Calculations for covar
     */

    wcovar = 1.0 / (sig12 + sig22);
    covar += diff1 * diff2 * wcovar;
    wsum12 += wcovar;

  }

  /*
   * Now normalize s1, s2, and covar properly
   */
  s1 *= npoints / ((npoints-1) * wsum1);
  s2 *= npoints / ((npoints-1) * wsum2);
  covar *= npoints / ((npoints-1) * wsum12);

  /*
   * Calculate the linear correlation coefficient and probability
   *  of getting that coefficient.  Put the coefficient in r->flux
   *  and the probability into r->err.
   */

  r->flux = covar / (sqrt(s1) * sqrt(s2));
  r->err = bevington_prob(*r,nind);

  return 0;
}

/*.......................................................................
 *
 * Function bevington_prob
 *
 * Calculates the probabilities of getting a certain value of linear
 *  correlation coefficient (r) from two uncorrelated curves.
 *
 *                  1       Gamma[(nu+1)/2]
 *  P_r (r,nu) = -------- * --------------- * (1 - r^2)^((nu-2)/2)
 *               sqrt(PI)     Gamma[nu/2]
 *
 *  where nu = N-2 is the number of degrees of freedom.  N is the number
 *  of independent points in the curves, which may differ from the number
 *  of points in the curve (N_tot).
 * This function gets Gamma(z) by calling the Numerical Recipes function 
 *  gammln. 
 * The probability of getting r or higher from two uncorrelated curves is:
 *
 *                (1
 *  P_c (r,N) = 2 |    P_r (rho,nu) d_rho
 *                )|r|
 *
 * NB: This only works if r has been calculated in corr_bevington.
 *
 * Inputs: Fluxrec maxcorr     container for r and N_tot (in flux and match
 *                              members)
 *         int nind            number of independent points.
 *
 * Output: double prob         probability of getting r or higher
 *
 * v12Sep98 CDF, Modification to let N differ from N_tot if necessary.
 */

double bevington_prob(Fluxrec maxcorr, int nind)
{
  int doprint=0;   /* Set to 1 for printout */
  float nu;        /* Number of degrees of freedom (maxcorr.match - 2) */
  double gammfac;  /* Factor from the ratio of Gamma functions */
  double rho;      /* Temporary value of r for numerical integration */
  double drho;     /* Numerical integration step width */
  double sum=0.0;  /* Running sum in numerical integration */
  double power;    /* Exponent in integrand */
  double prob;     /* Probability of getting r or greater */

  /*
   * Calculate the gamma values and include the 1/sqrt(PI) factor
   */

  nu = 1.0 * nind - 2.0;
  gammfac = exp(gammln((nu+1.0)/2.0) - gammln(nu/2.0))/sqrt(PI);
  if(doprint) {
    printf("bevington_prob: N_ind = %d, N_tot = %d\n",nind,maxcorr.match);
    printf("bevington_prob: nu = %3.0f, gammfac = %f\n",nu,gammfac);
  }

  /*
   * Now numerically integrate the (1 - rho^2)^((nu-2)/2) term
   */

  rho = maxcorr.flux;
  if(rho < 0)
    rho *= -1.0;
  if(doprint)
    printf("bevington_prob: rho = %f\n",rho);
  drho = (1.0 - rho)/100.0;
  rho += drho/2.0;
  while(rho <= 1.0) {
    power = (nu - 2.0)/2.0;
    sum += drho * pow((1 - rho*rho),power);
    rho += drho;
  }

  /*
   * Calculate the probability
   */

  prob = 2.0 * gammfac * sum;
  if(doprint)
    printf("bevington_prob: prob = %e\n\n",prob);

  return prob;
}

/*.......................................................................
 *
 * Function corr_white_peterson
 *
 * Calculates a correlation coefficient according to the normalization of
 *  White and Peterson (see Book 6, p. 53).
 *
 * Inputs: Fluxrec *flux1      first flux (x_j)
 *         Fluxrec *flux2      second flux (x_k)
 *         int npoints         number of points in each curve (N)
 *         float *r            linear correlation coefficient (set by this
 *                              function)
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int corr_white_peterson(Fluxrec *flux1, Fluxrec *flux2, int npoints, float *r)
{
  int i;              /* Looping variable */
  float mean1,mean2;  /* Straight means of each curve */
  float rms1,rms2;    /* RMS scatter about means of each curve */
  Fluxrec *fp1,*fp2;  /* Pointers to navigate flux1 and flux2 */

  /*
   * First calculate mean and rms of the two curves
   */

  if(calc_mean(flux1,npoints,&mean1,&rms1)) {
    fprintf(stderr,"ERROR: corr_white_peterson\n");
    return 1;
  }
  if(calc_mean(flux2,npoints,&mean2,&rms2)) {
    fprintf(stderr,"ERROR: corr_white_peterson\n");
    return 1;
  }

  /*
   * Initialize the linear correlation coefficient
   */

  *r = 0.0;

  /*
   * Now step through the curves, calculating the
   *  sum of the products of the fluxes decreased by the local mean.
   */

  for(i=0,fp1=flux1,fp2=flux2; i<npoints; i++,fp1++,fp2++)
    *r += (fp1->flux - mean1) * (fp2->flux - mean2);

  /*
   * Now normalize according to White and Peterson.
   */

  *r /= (npoints * rms1 * rms2);

  return 0;
}

/*.......................................................................
 *
 * Function zero_pad
 *
 * Zero-pads a light curve in preparation for cross-correlation.
 *  If the light curve has npoint data points, then the zero-padded
 *  curve will have 2*npoint points, with the first and last npoints/2 
 *  data points being set to zero.
 *
 * NB: This function assumes that the curve is regularly sampled.
 *
 * Inputs: Fluxrec *flux       unpadded light curve
 *         int npoints         number of points in the light curves
 *
 * Output: float *zpad         padded curve
 *
 */

float *zero_pad(Fluxrec *flux, int npoints)
{
  int i;             /* Looping variable */
  float *zpad=NULL;  /* Padded light curve */
  float *zptr;       /* Pointer to navigate zpad */
  Fluxrec *fptr;     /* Pointer to navigate flux */

  /*
   * Allocate memory for the padded array
   */

  if(!(zpad = new_array(npoints,ZEROFAC))) {
    fprintf(stderr,"ERROR: zero_pad.\n");
    return NULL;
  }

#if 0
  /*
   * Put the data from flux in the middle of the padded
   *  array container.  Note that new_array initializes the
   *  members of the array to zero, so there is no need to set the
   *  first and last (npoints/2) points to zero in this function.
   */

  for(i=0,fptr=flux,zptr=zpad+(npoints/2); i<npoints; i++,fptr++,zptr++) 
    *zptr = fptr->flux;
#endif
  /*
   * Put the data into the first npoints members in the array and put
   *  zeros in the last npoints members.  Note that since new_array
   *  sets all members of the new array to 0.0, all we have to do
   *  here is put the members of flux into the first npoints members
   *  of zpad.
   */

  for(i=0,fptr=flux,zptr=zpad; i<npoints; i++,fptr++,zptr++)
    *zptr = fptr->flux;

#if 0
  /*
   * Convert the padded array into Numerical Recipes' "wrap-around"
   *  order -- that is with points (npoints/2) to (npoints - 1) in
   *  positions 0 to (npoints/2 - 1) and points 0 to (npoints/2 - 1)
   *  in positionis (3*npoints/2) to (2*npoints - 1).  The middle
   *  npoints/2 points will be 0.0, since new_array sets all members
   *  of the new array to 0.0.  Thus we will have a zero-padded array.
   */

  for(i=0; i<npoints/2; i++)
    zpad[i + 3*npoints/2] = flux[i].flux;

  for(i=npoints/2; i<npoints; i++)
    zpad[i-npoints/2] = flux[i].flux;
#endif

  return zpad;
}

/*.......................................................................
 *
 * Function wt_hanning
 *
 * Takes a Fluxrec array and applies Hanning weighting to the flux values,
 *  returning a float array of the weighted fluxes.  The Hanning weighting
 *  is defined as:
 *
 *  w_i = 1/2 * [ 1 - cos(2*PI*i/npoints) ]
 *
 *  where i runs from 0 to npoints-1.
 *
 * Inputs: Fluxrec *flux       raw light curve
 *         int npoints         number of points in the light curves
 *
 * Output: float *zpad         weighted curve
 *
 */

float *wt_hanning(Fluxrec *flux, int npoints)
{
  int i;               /* Looping variable */
  float weight;        /* Weighting for any one point */
  float *hann=NULL;    /* Weighted curve */
  float *hptr;         /* Pointer to navigate hann */
  Fluxrec *fptr;       /* Pointer to navigate flux */

  /*
   * Allocate memory for weighted curve
   */

  if(!(hann = new_array(npoints,1))) {
    fprintf(stderr,"ERROR: wt_hann\n");
    return NULL;
  }

  /*
   * Loop through the Fluxrec array, applying the weighting
   */

  for(i=0,hptr=hann,fptr=flux; i<npoints; i++,hptr++,fptr++) {
    weight = 0.5 * (1 - cos(2*PI*i/npoints));
    *hptr = fptr->flux * weight;
  }

  return hann;

}

/*.......................................................................
 *
 * Function call_acorr
 *
 * Does a series of auto-correlations of light curves, by
 *  calling the function do_corr for each curve.  Puts the
 *  output into a file designated by filename.
 *
 * Inputs: Fluxrec *flux[]     light curves
 *         int size            number of points in each curve
 *         char *filename      name of output file
 *         int nbad            number of bad data points
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> failure
 *
 */

int call_acorr(Fluxrec *flux[], int size, char *filename, int nbad)
{
  int i;                       /* Looping variable */
  int no_error=1;              /* Flag set to 0 on error */
  int corsize;                 /* Length of correlation curves */
  Fluxrec *zmean[N08]={NULL};  /* Zero-mean light curves */
  Fluxrec *corra=NULL;         /* Auto-correlation for A */
  Fluxrec *corrb=NULL;         /* Auto-correlation for B */
  Fluxrec *corrc=NULL;         /* Auto-correlation for C */
  Fluxrec *corrd=NULL;         /* Auto-correlation for D */
  Fluxrec *aptr,*bptr;         /* Pointers to navigate auto-corr arrays */
  Fluxrec *cptr,*dptr;         /* Pointers to navigate auto-corr arrays */
  FILE *ofp=NULL;              /* Output file */

  /*
   * Normalize the curves since we are looking for correlations in
   *  fractional variations and then subtract 1 to create a zero-mean
   *  light curve.  This helps by getting rid of the DC component in
   *  the autocorrelation.
   */

  for(i=0; i<N08; i++) {
    if(no_error)
      if(!(zmean[i] = norm_zero_mean(flux[i],size)))
	no_error = 0;
  }

  /*
   * Do the autocorrelations
   */

  if(!(corra = do_corr(zmean[0],zmean[0],size,&corsize,nbad))) {
    fprintf(stderr,"ERROR: call_acorr\n");
    return 1;
  }

  if(!(corrb = do_corr(zmean[1],zmean[1],size,&corsize,nbad)))
     no_error = 0;

  if(no_error)
    if(!(corrc = do_corr(zmean[2],zmean[2],size,&corsize,nbad)))
       no_error = 0;

  if(no_error)
    if(!(corrd = do_corr(zmean[3],zmean[3],size,&corsize,nbad)))
      no_error = 0;

  if(no_error)
    if((ofp = fopen(filename,"w")) == NULL) {
      fprintf(stderr,"ERROR: call_acorr.  Cannot open %s\n",filename);
      no_error = 0;
    }

  if(no_error) {
    fprintf(ofp,"#  Lag   Acorr_a  Acorr_b  Acorr_c  Acorr_d\n");
    fprintf(ofp,"#------- -------- -------  -------  -------\n");
    for(i=0,aptr=corra,bptr=corrb,cptr=corrc,dptr=corrd; i<corsize; 
	i++,aptr++,bptr++,cptr++,dptr++)
      fprintf(ofp,"%8.3f %8.5f %8.5f %8.5f %8.5f\n",
	      aptr->day,aptr->flux,bptr->flux,cptr->flux,dptr->flux);
  }

  /*
   * Clean up
   */

  for(i=0; i<N08; i++)
    zmean[i] = del_fluxrec(zmean[i]);
  corra = del_fluxrec(corra);
  corrb = del_fluxrec(corrb);
  corrc = del_fluxrec(corrc);
  corrd = del_fluxrec(corrd);
  if(ofp)
    fclose(ofp);

  if(no_error)
    return 0;
  else {
    fprintf(stderr,"ERROR: call_acorr\n");
    return 1;
  }
}

