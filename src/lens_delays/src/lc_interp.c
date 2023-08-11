/*
 * lc_interp.c
 *
 * A library of functions to perform interpolation and smoothing of
 *  the 1608 light curves.
 *
 * 25Sep98 CDF,  Split off from lc_funcs.c.  For older history, see header
 *                in lc_funcs.c.
 * v12Sep00 CDF, Deleted the grid_size function and incorporated functionality
 *                into interp_switch.
 *               Deleted the interp_curve and smooth_1608 functions since
 *                they just existed to print call the actual interpolation
 *                and smoothing functions and then to print out the 
 *                interpolated curves.  The printing of the curves is now
 *                handled externally to those functions so they can be
 *                scrapped and the interpolation/smoothing functions 
 *                (lin_interp, csmooth, and smooth_in_place) can be called 
 *                directly.
 *               Got rid of the SMTHENINT and SMIPINT options since these
 *                can now be emulated by running interp.c twice -- the
 *                first time w/SMONLY or SMINPLACE and the second time
 *                with INTONLY.
 * v13Sep00 CDF, Got rid of *gridsize variable passed to interp_switch since
 *                it is now incorporated in the setup structure as ninterp.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "nr.h"
#include "nrutil.h"
#include "structdef.h"
#include "lc_funcs.h"
#include "lc_interp.h"

/*.......................................................................
 *
 * Function interp_switch
 *
 * Uses flag in setup structure to call the appropriate 
 *  smoothing/interpolation function
 *
 * Inputs: Fluxrec *raw        raw light curve
 *         int nraw            number of points in raw light curve
 *         Setup *setup        container with switch flag
 *         int nbad            number of bad days flagged
 *
 * Output: Fluxrec *interp     smoothed/interpolated light curve
 *
 * v12Sep00 CDF, Got rid of the chisq log and cross-correlation log file
 *                pointers, which weren't being used by any program.  Also
 *                got rid of the doprint flag since the printing is now
 *                done externally.
 *               Changed calls to interp_curve to direct calls to
 *                lin_interp and changed calls to smooth_1608 to calls
 *                to csmooth and smooth_in_place.
 * v13Sep00 CDF, Got rid of passed gridsize variable since it is now
 *                incorporated into the setup structure as the ninterp
 *                parameter.
 */

Fluxrec *interp_switch(Fluxrec *raw, int nraw, Setup *setup, int nbad)
{
  int no_error=1;        /* Flag set to 0 on error */
  Fluxrec *interp=NULL;  /* Smoothed/interpolated curve to be returned */
  Fluxrec *temp=NULL;    /* Temporary container for two-step interpolation */

  /*
   * Call the appropriate interpolation and/or smoothing function.
   */

  switch(setup->dosmooth) {
  case INTONLY:

    /*
     * Interpolate
     */

    if(!(interp = lin_interp(raw,nraw,setup->ninterp,setup->intstep,nbad)))
      no_error = 0;
    break;

  case SMONLY:

    /*
     * Smooth
     */

    if(!(interp = csmooth(raw,nraw,setup)))
      no_error = 0;
    break;

  case SMINPLACE:

    /*
     * Smooth
     */

    if(!(interp = smooth_in_place(raw,nraw,setup)))
      no_error = 0;
    break;


  default:
    fprintf(stderr,"ERROR: interp_switch: ");
    fprintf(stderr,"Invalid smoothing/interpolation choice.\n");
    no_error = 0;
  }

  /*
   * Clean up and exit
   */

  temp = del_fluxrec(temp);

  if(no_error)
    return interp;
  else {
    fprintf(stderr,"ERROR: interp_switch\n");
    return del_fluxrec(interp);
  }
}

/*.......................................................................
 *
 * Function lin_interp
 *
 * Given a function x and y does a linear interpolation onto a regularly
 *  sampled grid in x.
 *
 * Inputs: Fluxrec *data       array of x and y values
 *         int ndata           size of irregularly sampled array
 *         int ngrid           size of grid
 *         float dx            step size for regulary sampled grid
 *         int nbad            number of bad days flagged
 *
 * Output: Fluxrec *fluxrec    flux record with regular sampling
 *                             NULL on error.
 *
 */

Fluxrec *lin_interp(Fluxrec *data, int ndata, int ngrid, float dx, int nbad)
{
  int i;              /* Looping variable */
  int no_error=1;     /* Flag set to 0 on error */
  int ngood=0;        /* Number of good data points. Should = ndata - nbad */
  Fluxrec *tmp=NULL;  /* Temporary light curve containing unflagged data */
  Fluxrec *m=NULL;    /* Slopes between measured points, errors */
  Fluxrec *mptr;      /* Pointer used to navigate m array */
  Fluxrec *b=NULL;    /* y-intercepts of lines connecting measured points */
  Fluxrec *bptr;      /* Pointer used to navigate b array */
  Fluxrec *dataptr;   /* Pointer for navigating the data array */
  Fluxrec *grid=NULL; /* Regularly sampled grid */
  Fluxrec *flptr;     /* Pointer used to navigate fluxrec arrays */

  /*
   * Create a regularly sampled grid in x
   */

  if(!(grid = new_fluxrec(ngrid))) {
    return NULL;
  }

  /*
   * Allocate tmp array and transfer unflagged points from data array
   *  into it.
   */

  if(!(tmp = new_fluxrec(ndata-nbad)))
    no_error = 0;
  else {
    flptr = tmp;
    for(i=0,dataptr=data; i<ndata; i++,dataptr++) {
      if(dataptr->match >= 0) {
	ngood++;
	if(ngood > ndata-nbad) {
	  fprintf(stderr,"ERROR: lin_interp. Problem resampling curves.\n");
	  no_error = 0;
	}
	else {
	  *flptr = *dataptr;
	  flptr++;
	}
      }
    }
    if(ngood > ndata-nbad) {
      fprintf(stderr,"ERROR: lin_interp. ngood (%d) != ndata-nbad (%d).\n",
	      ngood,ndata-nbad);
      no_error = 0;
    }
    else
      printf(" lin_interp: There are %d good data points in the light curve.\n",
	     ngood);
  }

  /*
   * Allocate memory for slope and y-intercept arrays
   */

  if(no_error)
    if(!(m = new_fluxrec(ngood-1)))
      no_error = 0;

  if(no_error)
    if(!(b = new_fluxrec(ngood-1)))
      no_error = 0;

  /*
   * Fill in x part of grid array
   */

  if(no_error)
    for(i=0,flptr=grid; i<ngrid; i++,flptr++)
      flptr->day = data->day + i*dx;

  /*
   * Calculate slopes and y-intercepts, using straightforward method:
   *
   *        y2 - y1
   *    m = -------     and  b = y1 - m x1 
   *        x2 - x1 
   *
   */

  if(no_error) {
    for(i=0,dataptr=tmp,mptr=m,bptr=b; i<ngood-1; 
	i++,dataptr++,mptr++,bptr++) {

      /*
       * Error check to avoid dividing by zero
       */

      if(((dataptr+1)->day - dataptr->day) == 0) {
	fprintf(stderr,"ERROR: lin_interp. Two days with the same value\n");
	no_error = 0;
      }

      /*
       * Calculate slope
       */

      mptr->flux = ((dataptr+1)->flux - dataptr->flux) / 
	((dataptr+1)->day - dataptr->day);
      mptr->err = ((dataptr+1)->err - dataptr->err) / 
	((dataptr+1)->day - dataptr->day);

      /*
       * Calculate intercept
       */

      bptr->flux = dataptr->flux - mptr->flux * dataptr->day;
      bptr->err = dataptr->err - mptr->err * dataptr->day;

    }
  }

  /*
   * Calculate interpolated values
   */

  if(no_error) {
    dataptr = tmp;
    mptr = m;
    bptr = b;
    for(i=0,flptr=grid; i<ngrid; i++,flptr++) {

      /*
       * Increment the observed data pointer until the grid point
       *  is between two data points.
       */

      while((flptr->day > (dataptr+1)->day) && 
	    (dataptr < (tmp+ngood-2)) && no_error) {
	dataptr++;
	mptr++;
	bptr++;
      }

      /*
       * Make sure that we're not trying to calculate a grid point
       *  outside of our observed data.
       */

      if((dataptr == tmp+ngood-2) && (flptr->day > (dataptr+1)->day)) {
	fprintf(stderr,"ERROR: lin_interp. Grid outside data range\n");
	no_error = 0;
      }

      /*
       * Finally do the interpolation
       */

      else {
	flptr->flux = mptr->flux * flptr->day + bptr->flux;
	flptr->err = mptr->err * flptr->day + bptr->err;

	/*
	 * Mark point as an observed point if it is within +/-
	 *  DAYSTEP/2 of an observed point.
	 */

	if(fabs(flptr->day - dataptr->day) < (0.5*DAYSTEP) ||
	   fabs(flptr->day - (dataptr+1)->day) < (0.5*DAYSTEP)) {
	  flptr->match = 1;
	}
	else {
	  flptr->match = 0;
	}
      }
    }
  }

  /*
   * Clean up
   */

  m = del_fluxrec(m);
  b = del_fluxrec(b);
  tmp = del_fluxrec(tmp);

  if(no_error)
    return grid;
  else {
    fprintf(stderr,"ERROR: lin_interp\n");
    return NULL;
  }
}

/*.......................................................................
 *
 * Function csmooth
 *
 * Creates a smoothed and interpolated light curve with the specified smoothing 
 *  function.
 *
 * Inputs: Fluxrec *raw        unsmoothed light curve
 *         int nraw            number of points in input curve
 *         Setup *setup        container for width of smoothing window,
 *                              step size for gridded output curve and
 *                              type of smoothing function
 *
 * Output: Fluxrec *smooth     smoothed and gridded output curve.
 *
 * v30Aug98 CDF, Moved memory allocation for the smoothed curve and the
 *                initialization of the dates on the regularly stepped
 *                grid into this function from the various smoothing 
 *                functions.
 * v12Sep00 CDF, Moved determination of number of points in smoothed
 *                curve up into the calling function.
 * v13Sep00 CDF, Got rid of the passed nsmooth parameter since it is now
 *                contained in the setup structure as ninterp.
 */

Fluxrec *csmooth(Fluxrec *raw, int nraw, Setup *setup)
{
  int i;                 /* Looping variable */
  int no_error=1;        /* Flag set to 0 on error */
  Fluxrec *smooth=NULL;  /* Smoothed curve */
  Fluxrec *sptr;         /* Pointer to navigate smooth */

  /*
   * Allocate memory for smoothed curve
   */

  if(!(smooth = new_fluxrec(setup->ninterp))) {
    fprintf(stderr,"ERROR: csmooth\n");
    return NULL;
  }

  /*
   * Fill smoothed curve day pointer with the interpolated values
   */

  for(i=0,sptr=smooth; i<setup->ninterp; i++,sptr++)
      sptr->day = setup->intstart + i * setup->intstep;

  /*
   * Smooth the raw curve with the chosen smoothing function
   */

  switch(setup->smtype) {
  case SMBOXCAR:
    if(boxcar(raw,nraw,smooth,setup->ninterp,setup->smwidth,
	      setup->intstep))
      no_error = 0;
    break;
  case SMMEDIAN:
    if(medsmooth(raw,nraw,smooth,setup->ninterp,setup->smwidth,
		 setup->intstep))
      no_error = 0;
    break;
  case SMTRIANGLE:
    if(triangle(raw,nraw,smooth,setup->ninterp,setup->smwidth,
		setup->intstep))
      no_error = 0;
    break;
  case SMVARBOX:
    if(varbox(raw,nraw,smooth,setup->ninterp,setup->nvar,setup->intstep))
      no_error = 0;
    break;
  case SMVARTRI:
    if(vartri(raw,nraw,smooth,setup->ninterp,setup->nvar,setup->intstep))
      no_error = 0;
    break;
  case SMGAUSS:
    if(gaussian(raw,nraw,smooth,setup->ninterp,setup->smwidth,setup->intstep))
      no_error = 0;
    break;
  default:
    fprintf(stderr,"ERROR: csmooth. Invalid smoothing function.\n");
    no_error = 0;
  }

  /*
   * Exit
   */

  if(no_error)
    return smooth;
  else {
    fprintf(stderr,"ERROR: csmooth\n");
    return NULL;
  }
}

/*.......................................................................
 *
 * Function smooth_in_place
 *
 * Creates a smoothed light curve by replacing each point in the raw
 *  light curve with the value obtained by centering the smoothing
 *  window at that point and performing the specified smoothing 
 *  function.
 * NB: The output curve will have the same sampling pattern as the
 *  input curve, i.e. it will NOT be interpolated.
 *
 * Inputs: Fluxrec *raw        unsmoothed light curve
 *         int nraw            number of points in input curve
 *         Setup *setup        container for width of smoothing window,
 *                              step size for gridded output curve and
 *                              type of smoothing function
 *
 * Output: Fluxrec *smooth     smoothed and but ungridded output curve.
 *
 * 30Aug98 CDF, A variation of csmooth.
 */

Fluxrec *smooth_in_place(Fluxrec *raw, int nraw, Setup *setup)
{
  int i;                 /* Looping variable */
  int no_error=1;        /* Flag set to 0 on error */
  Fluxrec *smooth=NULL;  /* Smoothed curve */
  Fluxrec *rptr;         /* Pointer to navigate raw */
  Fluxrec *sptr;         /* Pointer to navigate smooth */

  /*
   * Allocate memory for smoothed curve
   */

  if(!(smooth = new_fluxrec(nraw))) {
    fprintf(stderr,"ERROR: smooth_in_place\n");
    return NULL;
  }

  /*
   * Fill smoothed curve day pointer with the values from the raw curve
   */

  for(i=0,sptr=smooth,rptr=raw; i<nraw; i++,sptr++,rptr++) {
      sptr->day = rptr->day;
#if 0
      sptr->match = rptr->match;
#endif
  }

  /*
   * Smooth the raw curve with the chosen smoothing function
   */

  switch(setup->smtype) {
  case SMBOXCAR:
    if(boxcar(raw,nraw,smooth,nraw,setup->smwidth,setup->intstep))
      no_error = 0;
    break;
  case SMMEDIAN:
    if(medsmooth(raw,nraw,smooth,nraw,setup->smwidth,setup->intstep))
      no_error = 0;
    break;
  case SMTRIANGLE:
    if(triangle(raw,nraw,smooth,nraw,setup->smwidth,setup->intstep))
      no_error = 0;
    break;
  case SMVARBOX:
    if(varbox(raw,nraw,smooth,nraw,setup->nvar,setup->intstep))
      no_error = 0;
    break;
  case SMVARTRI:
    if(vartri(raw,nraw,smooth,nraw,setup->nvar,setup->intstep))
      no_error = 0;
    break;
  case SMGAUSS:
    if(gaussian(raw,nraw,smooth,nraw,setup->smwidth,setup->intstep))
      no_error = 0;
    break;
  default:
    fprintf(stderr,"ERROR: smooth_in_place. Invalid smoothing function.\n");
    no_error = 0;
  }

  /*
   * Exit
   */

  if(no_error)
    return smooth;
  else {
    fprintf(stderr,"ERROR: smooth_in_place\n");
    return NULL;
  }
}

/*.......................................................................
 *
 * Function boxcar
 *
 * Performs boxcar smoothing on a light curve, producing an output
 *  curve sampled on a regular grid
 *
 * Inputs: Fluxrec *raw        unsmoothed and ungridded curve
 *         int nraw            number of points in unsmoothed curve
 *         Fluxrec *smooth     container for smoothed values.  The "day"
 *                              members of the smoothed curves are already set.
 *         int nsmooth         number of points in smoothed curve.
 *         float width         width of boxcar smoothing
 *         float step          step size of gridded output
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 * v13Jul98 CDF, Instead of cycling through raw again to calculate rms
 *                scatter about mean, just allocate a temporary array
 *                to hold the values used to find the mean value.
 * v14Jul98 CDF, Changed from straight mean to variance-weighted mean
 * v30Aug98 CDF, Now only set match to 1 if the smooth curve point is
 *                within 0.25 days of a raw curve point.
 *               Moved allocation of smoothed curve array up to csmooth.
 */

int boxcar(Fluxrec *raw, int nraw, Fluxrec *smooth, int nsmooth, float width, 
	   float step)
{
  int i,j;               /* Looping variables */
  int no_error=1;        /* Flag set to 0 on error */
  int meancount;         /* Number of points used for calculating mean */
  int nmatch=0;          /* Number of smooth points that match raw points */
  float rms;             /* RMS scatter about mean value */
  float *pweight=NULL;   /* Container for smoothing weights (all = 1.0) */
  float *pptr;           /* Pointer to navigate pweight */
  float sig2;            /* Variance of each point */
  float wtsum;           /* Sum of weights from calculating weighted mean */
  Fluxrec *temp=NULL;    /* Container for points used to calculate mean */
  Fluxrec *rptr;         /* Pointer for navigating raw lightcurve */
  Fluxrec *sptr;         /* Pointer for navigating smoothed lightcurve */
  Fluxrec *tptr;         /* Pointer for navigating temp */

  /*
   * Allocate container for points to calculate mean.  There will never be
   *  more points in this container than in the raw light curve, so
   *  set it to nraw and then just use the first meancount values.
   */

  if(!(temp = new_fluxrec(nraw)))
    no_error = 0;

  /*
   * Allocate memory for smoothing function weights (which in this case
   *  will all be set to 1.0.  Also set the size of this array to the
   *  maximum possible value, which is nraw.
   */

  if(!(pweight = new_array(nraw,1)))
    no_error = 0;

  /*
   * Step through the smoothed curve, calculating mean and error in the
   *  mean at each step.
   * Remember that an unbiassed estimator of the variance of a sample 
   *  is given by:
   *
   *        sum((x_i - mean)^2)
   *  s^2 = -------------------
   *             n - 1
   *
   *  and thus the error on the estimate of the mean value is:
   *
   *  sigma_mean = sqrt (s^2 / n).
   *
   * This value for sigma_mean will be compared to the error in the mean
   *  calculated using var_wmean and the larger will be used.
   */

  if(no_error) {
    for(i=0,sptr=smooth; i<nsmooth; i++,sptr++) {
      
      /*
       * Initialize pointers and running sums
       */

      meancount = 0;
      tptr = temp;
      pptr = pweight;
      wtsum = 0.0;

      /*
       * Get the points used in calculating the mean.
       */

      for(j=0,rptr=raw; j<nraw; j++,rptr++) {
	if(fabs(rptr->day - sptr->day) < width/2.0 && rptr->match >= 0) {
	  meancount++;
	  sig2 = rptr->err * rptr->err;
	  sptr->flux += rptr->flux / sig2;
	  sptr->err += sig2;
	  wtsum += 1/sig2;
	  *tptr = *rptr;
	  tptr++;
	  *pptr = 1.0;
	  pptr++;

	  /*
	   * Set match = 1 if sptr is within one fourth of a day of a
	   *  day on the raw curve
	   */

	  if(fabs(sptr->day - rptr->day) < 0.25) {
	    sptr->match = 1;
	    nmatch++;
	  }
	}
      }

      /*
       * Now calculate the mean and the uncertainty in the mean at each
       *  point.  There are three possible cases:
       *    1. No points in window -- Set mean and err to mean and err of
       *                              previous point or, if this is the
       *                              first point, set to flux density and 
       *                              err of first point in raw curve.
       *    2. One point in window -- Set mean and err to flux density and 
       *                              err of that point.
       *    3. Two or more points  -- Call var_wmean to calculate the variance
       *                              in the weighted mean and then calculate
       *                              the RMS scatter about the mean.  Take
       *                              the larger of these two values.
       */

      switch(meancount) {
      case 0:
	if(sptr > smooth) {
	  sptr->flux = (sptr-1)->flux;
	  sptr->err = (sptr-1)->err;
	}
	else {
	  sptr->flux = raw->flux;
	  sptr->err = raw->err;
	}
	break;
      case 1:
	sptr->flux /= wtsum;
	sptr->err = sqrt(sptr->err);
	break;
      default:
	sptr->flux /= wtsum;

	/*
	 * First calculate sigma of the weighted mean
	 */

	if(var_wmean(temp,pweight,meancount,&sptr->err)) {
	  no_error = 0;
	  break;
	}

	/*
	 * Now calculate rms scatter about mean value.  Note that sptr->flux
	 *  has been set to the proper variance-weighted mean above.
	 */

	rms = 0.0;
	for(j=0,tptr=temp; j<meancount; j++,tptr++)
	  rms += (tptr->flux - sptr->flux) * (tptr->flux - sptr->flux);
	rms = sqrt(rms/(meancount * (meancount - 1)));

	/*
	 * Compare the two values and take the larger
	 */
#if 0
	if(rms > sptr->err)
	  sptr->err = rms;
#endif
      }
    }
  }

  /*
   * Clean up and exit
   */

  printf("  boxcar: There are %d smooth points that match raw points.\n",
	 nmatch);
  temp = del_fluxrec(temp);
  pweight = del_array(pweight);

  if(no_error)
    return 0;
  else {
    fprintf(stderr,"ERROR: boxcar\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function medsmooth
 *
 * Smooths a light curve by finding the median value in a window that is
 *  shifted in regular steps, producing an output curve sampled on a regular
 *  grid.  The median is found by sorting the points in the window and
 *  taking either the middle point (odd number of points) or the average
 *  of the two middle points (even number of points).
 *
 * Inputs: Fluxrec *raw        unsmoothed and ungridded curve
 *         int nraw            number of points in unsmoothed curve
 *         Fluxrec *smooth     container for smoothed values.  The "day"
 *                              members of the smoothed curves are already set.
 *         int nsmooth         number of points in smoothed curve.
 *         float width         width of smoothing window
 *         float step          step size of gridded output
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 * v30Aug98 CDF, Now only set match to 1 if the smooth curve point is
 *                within 0.25 days of a raw curve point.
 *               Moved allocation of smoothed curve array up to csmooth.
 */

int medsmooth(Fluxrec *raw, int nraw, Fluxrec *smooth, int nsmooth, float width, 
	      float step)
{
  int i,j;               /* Looping variables */
  int meancount;         /* Number of points used for calculating mean */
  int nmatch=0;          /* Number of smooth points that match raw points */
  float rms;             /* RMS scatter about mean value */
  float *medarray=NULL;  /* Array containing values used in calculating median */
  float *medptr;         /* Pointer used to navigate medarray */
  Fluxrec *rptr;         /* Pointer for navigating raw lightcurve */
  Fluxrec *sptr;         /* Pointer for navigating smoothed lightcurve */

  /*
   * Allocate memory for medianing array.  The array will never (by definition)
   *  contain more points than the input curve, so just set the size to nraw 
   *  and when calling the medianing function (hpsort), just use the first n
   *  values in the array.
   */

  if(!(medarray = new_array(nraw,1))) {
    fprintf(stderr,"ERROR: medsmooth\n");
    return 1;
  }

  /*
   * Step through the smoothed curve, calculating median and uncertainty
   *  in the estimate of the median at each step.
   * Remember that an unbiassed estimator of the variance of a sample 
   *  is given by:
   *
   *        sum((x_i - mean)^2)
   *  s^2 = -------------------
   *             n - 1
   *
   *  and thus the error on the estimate of the mean value is:
   *
   *  sigma_mean = s^2 / sqrt(n).
   *
   * We can use this error estimate for now, even though we're taking
   *  a median rather than a mean.
   *
   */

  for(i=0,sptr=smooth; i<nsmooth; i++,sptr++) {
    meancount = 0;
    medptr = medarray;

    /*
     * Get the points used to calculate the median
     */

    for(j=0,rptr=raw; j<nraw; j++,rptr++) {
      if(fabs(rptr->day - sptr->day) < width/2.0 && rptr->match >= 0) {
	meancount++;
	*medptr = rptr->flux;
	sptr->err = rptr->err;
	medptr++;

	/*
	 * Set match = 1 if sptr is within one fourth of a day of a
	 *  day on the raw curve
	 */

	if(fabs(sptr->day - rptr->day) < 0.25) {
	  sptr->match = 1;
	  nmatch++;
	}
      }
    }

    /*
     * Now calculate the median and the uncertainty in the median at each
     *  point.  There are three possible cases:
     *    1. No points in window -- Set median and err to median and err of
     *                              previous point or, if this is the
     *                              first point, set to flux density and err 
     *                              of first point in raw curve.
     *    2. One point in window -- Set median and err to flux density and 
     *                              err of that point.
     *    3. Two or more points  -- Call var_wmean to calculate the variance
     *                               in the weighted mean and then calculate
     *                               the RMS scatter about the mean.  Take
     *                               the larger of these two values.
     */

    switch(meancount) {
    case 0:
      if(sptr > smooth) {
	sptr->flux = (sptr-1)->flux;
	sptr->err = (sptr-1)->err;
      }
      else {
	sptr->flux = raw->flux;
	sptr->err = raw->err;
      }
      break;
    case 1:
      sptr->flux = *medarray;
      break;
    default:
      
      /* 
       * Perform medianing by calling the Numerical Recipes function
       *  hpsort which sorts the array.  Claim that medarray contains only
       *  meancount values.  After sorting, take either the middle points
       *  if meancount is odd, or take the average of the two middle points
       *  if meancount is even.
       *
       * NB: The array pointer must be set to medarray-1 to use NR's
       *  Fortran-order functions.
       */

      hpsort(meancount,medarray-1);
      if((meancount/2.0) - (int) (meancount/2.0) == 0) {
	sptr->flux = (medarray[(meancount/2)-1] + medarray[(meancount/2)])/2.0;
      }
      else {
	sptr->flux = medarray[(meancount-1)/2];
      }
    }
      
    /*
     * Get the estimate on the uncertainty on the flux density.  Set this
     *  to the rms scatter about the median.
     */

    rms = 0.0;
    for(j=0,rptr=raw; j<nraw; j++,rptr++) {
      if(fabs(rptr->day - sptr->day) < width/2.0 && rptr->match >= 0) {
	rms += (rptr->flux - sptr->flux) * (rptr->flux - sptr->flux);
      }
    }
      
    sptr->err = rms/(meancount * (meancount - 1));
  }

  /*
   * Clean up and exit
   */

  printf("  medsmooth: There are %d smooth points that match raw points.\n",
	 nmatch);
  medarray = del_array(medarray);

  return 0;
}

/*.......................................................................
 *
 * Function triangle
 *
 * Performs triangle smoothing on a light curve, producing an output
 *  curve sampled on a regular grid.  The weighting scheme for the
 *  triangle smoothing is:
 *
 *               1 - (2*abs(day0 - day_i))/width
 *  weight_i = ------------------------------------
 *             sum(1 - (2*abs(day0 - day_i))/width)
 *
 * Inputs: Fluxrec *raw        unsmoothed and ungridded curve
 *         int nraw            number of points in unsmoothed curve
 *         Fluxrec *smooth     container for smoothed values.  The "day"
 *                              members of the smoothed curves are already set.
 *         int nsmooth         number of points in smoothed curve.
 *         float width         width of triangle smoothing
 *         float step          step size of gridded output
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 * v14Jul98 CDF, Include variance weighting in calculating mean by calling
 * v30Aug98 CDF, Now only set match to 1 if the smooth curve point is
 *                within 0.25 days of a raw curve point.
 *               Moved allocation of smoothed curve array up to csmooth.
 */

int triangle(Fluxrec *raw, int nraw, Fluxrec *smooth, int nsmooth, float width, 
	     float step)
{
  int i,j;               /* Looping variables */
  int no_error=1;        /* Flag set to 0 on error */
  int meancount;         /* Number of points used for calculating mean */
  int nmatch=0;          /* Number of smooth points that match raw points */
  float sig2;            /* Variance in flux density */
  float *pweight=NULL;   /* Weighting for points used in smoothing */
  float *pptr;           /* Pointer to navigate pweight */
  float wtsum;           /* Sum of weights used for each triangle */
  Fluxrec *temp=NULL;    /* Container for points used to calculate mean */
  Fluxrec *rptr;         /* Pointer for navigating raw lightcurve */
  Fluxrec *sptr;         /* Pointer for navigating smoothed lightcurve */
  Fluxrec *tptr;         /* Pointer to navigate temp */

  /*
   * Allocate container for points to calculate mean.  There will never be
   *  more points in this container than in the raw light curve, so
   *  set it to nraw and then just use the first meancount values.
   */

  if(!(temp = new_fluxrec(nraw)))
    no_error = 0;

  /*
   * Allocate memory for smoothing function weights.  Also set the size of 
   *  this array to the maximum possible value, which is nraw.
   */

  if(!(pweight = new_array(nraw,1)))
    no_error = 0;

  /*
   * Step through the smoothed curve, calculating mean and error in the
   *  mean at each step.  For the smoothing weights at each point, use the
   *  formula in the function descriptor.  Use these weights plus the
   *  variance weighting in calculating the uncertainties in the weighted mean.
   */

  if(no_error) {
    for(i=0,sptr=smooth; i<nsmooth; i++,sptr++) {

      /*
       * Initialize pointers and running sums
       */

      meancount = 0;
      wtsum = 0.0;
      tptr = temp;
      pptr = pweight;

      /*
       * Get the points used to calculate the mean value
       */

      for(j=0,rptr=raw; j<nraw; j++,rptr++) {
	if(fabs(rptr->day - sptr->day) < width/2.0 && rptr->match >= 0) {

	  /*
	   * Calculate the relative weighting for this point
	   */

	  sig2 = rptr->err * rptr->err;
	  *pptr = 1 - (2.0*(fabs(sptr->day - rptr->day))/width);
	  wtsum += *pptr / sig2;
	  sptr->flux += *pptr * rptr->flux / sig2;
	  sptr->err += sig2;
	  *tptr = *rptr;
	  meancount++;
	  pptr++;
	  tptr++;

	  /*
	   * Set match = 1 if sptr is within one fourth of a day of a
	   *  day on the raw curve
	   */

	  if(fabs(sptr->day - rptr->day) < 0.25) {
	    sptr->match = 1;
	    nmatch++;
	  }
	}
      }

      /*
       * Now calculate the mean and the uncertainty in the mean at each
       *  point.  There are three possible cases:
       *    1. No points in window -- Set mean and err to mean and err of
       *                              previous point or, if this is the
       *                              first point, set to mean and err of
       *                              first point in raw curve.
       *    2. One point in window -- Set mean and err to mean and err of
       *                              that point.
       *    3. Two or more points  -- Call var_wmean to calculate the variance
       *                               in the weighted mean.
       */

      switch(meancount) {
      case 0:
	if(sptr > smooth) {
	  sptr->flux = (sptr-1)->flux;
	  sptr->err = (sptr-1)->err;
	}
	else {
	  sptr->flux = raw->flux;
	  sptr->err = raw->err;
	}
	break;
      case 1:
	sptr->flux /= wtsum;
	sptr->err = sqrt(sptr->err);
	break;
      default:
	sptr->flux /= wtsum;
	
	/*
	 * Calculate sigma of the weighted mean
	 */
	  
	if(var_wmean(temp,pweight,meancount,&sptr->err)) {
	  no_error = 0;
	  break;
	}
      }
    }
  }

  /*
   * Clean up and exit
   */

  printf("  triangle: There are %d smooth points that match raw points.\n",
	 nmatch);
  temp = del_fluxrec(temp);
  pweight = del_array(pweight);

  if(no_error)
    return 0;
  else {
    fprintf(stderr,"ERROR: triangle\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function gaussian
 *
 * Performs smoothing on a light curve, weighting each point that contributes
 *  to the smoothing with a gaussian weighting funcion.  The output
 *  curve is sampled on a regular grid.  The weighting scheme for the
 *  gaussian smoothing is:
 *
 *                    1             -x^2 / (2 * sigma^2)
 *  weight_i = ----------------  * e
 *             sigma*sqrt(2 pi)
 *
 * Inputs: Fluxrec *raw        unsmoothed and ungridded curve
 *         int nraw            number of points in unsmoothed curve
 *         Fluxrec *smooth     container for smoothed values.  The "day"
 *                              members of the smoothed curves are already set.
 *         int nsmooth         number of points in smoothed curve.
 *         float sigma         1/e width of smoothing
 *         float step          step size of gridded output
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 * v14Jul98 CDF, Include variance weighting in calculating mean by calling
 *                var_wmean
 * v30Aug98 CDF, Now only set match to 1 if the smooth curve point is
 *                within 0.25 days of a raw curve point.
 *               Moved allocation of smoothed curve array up to csmooth.
 */

int gaussian(Fluxrec *raw, int nraw, Fluxrec *smooth, int nsmooth, float sigma, 
	     float step)
{
  int i,j;               /* Looping variables */
  int no_error=1;        /* Flag set to 0 on error */
  int meancount;         /* Number of points used for calculating mean */
  int nmatch=0;          /* Number of smooth points that match raw points */
  float xdiff;           /* Distance between smooth point and raw point */
  float gnorm;           /* Gaussian normalization factor */
  float *pweight=NULL;   /* Weighting for points used in smoothing */
  float *pptr;           /* Pointer to navigate pweight */
  float sig2;            /* Variance of an individual point */
  float wtsum;           /* Sum of weights used for each triangle */
  Fluxrec *temp=NULL;    /* Container for points used to calculate mean */
  Fluxrec *rptr;         /* Pointer for navigating raw lightcurve */
  Fluxrec *sptr;         /* Pointer for navigating smoothed lightcurve */
  Fluxrec *tptr;         /* Pointer to navigate temp */

  /*
   * Allocate container for points to calculate mean.  There will never be
   *  more points in this container than in the raw light curve, so
   *  set it to nraw and then just use the first meancount values.
   */

  if(!(temp = new_fluxrec(nraw)))
    no_error = 0;

  /*
   * Allocate memory for smoothing function weights.  Also set the size of 
   *  this array to the maximum possible value, which is nraw.
   */

  if(!(pweight = new_array(nraw,1)))
    no_error = 0;

  /*
   * Set the normalization for the gaussian weighting
   */
   
  gnorm = 1.0 / (sigma * sqrt(2*PI));

  /*
   * Step through the smoothed curve, calculating mean and error in the
   *  mean at each step.  For the weights at each point, use the
   *  formula in the function descriptor.  Use these weights when
   *  adding the variances as well, to get an idea of the error in
   *  the weighted mean.
   *
   */

  if(no_error) {
    for(i=0,sptr=smooth; i<nsmooth; i++,sptr++) {

      /*
       * Initialize pointers and running sums
       */

      meancount = 0;
      wtsum = 0.0;
      pptr = pweight;
      tptr = temp;

      /*
       * Get the mean value -- NB: All good points can be used for the
       *  gaussian-weighted mean.
       */

      for(j=0,rptr=raw; j<nraw; j++,rptr++) {
	if(rptr->match >= 0) {

	  /*
	   * Calculate the relative weighting for this point
	   */

	  *tptr = *rptr;
	  xdiff = rptr->day - sptr->day;
	  sig2 = rptr->err * rptr->err;
	  *pptr = gnorm * exp(-xdiff*xdiff/(2*sigma*sigma));
	  wtsum += *pptr / sig2;
	  sptr->flux += *pptr * rptr->flux / sig2;
	  sptr->err += rptr->err * rptr->err;
	  meancount++;
	  pptr++;
	  tptr++;

	  /*
	   * Set match = 1 if sptr is within one fourth of a day of a
	   *  day on the raw curve
	   */

	  if(fabs(sptr->day - rptr->day) < 0.25) {
	    sptr->match = 1;
	    nmatch++;
	  }
	}
      }
      /*
       * Now calculate the mean and the uncertainty in the mean at each
       *  point.  There are three possible cases:
       *    1. No points in window -- Set mean and err to mean and err of
       *                              previous point or, if this is the
       *                              first point, set to mean and err of
       *                              first point in raw curve.
       *    2. One point in window -- Set mean and err to mean and err of
       *                              that point.
       *    3. Two or more points  -- Call var_wmean to calculate the variance
       *                               in the weighted mean.
       */

      switch(meancount) {
      case 0:
	if(sptr > smooth) {
	  sptr->flux = (sptr-1)->flux;
	  sptr->err = (sptr-1)->err;
	}
	else {
	  sptr->flux = raw->flux;
	  sptr->err = raw->err;
	}
	break;
      case 1:
	sptr->flux /= wtsum;
	sptr->err = sqrt(sptr->err);
	break;
      default:
	sptr->flux /= wtsum;
	
	/*
	 * Calculate sigma of the weighted mean
	 */
	  
	if(var_wmean(temp,pweight,meancount,&sptr->err)) {
	  no_error = 0;
	  break;
	}
      }
    }
  }

  /*
   * Clean up and exit
   */

  printf("  gaussian: There are %d smooth points that match raw points.\n",
	 nmatch);
  temp = del_fluxrec(temp);
  pweight = del_array(pweight);

  if(no_error)
    return 0;
  else {
    fprintf(stderr,"ERROR: gaussian\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function varbox
 *
 * Performs a variable width boxcar smoothing on a light curve, producing 
 *  an output curve sampled on a regular grid.  The width of the boxcar
 *  at any point is determined by the density of points in the observed 
 *  light curve.  The function is passed the number of points to be included
 *  in each boxcar and it will change the width accordingly.
 *
 * Inputs: Fluxrec *raw        unsmoothed and ungridded curve
 *         int nraw            number of points in unsmoothed curve
 *         Fluxrec *smooth     container for smoothed values.  The "day"
 *                              members of the smoothed curves are already set.
 *         int nsmooth         number of points in smoothed curve.
 *         int nvar            number of points in the boxcar
 *         float step          step size of gridded output
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 * v30Aug98 CDF, Now only set match to 1 if the smooth curve point is
 *                within 0.25 days of a raw curve point.
 *               Moved allocation of smoothed curve array up to csmooth.
 */

int varbox(Fluxrec *raw, int nraw, Fluxrec *smooth, int nsmooth, int nvar, 
	   float step)
{
  int i;                 /* Looping variable */
  int no_error=1;        /* Flag set to 0 on error */
  int meancount;         /* Number of points used for calculating mean */
  int addcount;          /* Number of points used in current sum of fluxes */
  int nmatch=0;          /* Number of smooth points that match raw points */
  float mindist;         /* Distance of closest observed point */
  Fluxrec *rptr;         /* Pointer for navigating raw lightcurve */
  Fluxrec *sptr;         /* Pointer for navigating smoothed lightcurve */
  Fluxrec *centptr;      /* Pointer to the observed point at the box center */
  Fluxrec *tmpcent;      /* Temporary value for centptr */

  /*
   * Initialize the pointer for the center of the boxcar to be pointing
   *  at the first point in the observed curve.  Because the curves are
   *  arranged in chronological order, the position of the boxcar center
   *  will monotonically increase as we move through the curve.  We
   *  can take advantage of this to reduce the number of iterations
   *  through the observed curve
   */

  centptr = raw;

  /*
   * Step through the smoothed curve, calculating mean and error in the
   *  mean at each step.  For each point, find the nearest observed point
   *  and then take the (nvar-1)/2 observed points on either side of that 
   *  to contribute to the sum.
   *
   */

  for(i=0,sptr=smooth; i<nsmooth; i++,sptr++) {
    meancount = 0;
    mindist = 9999.0;

    /*
     * Step through the observed curve, starting from the current value
     *  of centptr, to find the observed point closest to the current
     *  point in the smoothed curve.
     */

    for(rptr=centptr; rptr<(raw+nraw-1); rptr++) {

      /*
       * Get central point
       */

      if(fabs(rptr->day - sptr->day) < mindist && rptr->match >= 0) {
	mindist = fabs(rptr->day - sptr->day);
	tmpcent = rptr;
	if(fabs(rptr->day - sptr->day) < 0.25) {
	  sptr->match = 1;
	  nmatch++;
	}
      }
    }

    /*
     * Get contribution from central point
     */

    centptr = tmpcent;
    sptr->flux += centptr->flux;
    sptr->err += centptr->err * centptr->err;
    meancount++;

    /*
     * Get contribution of (nvar-1)/2 points below centptr, if they exist
     */

    addcount = 0;
    rptr = centptr;
    while(addcount < (nvar-1)/2 && rptr >= raw+1) {
      rptr--;
      if(rptr->match >= 0) {
	addcount++;
	meancount++;
	sptr->flux += rptr->flux;
	sptr->err += rptr->err * rptr->err;
      }
    }

    /*
     * Get contribution of (nvar-1)/2 points above centptr, if they exist
     */

    addcount = 0;
    rptr = centptr;
    while(addcount < (nvar-1)/2 && rptr <= raw+nraw-2) {
      rptr++;
      if(rptr->match >= 0) {
	addcount++;
	meancount++;
	sptr->flux += rptr->flux;
	sptr->err += rptr->err * rptr->err;
      }
    }

    /*
     * Calculate the mean value, doing some error checking along the
     *  way.  Since this is a boxcar smoothing, the mean value should
     *  just be sptr->flux / meancount.
     */

    if(meancount == 0) {
      fprintf(stderr,"ERROR: varbox.  No valid points in smoothing box.\n");
      no_error = 0;
    } else if(meancount > nvar) {
      fprintf(stderr,"ERROR: varbox.  meancount (%d) > nvar (%d) in window.\n",
	      meancount,nvar);
      no_error = 0;
    } else {
      sptr->flux /= meancount;
      sptr->err = sqrt(sptr->err / meancount);
    }
  }

  printf("  varbox: There are %d smooth points that match raw points.\n",
	 nmatch);

  if(no_error)
    return 0;
  else {
    fprintf(stderr,"ERROR: varbox\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function vartri
 *
 * Performs a variable-width triangle smoothing on a light curve, producing 
 *  an output curve sampled on a regular grid.  The width of the window
 *  at any point is determined by the density of points in the observed 
 *  light curve.  The function is passed the number of points to be included
 *  in each window and it will change the width accordingly.
 *
 * Inputs: Fluxrec *raw        unsmoothed and ungridded curve
 *         int nraw            number of points in unsmoothed curve
 *         Fluxrec *smooth     container for smoothed values.  The "day"
 *                              members of the smoothed curves are already set.
 *         int nsmooth         number of points in smoothed curve.
 *         int nvar            number of points in the boxcar
 *         float step          step size of gridded output
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 * v30Aug98 CDF, Now only set match to 1 if the smooth curve point is
 *                within 0.25 days of a raw curve point.
 *               Moved allocation of smoothed curve array up to csmooth.
 */

int vartri(Fluxrec *raw, int nraw, Fluxrec *smooth, int nsmooth, int nvar, 
	   float step)
{
  int i,j;               /* Looping variables */
  int no_error=1;        /* Flag set to 0 on error */
  int meancount;         /* Number of points used for calculating mean */
  int addcount;          /* Number of points used in current sum of fluxes */
  int nmatch=0;          /* Number of smooth points that match raw points */
  float mindist;         /* Distance of closest observed point */
  float maxdist;         /* Half of width of triangle */
  float ptweight;        /* Weighting for a point used in smoothing */
  float wtsum;           /* Sum of weights used for each triangle */
  Fluxrec *tmpsmo=NULL;  /* Container for points used to calculate mean */
  Fluxrec *rptr;         /* Pointer for navigating raw lightcurve */
  Fluxrec *sptr;         /* Pointer for navigating smoothed lightcurve */
  Fluxrec *tmpptr;       /* Pointer for navigating tmpsmo */
  Fluxrec *centptr;      /* Pointer to the observed point at the box center */
  Fluxrec *tmpcent;      /* Temporary value for centptr */

  /*
   * Allocate an array to hold the points used to calculated the weighted
   *  mean at each step.  This array can have a maximum number of members
   *  equal to nvar.
   */

  if(!(tmpsmo = new_fluxrec(nvar))) 
    no_error = 0;

  /*
   * Initialize the pointer for the center of the boxcar to be pointing
   *  at the first point in the observed curve.  Because the curves are
   *  arranged in chronological order, the position of the boxcar center
   *  will monotonically increase as we move through the curve.  We
   *  can take advantage of this to reduce the number of iterations
   *  through the observed curve
   */

  centptr = raw;

  /*
   * Step through the smoothed curve, calculating mean and error in the
   *  mean at each step.  For each point, find the nearest observed point
   *  and then take the (nvar-1)/2 observed points on either side of that 
   *  to contribute to the sum.
   *
   */

  for(i=0,sptr=smooth; i<nsmooth; i++,sptr++) {
    meancount = 0;
    mindist = 9999.0;
    tmpptr = tmpsmo;

    /*
     * Step through the observed curve, starting from the current value
     *  of centptr, to find the observed point closest to the current
     *  point in the smoothed curve.
     */

    for(rptr=centptr; rptr<(raw+nraw-1); rptr++) {

      /*
       * Get central point
       */

      if(fabs(rptr->day - sptr->day) < mindist && rptr->match >= 0) {
	mindist = fabs(rptr->day - sptr->day);
	tmpcent = rptr;
	*tmpsmo = *rptr;
	if(fabs(rptr->day - sptr->day) < 0.25) {
	  sptr->match = 1;
	  nmatch++;
	}
      }
    }

    /*
     * Put central point into tmpsmo
     */

    centptr = tmpcent;
    meancount++;
    tmpptr++;

    /*
     * Put (nvar-1)/2 points below centptr into tmpsmo, if they exist
     */

    addcount = 0;
    rptr = centptr;
    while(addcount < (nvar-1)/2 && rptr >= raw+1) {
      rptr--;
      if(rptr->match >= 0) {
	addcount++;
	meancount++;
	*tmpptr = *rptr;
	tmpptr++;
      }
    }

    /*
     * Put (nvar-1)/2 points above centptr into tmpsmo, if they exist
     */

    addcount = 0;
    rptr = centptr;
    while(addcount < (nvar-1)/2 && rptr <= raw+nraw-2) {
      rptr++;
      if(rptr->match >= 0) {
	addcount++;
	meancount++;
	*tmpptr = *rptr;
	tmpptr++;
      }
    }

    /*
     * Calculate the maximum displacement from sptr->day and set the
     *  triangle size scale to twice that value.
     */

    maxdist = 0.0;
    for(j=0,tmpptr=tmpsmo; j<meancount; j++,tmpptr++)
      if(fabs(tmpptr->day - sptr->day) > maxdist)
	maxdist = fabs(tmpptr->day - sptr->day);
    maxdist *= 2.0;

    /*
     * Calculate the mean value using the triangle weighting scheme, doing 
     *  some error checking along the way.
     */

    if(meancount == 0) {
      fprintf(stderr,"ERROR: vartri.  No valid points in smoothing box.\n");
      no_error = 0;
    } else if(meancount > nvar) {
      fprintf(stderr,"ERROR: vartri.  meancount (%d) > nvar (%d) in ",
	      meancount,nvar);
      fprintf(stderr,"smoothing box.\n");
      no_error = 0;
    } else {

      /*
       * Get the mean value
       */

      wtsum = 0.0;
      for(j=0,tmpptr=tmpsmo; j<meancount; j++,tmpptr++) {

	/*
	 * Calculate the relative weighting for this point
	 */

	ptweight = 1 - fabs(sptr->day - tmpptr->day)/maxdist;
	wtsum += ptweight;
	sptr->flux += ptweight * tmpptr->flux;
	sptr->err += ptweight * tmpptr->err * tmpptr->err;
      }

      /*
       * Set mean = (weighted sum)/sum(weights) and error on estimate of mean 
       *  to sqrt(weighted sum of variances / sum(weights).  
       */

      sptr->flux /= wtsum;
      sptr->err = sqrt(sptr->err / wtsum);
    }
  }

  /*
   * Clean up and exit
   */

  printf("  vartri: There are %d smooth points that match raw points.\n",
	 nmatch);
  tmpsmo = del_fluxrec(tmpsmo);

  if(no_error)
    return 0;
  else {
    fprintf(stderr,"ERROR: vartri\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function var_wmean
 *
 * Calculates the one-sigma uncertainty on a weighted mean including
 *  effects from both the weighting function and the statistical
 *  (variance) weighting.  If each point x_i in the calculation of the
 *  weighted mean has a weight of w_i = f_i / (sigma_i)^2, where 
 *  f_i is the weight from the smoothing function, then the variance
 *  of the weighted mean is (see Book 6, p. 53):
 *
 *                        sum (f_i / sigma_i)^2
 *  (sigma_wmean)^2 = -----------------------------
 *                    ( sum (f_i / (sigma_i)^2) )^2
 *
 * Inputs: Fluxrec *flux       points used in calculating weighted mean,
 *                              including sigma_i information
 *         float *pweight      weighting (smoothing function only) for each
 *                              point
 *         int npoints         number of points used in calculating mean
 *         float *sig_wmean    uncertainty in weighted mean (set by this
 *                              function)
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> errors
 *
 */

int var_wmean(Fluxrec *flux, float *pweight, int npoints, float *sig_wmean)
{
  int i;             /* Looping variable */
  float *wptr;       /* Pointer to navigate pweight */
  float num=0.0;     /* Numerator in above equation */
  float den=0.0;     /* Denominator in above equation */
  float sig2;        /* Variance of points in flux */
  Fluxrec *fptr;     /* Pointer to navigate flux */

  /*
   * Loop through the points, calculating the numerator and denominator in
   *  the above equation.
   */

  for(i=0,wptr=pweight,fptr=flux; i<npoints; i++,wptr++,fptr++) {
    sig2 = fptr->err * fptr->err;
    num += (*wptr * *wptr / sig2);
    den += (*wptr / sig2) * (*wptr / sig2);
  }

  /*
   * Set sig_wmean to the square root of var_wmean
   */

  *sig_wmean = sqrt(num/den);

  return 0;
}

