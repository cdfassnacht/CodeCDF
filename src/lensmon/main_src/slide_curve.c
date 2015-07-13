/*
 * slide_curve.c
 *
 * Usage: slide_curve [lightcurve filename] [polynomial filename]
 *
 * Reads in a light curve and a polynomial function defined by its
 *  coefficients.  The program then slides the polynomial, which is defined 
 *  in the range 0 +/- MAXDAY days, along the light curve.  At each
 *  step, the points in the light curve within the +/- MAXDAY range are
 *  used to compute a reduced chisq.  The step at which the smallest
 *  chisq is computed is recorded.
 *
 * 16Aug98 CDF
 * v17Aug98 CDF, Changed limits on x in step_curve.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "structdef.h"
#include "dataio.h"
#include "lc_funcs.h"

#define MAXDAY 69.0

/*.......................................................................
 *
 * Function declarations
 *
 */

float *read_poly(char *polyname, int *ncoeff);
int step_curve(Fluxrec *flux, int npoints, float *a, int ncoeff);

/*.......................................................................
 *
 * Main program
 *
 */

int main(int argc, char *argv[])
{ 
  int i;               /* Looping variable */
  int no_error=1;      /* Flag set to 0 on error */
  int npoints;         /* Number of points in light curve */
  int ncoeff;          /* Number of coefficients in polynomial fit */
  float *a=NULL;       /* Container for polynomial coefficients */
  float *polyx=NULL;   /* Polynomial in x */
  char inname[MAXC];   /* Name of light curve file */
  char polyname[MAXC]; /* Name of polynomial file */
  Fluxrec *flux=NULL;  /* Light curve */

  /*
   * Check command line
   */

  if(argc != 3) {
    fprintf(stderr,"\nUsage: slide_curve [light curve filename] ");
    fprintf(stderr,"[polynomial filename].\n\n");
    return 1;
  }

  /*
   * Get lightcurve filename from command line in case it needs to be
   *  modified later
   */

  strcpy(inname,argv[1]);

  /*
   * Read light curve file
   */

  if(!(flux = read_fluxrec(inname,'#',&npoints))) {
    fprintf(stderr,"ERROR.  Exiting program\n");
    return 1;
  }

  /*
   * Get polynomial filename from command line in case it needs to be
   *  modified later
   */

  strcpy(polyname,argv[2]);

  /*
   * Read light curve file
   */

  if(!(a = read_poly(polyname,&ncoeff)))
    no_error = 0;
  else {
    printf("Polynomial coefficients:\n");
    for(i=0; i<ncoeff; i++)
      printf("  a_%d = %14.10f\n",i,a[i]);
    printf("\n");
  }

  /*
   * Call stepping function
   */

  if(no_error)
    if(step_curve(flux,npoints,a,ncoeff))
      no_error = 0;

#if 0
  /*
   * Call polynomial fitting function
   */

  if(no_error)
    if(!(a = fit_poly(flux,npoints,ndeg+1)))
      no_error = 0;

  /*
   * Allocate memory for container for polynomial values
   */

  if(no_error)
    if(!(polyx = new_array(ndeg+1,1)))
      no_error = 0;

  /*
   * Calculate chisq for the fit
   */

  if(no_error)
    if((chisq = chisq_fit(flux,npoints,polyx,a,ndeg+1)) < 0.0)
      no_error = 0;
#endif
  /*
   * Clean up and exit
   */

  flux = del_fluxrec(flux);
  a = del_array(a);
  polyx = del_array(polyx);

  if(no_error)
    return 0;
  else {
    fprintf(stderr,"ERROR: Exiting program\n\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function read_poly
 *
 * Reads the file containing the coefficients of the polynomial fit.
 *
 * Inputs: char *polyname      name of polynomial file
 *         int *ncoeff         number of coefficients (set by this function)
 *
 * Output: float *a            array containing polynomial coefficients
 *
 */

float *read_poly(char *polyname, int *ncoeff)
{
  float *a=NULL;   /* Container for polynomial coefficients */
  float *aptr;     /* Pointer used to navigate a */
  char line[MAXC]; /* String to read input file lines */
  char junk[MAXC]; /* Used to read unwanted info from input file */
  FILE *ifp=NULL;  /* Input file pointer */

  /*
   * Open the input file, if possible
   */

  while(!(ifp = fopen(polyname,"r"))) {
    fprintf(stderr,"ERROR: read_poly.  Cannot open %s.\n",polyname);
    printf("Enter filename again:  ");
    gets(polyname);
  }

  /*
   * Count number of data lines in input file, which should be the
   *  number of coefficients in the polynomial fit.
   */

  if((*ncoeff = n_lines(ifp,'#')) <= 0) {
    fprintf(stderr,"ERROR: read_poly.  No data lines in input file\n");
    return NULL;
  }

  /*
   * Allocate memory for the coefficient array
   */

  if(!(a = new_array(*ncoeff,1))) {
    fprintf(stderr,"ERROR: read_poly.\n");
    return NULL;
  }

  /*
   * Read the coefficients from the input file
   */

  rewind(ifp);
  aptr = a;
  while(fgets(line,MAXC,ifp) != NULL) {
    if(line[0] != '#') {
      if(sscanf(line,"%s %f",junk,aptr) != 2) {
	fprintf(stderr,"ERROR: read poly.  Bad input format in %s.\n",
		polyname);
	return del_array(a);
      }
      else {
	aptr++;
      }
    }
  }

  /*
   * Clean up and exit
   */

  if(ifp)
    fclose(ifp);

  return a;
}

/*.......................................................................
 *
 * Function step_curve
 *
 * Creates an evenly stepped series of points along the time axis of
 *  the light curve.
 *
 * Inputs: Fluxrec *flux       light curve
 *         int npoints         number of points in light curve
 *         float *a            coefficients of polynomial fit
 *         int ncoeff          number of coefficients
 *
 * Output: int (PROBABLY WILL CHANGE)
 *
 */

int step_curve(Fluxrec *flux, int npoints, float *a, int ncoeff)
{
  int i,j,k;           /* Looping variables */
  int no_error=1;      /* Flag set to 0 on error */
  int ncalc;           /* Number of points used in chisq calculation */
  float ymean;         /* Mean value of flux along light curve */
  float yrms;          /* RMS scatter about ymean */
  float yscale;        /* Scale factor for polynomial y values */
  float x0;            /* Current step along light curve */
  float x;             /* Day relative to x0 */
  float y;             /* Value of the polynomial at x */
  float chisq;         /* Reduced chisq at x0 */
  float chimin=999.0;  /* Minimum value of chisq found */
  float bestx,besty;   /* Values of x and y associated with chimin */
  float *polyx=NULL;   /* x values in polynomial */
  Fluxrec *fptr;       /* Pointer to navigate flux */

  /*
   * Allocate memory for polyx
   */

  if(!(polyx = new_array(ncoeff,1))) {
    fprintf(stderr,"ERROR: step_curve.\n");
    return 1;
  }

  /*
   * Find ymean
   */

  if(no_error) {
    if(calc_mean(flux,npoints,&ymean,&yrms))
      no_error = 0;
    else
      printf("ymean = %6.3f\n\n",ymean);
  }

  /*
   * Do outer loop on yscale
   */

  for(k=-50; k<51; k++) {
    
    /*
     * Set normalization factor for this loop
     */

    yscale = ymean * (1.0 + k*FLUXSTEP/5.0);

    /*
     * Set x0 to first point in light curve + MAXDAY
     */

    x0 = flux->day + MAXDAY;

    /*
     * Loop through light curve, increasing x0 by DAYSTEP each time
     */

    while(x0 <= (flux+npoints-1)->day - MAXDAY) {

      /*
       * Initialize chisq and ncalc
       */

      chisq = 0.0;
      ncalc = 0;

      /*
       * Loop over points in the light curve
       */

      for(i=0,fptr=flux; i<npoints; i++,fptr++) {

	/*
	 * Set x to the day relative to x0
	 */

	x = fptr->day - x0;

	/*
	 * Check to see if x is within +/- MAXDAY of x0.  If it is,
	 *  do the calculations.
	 */

	if(fabs(x) < MAXDAY) {

	  /*
	   * Calculate "x" values of the polynomial
	   */

	  poly(x,polyx-1,ncoeff);

	  /*
	   * Multiply x values by coefficients to get y
	   */

	  y = 0;
	  for(j=0; j<ncoeff; j++)
	    y += a[j] * polyx[j];

	  /*
	   * Now multiply y by its proper scale factor
	   */

	  y = yscale * (1.0 + y);

	  /*
	   * Add chisq for this point to running sum
	   */

	  chisq += (fptr->flux - y) * (fptr->flux - y) /
	    (fptr->err * fptr->err);
	  ncalc++;
	}
      }

      /*
       * Convert chisq into reduced chisq
       */

      if(ncalc <= ncoeff) {
	fprintf(stderr,"ERROR: step_curve.\n");
	fprintf(stderr,"  ncalc (%d) < ncoeff (%d) for day = %5.1f\n",
		ncalc,ncoeff,x0);
	fprintf(stderr,"  Setting reduced chisq = 999.0\n");
	chisq = 999.0;
      }
      else {
	chisq /= (ncalc - ncoeff);
	if(chisq < chimin) {
	  chimin = chisq;
	  bestx = x0;
	  besty = yscale;
	}
      }

      /*
       * Increment x0
       */

      x0 += DAYSTEP;
    }
  }

  /*
   * Print out results from fitting
   */

  printf("\nBest fit -- day = %5.1f, scaling %6.3f.  Reduced chisq = %f\n\n",
	 bestx,besty,chimin);

  /*
   * Clean up and return
   */

  polyx = del_array(polyx);

  return 0;
}
