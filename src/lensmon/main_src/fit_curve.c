/*
 * fit_curve.c
 *
 * Usage: fit_curve [input_filename] [degree of polynomial fit]
 *
 * Reads in a data set consisting of 3 columns - x, y, sigma and fits
 *  a polynomial to the data points.
 *
 * 15Aug98 CDF
 * v16Aug98 CDF, Moved plotting into new functions plot_points and 
 *                plot_fit_curve.
 *               Moved the chisq calculation function (chisq_fit) into
 *                lc_funcs.c
 *               Added print_fit_info to print parameters to output file.
 * v30Sep98 CDF, Moved plot_points into plotfuncs.c where it is now called 
 *                plot_lcurve.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "structdef.h"
#include "plotfuncs.h"
#include "cpgplot.h"
#include "lc_funcs.h"

/*.......................................................................
 *
 * Function declarations
 *
 */

int plot_fit_curve(Fluxrec *flux, int npoints, float *a, int ncoeff);
void print_fit_info(float *a, int ncoeff);

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
  int ndeg;            /* Degree of polynomial fit */
  float xmid;          /* Value of x in middle of light curve */
  float ymean;         /* Mean value of y */
  float yrms;          /* RMS scatter about ymean */
  float *a=NULL;       /* Container for polynomial coefficients */
  float chisq=0.0;     /* Chi-squared for the fit */
  char inname[MAXC];   /* Name of input file */
  Fluxrec *flux=NULL;  /* Light curve to which polynomial is to be fit */
  Fluxrec *fptr;       /* Pointer to navigate flux */

  /*
   * Check command line
   */

  if(argc != 3) {
    fprintf(stderr,"\nUsage: fit_curve [input_filename] ");
    fprintf(stderr,"[degree of polynomial fit].\n\n");
    return 1;
  }

  /*
   * Read degree of fit from command line
   */

  if(sscanf(argv[2],"%d",&ndeg) != 1 && ndeg < 1) {
    fprintf(stderr,"\nERROR: Bad input.  Exiting program.\n\n");
    return 1;
  }

  /*
   * Get input file name from command line in case it needs to be
   *  modified later
   */

  strcpy(inname,argv[1]);

  /*
   * Read input file
   */

  if(!(flux = read_fluxrec(inname,'#',&npoints))) {
    fprintf(stderr,"ERROR.  Exiting program\n");
    return 1;
  }

  /*
   * Find xmid and ymean
   */

  xmid = flux->day + ((flux+npoints-1)->day - flux->day)/2.0;
  if(calc_mean(flux,npoints,&ymean,&yrms))
    no_error = 0;
  else
    printf("xmid = %6.2f, ymean = %6.3f\n",xmid,ymean);

  /*
   * Convert light curve to a light curve centered at xmid
   */

  if(no_error) {
    for(i=0,fptr=flux; i<npoints; i++,fptr++) {
#if 0
      fptr->day -= xmid;
      fptr->flux = fptr->flux / ymean - 1.0;
      fptr->err = fptr->err / ymean;
#endif
    }
  }

  /*
   * Call polynomial fitting function
   */

  if(no_error)
    if(!(a = fit_poly(flux,npoints,ndeg+1,1)))
      no_error = 0;

  /*
   * Calculate chisq for the fit
   */

  if(no_error)
    if((chisq = chisq_fit(flux,npoints,a,ndeg+1,0,1)) < 0.0)
      no_error = 0;

  /*
   * Open the plot
   */

  plot_open(1.5);

  /*
   * Plot points
   */

  if(no_error)
    if(plot_lcurve(flux,npoints,"Days","Flux Density",""))
      no_error = 0;

  /*
   * Plot fitted polynomial and close plot
   */

  if(no_error)
    if(plot_fit_curve(flux,npoints,a,ndeg+1))
      no_error = 0;

  plot_close();

  /*
   * Print important info to output file
   */

  if(no_error)
    print_fit_info(a,ndeg+1);

  /*
   * Clean up and exit
   */

  flux = del_fluxrec(flux);
  a = del_array(a);

  if(no_error)
    return 0;
  else {
    fprintf(stderr,"ERROR: Exiting program\n\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function plot_fit_curve
 *
 * Plots the polynomial that has been fit to the light curve
 *
 * Inputs: Fluxrec *flux       input light curve
 *         int npoints         number of points in light curve
 *         float *a            polynomial coefficients
 *         int ncoeff          number of coefficients
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int plot_fit_curve(Fluxrec *flux, int npoints, float *a, int ncoeff)
{
  int i,j;            /* Looping variables */
  int nstep=500;      /* Number of steps used to compute curve */
  float dmin,dmax;    /* Min and max values of day, in case flux is not sorted */
  float xtmp,ytmp;    /* Values of x,y used in constructing the curve */
  float *polyx=NULL;  /* Polynomial in x */
  Fluxrec *fptr;      /* Pointer to navigate flux */

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
   * Loop through curve, finding min and max values of day
   */

  dmin = dmax = flux->day;
  for(i=0,fptr=flux; i<npoints; i++,fptr++) {
    if(fptr->day > dmax)
      dmax = fptr->day;
    if(fptr->day < dmin)
      dmin = fptr->day;
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

    xtmp = flux->day + i*(dmax - dmin)/(1.0 * nstep);
    poly(xtmp,polyx-1,ncoeff);
    ytmp = 0;
    for(j=0; j<ncoeff; j++)
      ytmp += a[j] * polyx[j];
    printf("%f %f\n",xtmp,ytmp);

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
 * Function print_fit_info
 *
 * Prints components of fit to an output file.
 *
 * Inputs: float *a            coefficients of fit
 *         int ncoeff          number of coefficients
 *
 * Output: none
 *
 */

void print_fit_info(float *a, int ncoeff)
{
  int i;
  char outname[MAXC];  /* Output file name */
  FILE *ofp=NULL;      /* Output file pointer */

  /*
   * Get output file name and open file
   */

  printf("Enter output file name: [no output file] ");
  gets(outname);
  if(strcmp(outname,"") == 0) {
    return;
  }
  while(!(ofp = fopen(outname,"w"))) {
    fprintf(stderr,"ERROR: print_fit_info.  %s cannot be opened.\n",outname);
    printf("Enter new filename:  ");
    gets(outname);
  }

  /*
   * Print output info
   */

  for(i=0; i<ncoeff; i++)
    fprintf(ofp,"a_%d %14.10f\n",i,a[i]);

  /*
   * Clean up
   */

  if(ofp)
    fclose(ofp);
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
