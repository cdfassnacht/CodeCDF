/*
 * smooth_curve.c
 *
 * Usage: smooth_curve [input filename] [output filename] ([setup file])
 *
 * Reads in a light curve and smooths it according to the data in the
 *  setup file or from data entered interactively.
 *
 * 28Jun98 CDF, Basically a modification of lcurve.c
 * v24Jul98 CDF, Replaced calls to smoothing functions by a call to csmooth
 *                (which is in lc_interp.c).
 * v30Aug98 CDF, Added smooth-in-place option.
 *               Added plotting of raw and smoothed curves.
 * v01Sep98 CDF, Added bad day flagging.
 * v25Sep98 CDF, Modified smoothing to call functions in lc_interp.c rather
 *                than doing smoothing in the main program.
 * v12Sep00 CDF, Changed some function calls in response to changes in
 *                lc_interp.c and older changes in lc_setup.c.
 * v04Oct00 CDF, Got rid of bad-day flagging option, since any bad-day
 *                flagging should take place in fluxplot.sm.
 *               Added a necessary call to set_grid_params.
 * v19Feb02 CDF, Added plotting of residuals if smooth-in-place option
 *                was chosen.
 * v03Apr02 CDF, Moved writing of output files into new write_fluxrec
 *                function in lc_funcs.c
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "structdef.h"
#include "dataio.h"
#include "plotfuncs.h"
#include "cpgplot.h"
#include "nr.h"
#include "lc_funcs.h"
#include "lc_interp.h"

/*.......................................................................
 *
 * Function declarations
 *
 */

int plot_points(Fluxrec *flux,  int npoints, char *xlab, char *ylab, 
		char *title);

/*.......................................................................
 *
 * Main program
 *
 */

int main(int argc, char *argv[])
{
  int i;                /* Looping variable */
  int no_error=1;       /* Flag set to 0 on error */
  int nlines;           /* Number of data lines in input file */
  float matchtol=0.04;  /* Tolerance for determining if days are the same */
  char infile[MAXC];    /* Input file name */
  char outfile[MAXC];   /* Output file name */
  char setupfile[MAXC]; /* Name of setup file */
  Fluxrec *flux=NULL;   /* Input light curve */
  Fluxrec *diff=NULL;   /* Difference: smoothed - unsmoothed */
  Fluxrec *smooth=NULL; /* Smoothed light curve */
  Fluxrec *sptr;        /* Pointer to navigate smooth */
  Fluxrec *fptr;        /* Pointer to navigate flux */
  Fluxrec *dptr;        /* Pointer to navigate diff */
  Setup *setup=NULL;    /* Container for setup information */

  /*
   * Check input line
   */

  if(argc < 3) {
    fprintf(stderr,"\nUsage: smooth_curve [input filename] [output filename]");
    fprintf(stderr," ([setup filename])\n");
    fprintf(stderr," Where the setup file is optional.\n\n");
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
    if(setup_interp(setup))
      no_error = 0;

  /*
   * Summarize setup parameters
   */

  if(no_error)
    setup_interp_summary(setup);

  /*
   * Read data from input file
   */

  if(no_error) {
    strcpy(infile,argv[1]);
    if(!(flux = read_fluxrec(infile,'#',&nlines)))
      no_error = 0;
    else
      printf("%d data lines in input file.\n\n",nlines);
  }

  /*
   * Set parameters for interpolated curves.
   */

  if(no_error) {
    set_grid_params(setup,flux,nlines);
    printf("\n");
  }

  /*
   * Smooth or linearly interpolate the flattened curves, depending on
   *  flag in setup.
   */

  if(setup->dosmooth != NOSMOOTH && no_error)
    if(!(smooth = interp_switch(flux,nlines,setup,0)))
	no_error = 0;

  /*
   * Print out curve to output file.
   */


  if(no_error) {
    printf("\n");
    strcpy(outfile,argv[2]);
    if(write_fluxrec(smooth,setup->ninterp,outfile,1,matchtol))
      no_error = 0;
  }

  /*
   * Print out number of days discarded during the printing to file phase
   */

  /*
   * If the smooth-in-place option was chosen, calculate the difference
   *  between the smoothed and unsmoothed light curves.
   */

  if(setup->dosmooth == SMINPLACE && no_error) {
    if(!(diff = new_fluxrec(nlines)))
      no_error = 0;
    else {
      for(i=0,fptr=flux,sptr=smooth,dptr=diff; i<nlines; 
	  i++,fptr++,sptr++,dptr++) {
	*dptr = *fptr;
	dptr->flux = fptr->flux - sptr->flux;
      }
    }
    printf("\nWriting difference data to diff.dat.\n");
    if(write_fluxrec(diff,nlines,"diff.dat",0,0.0))
      no_error = 0;
  }

  /*
   * Open the plot
   */

  printf("\n");
  if(cpgbeg(0, "?", 1, 2) != 1)
    return 1;
  cpgpap(0.0,1.0);
  cpgsch(1.5);
  cpgvstd();
  cpgscf(2);

  /*
   * Plot the raw and smoothed curves
   */

  if(no_error)
    if(plot_points(flux,nlines,"Day","Flux Density","Raw Curve"))
      no_error = 0;

  
  if(no_error)
    if(plot_points(smooth,setup->ninterp,"Day","Flux Density",
		   "Smoothed Curve"))
      no_error = 0;

  /*
   * If the smooth-in-place option was chosen, plot out smoothed
   *  curve and difference curve.
   */

  if(setup->dosmooth == SMINPLACE && no_error) {
    if(plot_points(smooth,setup->ninterp,"Day","Flux Density",
		   "Smoothed Curve"))
      no_error = 0;

    if(no_error)
      if(plot_points(diff,nlines,"Day","Residual Flux Density","Residuals"))
	no_error = 0;
  }

  /*
   * Clean up and exit
   */

  if(no_error)
    printf("\nCleaning up\n");

  plot_close();
  flux = del_fluxrec(flux);
  diff = del_fluxrec(diff);
  smooth = del_fluxrec(smooth);
  setup = del_setup(setup);

  if(no_error) {
    printf("Program finished.\n\n");
    return 0;
  }
  else {
    fprintf(stderr,"ERROR: Exiting program\n\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function plot_points
 *
 * Given input lightcurve, draws a box and then plots
 *  points with errorbars.
 *
 * Inputs: Fluxrec *flux       input lightcurve
 *         int npoints         number of points in curve
 *         char *xlab          label for x axis
 *         char *ylab          label for y axis
 *         char *title         title
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int plot_points(Fluxrec *flux,  int npoints, char *xlab, char *ylab, 
		char *title)
{
  int i;             /* Looping variable */
  float xmin,xmax;   /* Min and max x values */
  float ymin,ymax;   /* Min and max y values */
  float size;        /* Temporary difference between min and max values */
  Fluxrec *fptr;     /* Pointer for navigating flux */

  /*
   * Initialize
   */

  xmin = xmax = flux->day;
  ymin = ymax = flux->flux;

  /*
   * Find min and max values
   */

  for(i=0,fptr=flux; i<npoints; i++,fptr++) {
    if(fptr->day > xmax)
      xmax = fptr->day;
    if(fptr->day < xmin)
      xmin = fptr->day;
    if(fptr->flux > ymax)
      ymax = fptr->flux;
    if(fptr->flux < ymin)
      ymin = fptr->flux;
  }

  /*
   * Increase size of box a bit
   */

  size = xmax - xmin;
  xmin -= 0.15*size;
  xmax += 0.15*size;
  size = ymax - ymin;
  ymin -= 0.15*size;
  ymax += 0.15*size;

  /*
   * Draw box and label axes
   */

  cpgenv(xmin,xmax,ymin,ymax,0,0);
  cpglab(xlab,ylab,title);

  /*
   * Put in the points
   */

  for(i=0,fptr=flux; i<npoints; i++,fptr++) {
    ymin = fptr->flux - fptr->err;
    ymax = fptr->flux + fptr->err;
    cpgpt(1,&fptr->day,&fptr->flux,16);
    cpgerry(1,&fptr->day,&ymax,&ymin,1.0);
  }

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
