/*
 * add_disp.c
 *
 * Usage: 
 *
 * Takes as input two sets of identically sampled gridpoints in 
 *  (delay,magnification,dispersion) space, as produced by running delays.c
 *  on separate sets of light curves.  The program will then find the
 *  minimum COMBINED dispersion, allowing for shifts in magnification between
 *  seasons, to deal with possible changes in magnifications caused by 
 *  microlensing.  The input data sets have 3 parameter values: delay, 
 *  magnification, and dispersion, and will be contained in the dispb*.dat 
 *  files produced by choosing the dispersion option in delays.c.  
 *
 * 10Oct00 CDF,  A direct modification of add_chisq.c
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
  int i,j,k;                /* Looping variables */
  int no_error=1;           /* Flag set to 0 on error */
  int nfiles;               /* Number of input files */
  int minpos;               /* Array position associated with dispmin */
  int *mu_shift=NULL;       /* Grid shifts between datasets */
  int maxshift=0;           /* Maximum value of mu_shift */
  float dispmin;            /* Minimum combined reduced disp */
  Axisinfo tau,tauchk;      /* Tau axis info (and check) */
  Axisinfo mu,muchk;        /* Mu axis info (and check) */
  char file1[MAXC];         /* Input file name for first file */
  char newfile[MAXC];       /* Input file name for other input file(s) */
  LCdisp mincheck;          /* Container used in min disp search */
  LCdisp **grid={NULL};     /* Grid values */
  LCdisp **mindisp={NULL};  /* Vectors containing minimum disp info */
  LCdisp *finaldisp=NULL;    /* Vectors containing minimum disp info */
  LCdisp *gptr;             /* Pointer for navigating the grids */
  LCdisp *dptr;             /* Pointer for navigating LCdisp arrays */
  FILE *ofp=NULL;           /* Output file pointer for location of minimum */
 
  /*
   * Check input line
   */

  if(argc < 3) {
    fprintf(stderr,"\nUsage: add_disp [file1] [file2] ... [fileN]\n");
    fprintf(stderr," where the input files are disp files produced by");
    fprintf(stderr," separate runs of delays.c.\n\n");
    fprintf(stderr,"*** NB: The input parameter-space grids must have ");
    fprintf(stderr,"identical sizes and spacings.\n\n");
    return 1;
  }
  else {
    nfiles = argc - 1;
    printf("\nRunning add_disp.c with %d input files.\n",nfiles);
  }

  /*
   * Allocate first level of grid pointers and mindisp pointers
   */

  grid = (LCdisp **) malloc(sizeof(LCdisp *) * nfiles);
  if(!grid) {
    fprintf(stderr,"ERROR:  Insufficient memory for grid array.\n");
    return 1;
  }

  mindisp = (LCdisp **) malloc(sizeof(LCdisp *) * nfiles);
  if(!mindisp) {
    fprintf(stderr,"ERROR:  Insufficient memory for mindisp array.\n");
    no_error = 0;
  }

  /*
   * Read data from first input file
   */

  if(no_error) {
    strcpy(file1,argv[1]);
    if(!(grid[0] = read_lcdisp(file1,&tau,&mu))) {
      fprintf(stderr,"ERROR: Exiting add_disp.c.\n");
      return 1;
    }
  }

  /*
   * Read data from other input file(s) and check that grid sizes
   *  match those of first input file.
   */

  if(no_error) {
    for(j=1; j<nfiles; j++) {
      strcpy(newfile,argv[j+1]);
      if(!(grid[j] = read_lcdisp(newfile,&tauchk,&muchk))) 
	no_error = 0;
      if(tauchk.nval != tau.nval || muchk.nval != mu.nval) {
	fprintf(stderr,"ERROR: Sizes of grids do not match (%s).\n",
		newfile);
	fprintf(stderr,"       ntau = %d and %d, nmu = %d and %d.\n",
		tau.nval,tauchk.nval,mu.nval,muchk.nval);
	no_error = 0;
      }
    }
  }

  /*
   * Allocate memory for mindisp and finaldisp arrays
   */

  if(no_error) {
    for(j=0; j<nfiles; j++)
      if(!(mindisp[j] = new_lcdisp(tau.nval)))
	no_error = 0;
    if(!(finaldisp = new_lcdisp(tau.nval)))
      no_error = 0;
  }

  /*
   * Loop through all grid values associated with each value of tau
   *  and find the one with the minimum value of disp.  Store the positions
   *  of these points in the mindisp arrays.
   */

  if(no_error) {
    for(j=0; j<nfiles; j++) {
      for(i=0,dptr=mindisp[j]; i<tau.nval; i++,dptr++) {

	/*
	 * Load the tau values into the mindisp array for this
	 *  position (set by i).
	 */

	*dptr = *(grid[j]+i);

	/*
	 * Initialize the disp, mu, and arraypos values.
	 */

	dptr->disp = 9999.9;
	dptr->mu = 0.0;
	dptr->arraypos = 0;

	/*
	 * Loop through grid with a fixed value of tau (set by the i loop).
	 *  Look for the minimum value of disp at this value of tau, and
	 *  record the value, the mu value at the grid point, and the
	 *  grid position in the disp, mu, and arraypos variables.
	 */

	for(k=0; k<mu.nval; k++) {
	  gptr = grid[j]+i+(k*tau.nval);
	  dptr->tau = gptr->tau;
	  if(gptr->disp < dptr->disp) {
	    dptr->disp = gptr->disp;
	    dptr->mu = gptr->mu;
	    dptr->arraypos = k;
	  }
	}
      }
    }
  }

  /*
   * Combine the mindisp information to find the lowest total dispersion.
   */

  if(no_error) {

    /*
     * Initialize
     */
    dispmin = 9999.9;
    mincheck.arraypos = 0;

    /*
     * Loop through, checking combined disp for min value.
     */

    for(i=0; i<tau.nval; i++) {
      mincheck.disp = 0.0;
      for(j=0; j<nfiles; j++)
	mincheck.disp += (mindisp[j]+i)->disp;
      if(mincheck.disp < dispmin) {
	dispmin = mincheck.disp;
	mincheck.arraypos = i;
      }
    }

    /*
     * Set mincheck to best value
     */

    i = mincheck.arraypos;
    mincheck.disp = 0.0;
    mincheck.tau = (mindisp[0]+i)->tau;
    for(j=0; j<nfiles; j++)
      mincheck.disp += (mindisp[j]+i)->disp;

    /*
     * Print out results
     */

    if(!(ofp = open_writefile("minpos.dat")))
      no_error = 0;
    else {
      minpos = mincheck.arraypos;
      printf("\n");
      printf("Minimum dispersion=%7.4f at:\n",mincheck.disp);
      fprintf(ofp,"%6.2f ",mincheck.tau);
      for(j=0; j<nfiles; j++) {
	printf("  mu_%d = %6.4f\n",j+1,(mindisp[j]+minpos)->mu);
	fprintf(ofp,"%6.4f ",(mindisp[j]+minpos)->mu);
      }
      fprintf(ofp,"%7.4f\n",dispmin);
      printf("  tau = %6.2f\n",mincheck.tau);
    }
  }

  /*
   * Calculate shifts in mu between the data sets
   */

  if(no_error) {
    if(!(mu_shift = new_intarray(nfiles,1)))
      no_error = 0;
    else {
      minpos = (mindisp[0] + mincheck.arraypos)->arraypos;
      for(j=0; j<nfiles; j++) {
	mu_shift[j] = (mindisp[j] + mincheck.arraypos)->arraypos - minpos;
	if(fabs(mu_shift[j]) > maxshift)
	  maxshift = fabs(mu_shift[j]);
      }
      printf("\nMaximum grid shift in mu is %d.\n",maxshift);
    }
  }

  /*
   * Add dispersion grids with appropriate shifts to get final (tau,disp)
   *  curve and print them out with a call to print_disp_slice.
   */

  if(no_error) {
    for(j=0; j<nfiles; j++) {
      gptr = grid[j] + 
	((mindisp[j] + mincheck.arraypos)->arraypos)*tau.nval;
      for(i=0,dptr=finaldisp; i<tau.nval; i++,dptr++,gptr++) {
	dptr->tau = gptr->tau;
	dptr->mu = gptr->mu;
	dptr->disp += gptr->disp;
      }
    }
    if(print_disp_slice(finaldisp,tau.nval,"disp_comb.dat_slice"))
      no_error = 0;
  }

  /*
   * Clean up
   */

  if(no_error)
    printf("\nCleaning up\n");

  if(ofp)
    fclose(ofp);
  mu_shift = del_intarray(mu_shift);
  finaldisp = del_lcdisp(finaldisp);
  for(i=0; i<nfiles; i++) {
    grid[i] = del_lcdisp(grid[i]);
    mindisp[i] = del_lcdisp(mindisp[i]);
  }
  if(grid)
    free(grid);
  if(mindisp)
    free(mindisp);

  /*
   * Exit
   */

  if(no_error) {
    printf("\nFinished with program add_disp.c\n\n");
    return 0;
  }
  else {
    fprintf(stderr,"ERROR: Exiting program add_disp.c\n\n");
    return 1;
  }
}
