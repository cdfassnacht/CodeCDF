/*
 * add_chisq.c
 *
 * Usage: 
 *
 * Takes as input two sets of identically sampled gridpoints in 
 *  (delay,magnification,chisq) space, as produced by running delays.c on
 *  two separate sets of light curves.  The program will then find the
 *  minimum COMBINED chisq, allowing for shifts in magnification between
 *  seasons, to deal with possible changes in magnifications caused by 
 *  microlensing.  The input data sets have 4 parameter values: delay, 
 *  magnification, chisq, and N_DOF, and will be contained in the chib*.dat 
 *  files produced by choosing the chisq-minimization option in delays.c.  
 *
 * 11Sep00 CDF,  First working version.
 * v10Oct00 CDF, Fixed bug in error checking.
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
  int minpos;               /* Array position associated with rchimin */
  int *mu_shift=NULL;       /* Grid shifts between datasets */
  int maxshift=0;           /* Maximum value of mu_shift */
  float rchimin;            /* Minimum combined reduced chisq */
  Axisinfo tau,tauchk;      /* Tau axis info (and check) */
  Axisinfo mu,muchk;        /* Mu axis info (and check) */
  char file1[MAXC];         /* Input file name for first file */
  char newfile[MAXC];       /* Input file name for other input file(s) */
  LCchisq mincheck;         /* Container used in min chisq search */
  LCchisq **grid={NULL};    /* Grid values */
  LCchisq **minchi={NULL};  /* Vectors containing minimum chisq info */
  LCchisq *finalchi=NULL;   /* Vectors containing minimum chisq info */
  LCchisq *gptr;            /* Pointer for navigating the grids */
  LCchisq *chiptr;          /* Pointer for navigating LCchisq arrays */
  FILE *ofp=NULL;           /* Output file pointer for location of minimum */
 
  /*
   * Check input line
   */

  if(argc < 3) {
    fprintf(stderr,"\nUsage: add_chisq [file1] [file2] ... [fileN]\n");
    fprintf(stderr," where the input files are chisq files produced by");
    fprintf(stderr," separate runs of delays.c.\n\n");
    fprintf(stderr,"*** NB: The input parameter-space grids must have ");
    fprintf(stderr,"identical sizes and spacings.\n\n");
    return 1;
  }
  else {
    nfiles = argc - 1;
    printf("\nRunning add_chisq.c with %d input files.\n",nfiles);
  }

  /*
   * Allocate first level of grid pointers and minchi pointers
   */

  grid = (LCchisq **) malloc(sizeof(LCchisq *) * nfiles);
  if(!grid) {
    fprintf(stderr,"ERROR:  Insufficient memory for grid array.\n");
    return 1;
  }

  minchi = (LCchisq **) malloc(sizeof(LCchisq *) * nfiles);
  if(!minchi) {
    fprintf(stderr,"ERROR:  Insufficient memory for minchi array.\n");
    no_error = 0;
  }

  /*
   * Read data from first input file
   */

  if(no_error) {
    strcpy(file1,argv[1]);
    if(!(grid[0] = read_lcchisq(file1,&tau,&mu))) {
      fprintf(stderr,"ERROR: Exiting add_chisq.c.\n");
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
      if(!(grid[j] = read_lcchisq(newfile,&tauchk,&muchk))) 
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
   * Allocate memory for minchi and finalchi arrays
   */

  if(no_error) {
    for(j=0; j<nfiles; j++)
      if(!(minchi[j] = new_lcchisq(tau.nval)))
	no_error = 0;
    if(!(finalchi = new_lcchisq(tau.nval)))
      no_error = 0;
  }

  /*
   * Loop through all grid values associated with each value of tau
   *  and find the one with the minimum value of chisq.  Store the positions
   *  of these points in the minchi arrays.
   */

  if(no_error) {
    for(j=0; j<nfiles; j++) {
      for(i=0,chiptr=minchi[j]; i<tau.nval; i++,chiptr++) {

	/*
	 * Load the tau and N_dof values into the minchi array for this
	 *  position (set by i).
	 */

	*chiptr = *(grid[j]+i);

	/*
	 * Initialize the chisq, mu, and arraypos values.
	 */

	chiptr->chisq = 9999.9;
	chiptr->mu = 0.0;
	chiptr->arraypos = 0;

	/*
	 * Loop through grid with a fixed value of tau (set by the i loop).
	 *  Look for the minimum value of chisq at this value of tau, and
	 *  record the value, the mu value at the grid point, and the
	 *  grid position in the chisq, mu, and arraypos variables.
	 */

	for(k=0; k<mu.nval; k++) {
	  gptr = grid[j]+i+(k*tau.nval);
	  chiptr->tau = gptr->tau;
	  chiptr->ndof = gptr->ndof;
	  if(gptr->chisq < chiptr->chisq) {
	    chiptr->chisq = gptr->chisq;
	    chiptr->mu = gptr->mu;
	    chiptr->arraypos = k;
	  }
	}
      }
    }
  }

  /*
   * Combine the minchi information to find the lowest reduced chi^2.
   */

  if(no_error) {

    /*
     * Initialize
     */
    rchimin = 9999.9;
    mincheck.arraypos = 0;

    /*
     * Loop through, checking combined reduced chisq for min value.
     */

    for(i=0; i<tau.nval; i++) {
      mincheck.chisq = 0.0;
      mincheck.ndof = 0;
      for(j=0; j<nfiles; j++) {
	mincheck.chisq += (minchi[j]+i)->chisq;
	mincheck.ndof  += (minchi[j]+i)->ndof;
      }
      if((mincheck.rchisq = mincheck.chisq / mincheck.ndof) < rchimin) {
	rchimin = mincheck.rchisq;
	mincheck.arraypos = i;
      }
    }

    /*
     * Set mincheck to best value
     */

    i = mincheck.arraypos;
    mincheck.chisq = 0.0;
    mincheck.ndof = 0;
    mincheck.tau = (minchi[0]+i)->tau;
    for(j=0; j<nfiles; j++) {
      mincheck.chisq += (minchi[j]+i)->chisq;
      mincheck.ndof  += (minchi[j]+i)->ndof;
    }

    /*
     * Print out results
     */

    if(!(ofp = open_writefile("minpos.dat")))
      no_error = 0;
    else {
      minpos = mincheck.arraypos;
      printf("\n");
      printf("Minimum Chisq=%8.2f, N_dof=%d, ",mincheck.chisq,mincheck.ndof);
      printf("Red_chisq=%7.3f at:\n",rchimin);
      fprintf(ofp,"%6.2f ",mincheck.tau);
      for(j=0; j<nfiles; j++) {
	printf("  mu_%d = %6.4f\n",j+1,(minchi[j]+minpos)->mu);
	fprintf(ofp,"%6.4f ",(minchi[j]+minpos)->mu);
      }
      fprintf(ofp,"%8.2f %d %7.3f\n",mincheck.chisq,mincheck.ndof,rchimin);
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
      minpos = (minchi[0] + mincheck.arraypos)->arraypos;
      for(j=0; j<nfiles; j++) {
	mu_shift[j] = (minchi[j] + mincheck.arraypos)->arraypos - minpos;
	if(fabs(mu_shift[j]) > maxshift)
	  maxshift = fabs(mu_shift[j]);
      }
      printf("\nMaximum grid shift in mu is %d.\n\n",maxshift);
    }
  }

  /*
   * Add chi^2 grids with appropriate shifts to get final (tau,rchisq)
   *  curve and print them out with a call to print_chisq_slice.
   */

  if(no_error) {
    for(j=0; j<nfiles; j++) {
      gptr = grid[j] + 
	((minchi[j] + mincheck.arraypos)->arraypos)*tau.nval;
      for(i=0,chiptr=finalchi; i<tau.nval; i++,chiptr++,gptr++) {
	chiptr->tau = gptr->tau;
	chiptr->mu = gptr->mu;
	chiptr->chisq += gptr->chisq;
	chiptr->ndof += gptr->ndof;
	chiptr->rchisq = chiptr->chisq / chiptr->ndof;
      }
    }
    if(print_chisq_slice(finalchi,tau.nval,"chi_comb.dat_slice"))
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
  finalchi = del_lcchisq(finalchi);
  for(i=0; i<nfiles; i++) {
    grid[i] = del_lcchisq(grid[i]);
    minchi[i] = del_lcchisq(minchi[i]);
  }
  if(grid)
    free(grid);
  if(minchi)
    free(minchi);

  /*
   * Exit
   */

  if(no_error) {
    printf("\nFinished with program add_chisq.c\n\n");
    return 0;
  }
  else {
    fprintf(stderr,"ERROR: Exiting program add_chisq.c\n\n");
    return 1;
  }
}
