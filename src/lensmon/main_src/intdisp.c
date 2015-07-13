/*
 * intdisp.c
 *
 * Usage: intdisp [1608 data file]
 *   where the 1608 file contains lightcurves and errors for the four
 *   components, in the format:
 *
 *      day S_A S_B S_C S_D err_A err_B err_C err_D
 *
 *
 * 19Apr01 CDF,  A modification of delays.c
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

/*.......................................................................
 *
 * Main program
 *
 */

int main(int argc, char *argv[])
{
  int i,j;                      /* Looping variables */
  int no_error=1;               /* Flag set to 0 on error */
  int nflux=0;                  /* Number of points in light curves */
  float disp[N08];              /* Internal dispersions */
  char infile[MAXC];            /* Input file name */
  Fluxrec *flux[N08]={NULL};    /* Flux 1608 light curves */
  Fluxrec *fptr;                /* Pointer to navigate flux arrays */
  FILE *ifp=NULL;               /* File pointer for interpolated data file */
 
  /*
   * Check input line
   */

  if(argc < 2) {
    fprintf(stderr,"\nUsage: intdisp [1608_data_file]\n\n");
    return 1;
  }

  /*
   * Open lens light curve file 
   */

  strcpy(infile,argv[1]);
  if(!(ifp = open_readfile(infile)))
    no_error = 0;

  /*
   * Count the number of data input lines in order to set
   *  sizes of arrays.
   */

  if(no_error) {
    nflux = n_lines(ifp,'#');
  }

  if(nflux == 0) {
    fprintf(stderr,"ERROR. No valid data lines in %s.\n",infile);
    no_error = 0;
  }
  else {
    rewind(ifp);
  }

  /*
   * Allocate arrays
   */

  if(no_error) {
    for(i=0; i<N08; i++) {
      if(!(flux[i] = new_fluxrec(nflux)))
	no_error = 0;
    }
  }

  /*
   * Read in data
   */

  if(no_error) {
    if(read_1608(flux,nflux,ifp,infile))
      no_error = 0;
  }

  /*
   * Set all the "match" values in the component light curves to their
   *  array index value.  This will make all of the match values different
   *  and will allow disp_d1 to work on a single light curve.
   */

  if(no_error) {
    printf("Re-setting match variables.\n");
    for(i=0,fptr=flux[0]; i<nflux; i++,fptr++) {
      (flux[0]+i)->match = i;
      (flux[1]+i)->match = i;
      (flux[2]+i)->match = i;
      (flux[3]+i)->match = i;
    }
  }

  /*
   * Calculate the internal dispersion of all four light curves
   *  through calls to disp_d1 or disp_d2.
   */

  if(no_error) {
    printf("\nCalculating dispersions.\n");
    for(i=0; i<N08; i++) {
      disp[i] = disp_d2(flux[i],nflux,10.0);
    }
  }

  /*
   * Print out results
   */

  printf("\nInternal Dispersion Results\n");
  printf("-----------------------------------\n");
  printf("%7.5f %7.5f %7.5f %7.5f\n\n",disp[0],disp[1],disp[2],disp[3]);

  /*
   * Clean up
   */

  if(no_error)
    printf("\nCleaning up\n\n");

  for(i=0; i<4; i++) {
    flux[i] = del_fluxrec(flux[i]);
  }

  if(ifp)
    fclose(ifp);

  if(no_error)
    return 0;
  else {
    fprintf(stderr,"ERROR: Exiting program intdisp.c\n");
    return 1;
  }
}
