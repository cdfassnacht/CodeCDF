/* Library: fitsim_contours.c
 *
 * Contains contouring section of fits plotting package.  The functions
 *  in this library include those which get contouring information and
 *  actually plot the contours.
 *
 * 04Dec99 CDF,  Split out of fitsim.c
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "structdef.h"
#include "libfits.h"
#include "cpgplot.h"
#include "plotfuncs.h"
#include "fitsim.h"

/*.......................................................................
 *
 * Function get_contours
 *
 * Gets number of contours and contour levels, if needed.  If contours
 *  are desired, this function allocates memory for the contour array.
 *
 * Inputs: Setup *setup         image display information
 *
 * Output: int (0 or 1)         0 ==> success, 1 ==> error
 *
 * v06Nov99 CDF, Moved contour color selection from display_image
 */

int get_contours(Setup *setup)
{
  int cchoice=3;    /* Choice of contour definition */
  char line[MAXC];  /* General string for getting input */

  if(setup->docontour == UNSET) {

    /*
     * Initialize
     */

    setup->ncont = 0;

    /*
     * Get the desired number of contours
     */

    printf("\nEnter the number of contours to plot: [%d] ",setup->ncont);
    fgets(line,MAXC,stdin);
    if(line[0] != '\n') {
      while(sscanf(line,"%d",&setup->ncont) != 1 && setup->ncont < 0) {
	fprintf(stderr,"ERROR: bad input.  Enter number of contours:  ");
	fgets(line,MAXC,stdin);
      }
    }

    /*
     * If no contours are desired, just return
     */

    if(setup->ncont == 0)
      return 0;

    /*
     * Otherwise, see how the contours are going to be defined
     */

    setup->docontour = TRUE;
    printf("\nContour plotting choices\n");
    printf("  1. Percentage of peak\n");
    printf("  2. Multiples of a given value\n");
    printf("  3. Multiples of a given value -- log spacing\n");
    printf("---------------------------------\n");
    printf("Enter choice [%d]:  ",cchoice);
    fgets(line,MAXC,stdin);
    if(line[0] != '\n') {
      while(sscanf(line,"%d",&cchoice) != 1 && cchoice < 0 && cchoice > 20) {
	fprintf(stderr,"ERROR: bad input.  Enter choice:  ");
	fgets(line,MAXC,stdin);
      }
    }

    /*
     * Allocate the contour array
     */

    if(!(setup->cont = new_array(setup->ncont,1))) {
      fprintf(stderr,"ERROR: get_contours\n");
      return 1;
    }

    /*
     * Fill in contour array using choice defined above
     */

    switch(cchoice) {
    case 1:
      if(pcontour(setup->cont,setup->ncont,setup->hi)) {
	fprintf(stderr,"ERROR: get_contours\n");
	return 1;
      }
      break;
    case 2:
      if(mcontour(setup->cont,setup->ncont,setup->hi)) {
	fprintf(stderr,"ERROR: get_contours\n");
	return 1;
      }
      break;
    case 3:
      if(lcontour(setup->cont,setup->ncont,setup->hi)) {
	fprintf(stderr,"ERROR: get_contours\n");
	return 1;
      }
      break;
    default:
      fprintf(stderr,"ERROR: Not a valid option.  ");
      fprintf(stderr,"No contours will be plotted\n");
      setup->ncont = 0;
    }
  }

  /*
   * Finally see in which color the contours are going to be plotted.
   */

  if(setup->ncont > 0) {
    setup->contcolor = 1;
    printf("\nWhite (0) or Black (1) contours?:  [%d] ",setup->contcolor);
    fgets(line,MAXC,stdin);
    if(line[0] != '\n') {
      while(sscanf(line,"%d",&setup->contcolor) != 1 && setup->contcolor < 0
	    && setup->contcolor > 1) {
	fprintf(stderr,"ERROR: Bad input.  Enter color choice again:  ");
	fgets(line,MAXC,stdin);
      }
    }
  }

  return 0;
}

/*.......................................................................
 *
 * Function pcontour
 *
 * Fills in an array of contour values by setting the values as percentage
 *  of the peak value in the map.
 *
 * Inputs: float *cont         array to be filled in
 *         int ncont           number of contour levels
 *         float max           max value in image being displayed
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int pcontour(float *cont, int ncont, float max)
{
  int i;            /* Looping variable */
  char line[MAXC];  /* General string for getting input */

  /*
   * Get the contour values
   */

  printf("\n");
  for (i=0; i<ncont; i++) {
    printf("Enter contour %d, as a percentage of peak (10 gives 10 percent): ",
	   i+1);
    fgets(line,MAXC,stdin);
    while(sscanf(line,"%f",&cont[i]) != 1) {
      fprintf(stderr,"ERROR: bad input.  Enter contour level:  ");
      fgets(line,MAXC,stdin);
    }

    cont[i] *= max / 100;
  }

  return 0;
}

/*.......................................................................
 *
 * Function mcontour
 *
 * Fills in an array of contour values by setting the values as multiples
 *  of an entered value (e.g. the rms noise in the map).
 *
 * Inputs: float *cont         array to be filled in
 *         int ncont           number of contour levels
 *         float max           max value in image being displayed
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int mcontour(float *cont, int ncont, float max)
{
  int i;            /* Looping variable */
  float cmul;       /* Base contour level */
  char line[MAXC];  /* General string for getting input */

  /*
   * Get the base contour level
   */

  printf("\nEnter the base contour level:  ");
  fgets(line,MAXC,stdin);
  while(sscanf(line,"%f",&cmul) != 1) {
    fprintf(stderr,"ERROR: bad input.  Enter base contour level:  ");
    fgets(line,MAXC,stdin);
  }

  /*
   * Get the contour values
   */

  for (i=0; i<ncont; i++) {
    printf("Enter contour %d, as a factor by which to multiply the ",i+1);
    printf("base level:  ");
    fgets(line,MAXC,stdin);
    while(sscanf(line,"%f",&cont[i]) != 1) {
      fprintf(stderr,"ERROR: bad input.  Enter contour level:  ");
      fgets(line,MAXC,stdin);
    }

    cont[i] *= cmul;
  }

  return 0;
}

/*.......................................................................
 *
 * Function lcontour
 *
 * Fills in an array of contour values by setting the values as multiples
 *  of an entered value, separated by log spacing
 *
 * Inputs: float *cont         array to be filled in
 *         int ncont           number of contour levels
 *         float max           max value in image being displayed
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 * v25Mar99 CDF, Fixed up output format.
 */

int lcontour(float *cont, int ncont, float max)
{
  int i;            /* Looping variable */
  float cmul;       /* Base contour level */
  float first;      /* First contour */
  float *fptr;      /* Pointer to navigate cont */
  char line[MAXC];  /* General string for getting input */

  /*
   * Get the base contour level
   */

  printf("\nEnter the base contour level:  ");
  fgets(line,MAXC,stdin);
  while(sscanf(line,"%f",&cmul) != 1) {
    fprintf(stderr,"ERROR: bad input.  Enter base contour level:  ");
    fgets(line,MAXC,stdin);
  }

  /*
   * Get the first contour
   */

  printf("\nEnter the first value by which to multiply the ");
  printf("base level:  ");
  fgets(line,MAXC,stdin);
  while(sscanf(line,"%f",&first) != 1) {
    fprintf(stderr,"ERROR: Bad input.  Enter first value again:  ");
    fgets(line,MAXC,stdin);
  }

  /*
   * Set the first two contours to be -(first) and first
   */

  fptr = cont;
  *fptr = -1.0 * first * cmul;
  *(fptr+1) = first * cmul;
  printf("  Level  1 = %11.5f (%9.3f * cmul)\n",-first*cmul,-first);
  printf("  Level  2 = %11.5f (%9.3f * cmul)\n",first*cmul,first);

  /*
   * Set the rest of the contour values
   */

  for (i=2,fptr=cont+2; i<ncont; i++,fptr++) {
    *fptr = 2.0 * *(fptr-1);
    printf("  Level %2d = %11.5f (%9.3f * cmul)\n",i+1,*fptr,*fptr/cmul);
  }

  return 0;
}

