/*
 * xy2spos.c
 *
 * Converts an (x,y) CCD position to an (RA,Dec) position using the outputs
 *  of the iraf function "geomap", which has been run with fitgeom=rxyscale.
 *
 * The conversion is:
 *
 *   alpha = a + b*x + c*y
 *   delta = d + e*x + f*y
 *
 * where
 *
 *  a = xshift
 *  b = xmag * cos(xrotation)
 *  c = ymag * sin(yrotation)
 *  d = yshift
 *  e = -xmag * sin(xrotation)
 *  f = ymag * cos(yrotation)
 *
 * and xshift, yshift, xmag, ymag, xrotation, and yrotation are the outputs
 * of geomap.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "structdef.h"
#include "dataio.h"
#include "coords.h"

/*.......................................................................
 *
 * Function declarations
 *
 */

int find_conv_const(double *a, double *b, double *c, double *d, double *e,
		    double *f);

/*.......................................................................
 *
 * Main program
 *
 */


int main(int argc, char *argv[])
{
  int no_error=1;         /* Flag set to 0 on error */
  int nlines=0;           /* Number of lines in the input file */
  double x,y;             /* Input (x,y) positions from input file */
  double a,b,c,d,e,f;     /* Conversion constants from x,y to alpha,delta */
  double alpha,delta;     /* Sky positions in decimal degrees */
  char line[MAXC];        /* General string variable for getting input */
  Skypos spos;            /* Sky position in hh:mm:ss dd:mm:ss format */
  FILE *ifp=NULL;         /* Input file pointer */

  /*
   * Open the input file
   */

  if(!(ifp = open_readfile(argv[1]))) {
    fprintf(stderr,"ERROR: calc_cent.\n");
    return 1;
  }

  /*
   * Count the number of lines in the input file
   */

  if(no_error) {
    if((nlines = n_lines(ifp,'#')) == 0) {
      fprintf(stderr,"ERROR.  %s has no valid data lines.\n",argv[1]);
      no_error = 0;
    }
    else {
      rewind(ifp);
      printf("\nRead %d lines from %s\n",nlines,argv[1]);
    }
  }

  /*
   * Calculate the conversion constants
   */

  if(find_conv_const(&a,&b,&c,&d,&e,&f))
    no_error = 0;
  else
    printf("%f %f %f %f %f %f\n",a,b,c,d,e,f);

  /*
   * Step through the input file and convert each line
   */

  while(fgets(line,MAXC,ifp) != NULL && no_error) {
    if(line[0] != '#') {
      if(sscanf(line,"%lf %lf",&x,&y) != 2) {
	fprintf(stderr,"ERROR: bad format in input file.\n");
	fprintf(stderr," in line: %s",line);
	no_error = 0;
      }
      else {
	alpha = a + b*x + c*y;
	delta = d + e*x + f*y;
	deg2spos(alpha,delta,&spos);
	printf("%7.2f %7.2f   %02d %02d %07.4f %+03d %02d %06.3f\n",
	       x,y,spos.hr,spos.min,spos.sec,spos.deg,spos.amin,spos.asec);
      }
    }
  }

  /*
   * Clean up and exit
   */

  if(ifp)
    fclose(ifp);

  return 0;
}

/*.......................................................................
 *
 * Function find_conv_const
 *
 * Converts the values from geomap to the conversion constants a->f.
 *
 * Inputs: double *a           Conversion constants, set by this function
 *         double *b                               |
 *         double *c                               |
 *         double *d                               |
 *         double *e                             \ | /
 *         double *f                              \|/
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int find_conv_const(double *a, double *b, double *c, double *d, double *e,
		    double *f)
{

  double xshift,yshift;   /* Translations output from geomap */
  double xmag,ymag;       /* Magnifications output from geomap */
  double xrot,yrot;       /* Rotations output from geomap */
  char line[MAXC];        /* General string variable for getting input */

  /*
   * Initialize the variables
   */

  xshift = yshift = 0.0;
  xmag = ymag = 1.0;
  xrot = yrot = 0.0;

  /*
   * Get the variables from the user
   */

  printf("\nEnter xshift and yshift, separated by spaces:  ");
  fgets(line,MAXC,stdin);
  while(sscanf(line,"%lf %lf",&xshift,&yshift) != 2) {
    fprintf(stderr,"ERROR: bad input.  Enter values again:  ");
    fgets(line,MAXC,stdin);
  }

  printf("\nEnter xmag and ymag, separated by spaces:  ");
  fgets(line,MAXC,stdin);
  while(sscanf(line,"%lf %lf",&xmag,&ymag) != 2) {
    fprintf(stderr,"ERROR: bad input.  Enter values again:  ");
    fgets(line,MAXC,stdin);
  }

  printf("\nEnter xrotation and yrotation, separated by spaces:  ");
  fgets(line,MAXC,stdin);
  while(sscanf(line,"%lf %lf",&xrot,&yrot) != 2) {
    fprintf(stderr,"ERROR: bad input.  Enter values again:  ");
    fgets(line,MAXC,stdin);
  }

  /*
   * Do the calculations
   */

  *a = xshift;
  *b = xmag * cos(xrot * PI / 180.0);
  *c = ymag * sin(yrot * PI / 180.0);
  *d = yshift;
  *e = -xmag * sin(xrot * PI / 180.0);
  *f =  ymag * cos(yrot * PI / 180.0);

  return 0;
}
