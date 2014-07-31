/*
 * randgen.c
 *
 * Usage: randgen [position_filename] [input_catalog] [output_filename] 
 *   [calcmethod] ([format])
 *
 * This program generates a list of uniformly distributed random numbers.
 * The program uses the Numerical Recipes ran2 function (included here as
 *  myran2) to create the uniform random sample.
 * The output can then be used in SM to transform to other distributions.
 *
 * 09May2006 CDF, First working version
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "structdef.h"
#include "dataio.h"


/*.......................................................................
 *
 * Definitions for myran2
 *
 */

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

/*.......................................................................
 *
 * Function declaration
 *
 */

float myran2(long *idum);

/*.......................................................................
 *
 * Main program
 *
 */

int main(int argc, char *argv[])
{
  int i;                    /* Looping variable */
  int no_error=1;           /* Flag set to 0 on error */
  int nrand=0;              /* Number of random numbers to be generated. */
  long randseed;            /* Seed for random number generator. */
  float tmprand;            /* Temporary holder for the random number */
  char outfile[MAXC];       /* Output file name */
  FILE *ofp=NULL;           /* Pointer for output file */

  
  /*
   * Check the command line invocation
   */

  if(argc < 3) {
    fprintf(stderr,"\nUsage: randgen [outfile] [n]\n\n");
    fprintf(stderr,"  outfile is the requested output file name.\n");
    fprintf(stderr,"  n is the number of random numbers to be generated\n\n");
    return 1;
  }
  printf("\n");

  /*
   * Get the inputs from the command line.
   */

  strcpy(outfile,argv[1]);
  if(sscanf(argv[2],"%d",&nrand) != 1 || nrand < 1) {
    fprintf(stderr,"ERROR: n must be a positive integer.\n\n");
    return 1;
  }
  printf("\nWill generate %d random numbers.\n",nrand);

  /*
   * Set the initial seed for the random number generator
   */

  randseed = -97867564;

  /*
   * Open the output file for writing
   */

  if(!(ofp = open_writefile(outfile))) {
    fprintf(stderr,"ERROR: Cannot open output file.\n\n");
    return 1;
  }

  /*
   * Fill the day members of the four light curves with uniform
   *  random deviates between 0 and 1.
   */

  for(i=0; i<nrand; i++) {
    tmprand = myran2(&randseed);
    fprintf(ofp,"%f\n",tmprand);
  }

  /*
   * Close the output file
   */

  if(ofp)
    fclose(ofp);

  /*
   * Clean up and exit;
   */

  if(no_error) {
    printf("\nProgram randgen finished.\n\n");
    return 0;
  }
  else {
    fprintf(stderr,"\nERROR.  Exiting randgen.\n\n");
    return 1;
  }
}


/*.......................................................................
 *
 * Function myran2
 *
 * A transcription of the Numerical Recipes ran2 function.
 *
 */

float myran2(long *idum)
{
  int j;
  long k;
  static long idum2=123456789;
  static long iy=0;
  static long iv[NTAB];
  float temp;

  /*
   * Initialize
   */

  if (*idum <= 0) {
    if (-(*idum) < 1)
      *idum = 1;
    else
      idum2 = (*idum);

    /*
     * Load the shuffle table after 8 warm-ups.
     */

    for (j=NTAB+7; j>=0; j--) {
      k = (*idum)/IQ1;
      *idum = IA1 * (*idum - k * IQ1) - k*IR1;
      if (*idum < 0)
	*idum += IM1;
      if (j<NTAB)
	iv[j] = *idum;
    }
    iy = iv[0];
  }

  /*
   * Start here when not initializing
   */

  k = (*idum) / IQ1;

  /*
   * Compute idum = (IA1*idum)%IM1 without overflows by Schrage's method
   */

  *idum= IA1 * (*idum - k * IQ1) - k*IR1;
  if (*idum<0)
    *idum += IM1;
  k = idum2/IQ2;

  /*
   * Calculate idum2 in the same way
   */

  idum2 = IA2*(idum2-k*IQ2)-k*IR2;
  if(idum2<0)
    idum2 += IM2;

  /*
   * Generate j in the range 0 to NTAB-1
   */

  j = iy/NDIV;

  /*
   * Shuffle idum and then combine idum and idum2 to give output.
   */

  iy = iv[j] - idum2;
  iv[j] = *idum;
  if(iy<1)
    iy+= IMM1;

  if((temp=AM*iy) > RNMX) 
    return RNMX;
  else
    return temp;
}
