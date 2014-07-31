/*
 * fitline.c
 *
 * Usage: fitline [infile] [outfile]
 *
 * This program fits a line (i.e., y = a1 + a2 x) to a data set
 *  consisting of (x,y) pairs that have errors on the y values.
 * The input file should have 3 columns: x y err_y
 * The output file will have 4 columns:  x y err_y y_mod
 *
 * 07Jun2006 CDF
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
 * Main program
 *
 */

int main(int argc, char *argv[])
{
  int i;                     /* Looping variable */
  int no_error=1;            /* Flag set to 0 on error */
  int logx=0;                /* Flag set to 1 if x --> logx before fitting */
  int logy=0;                /* Flag set to 1 if y --> logy before fitting */
  int ndata;                 /* Number of data points */
  int ndof;                  /* Number of degrees of freedom */
  double x,y,sig;            /* Temporary values of input data */
  double sig2;               /* Square of current value of err_y */
  double S,Sx,Sxx;           /* Temporary values for the fitting */
  double Sy,Sxy, Delta;      /* Temporary values for the fitting */
  double a[2];               /* Fitted parameters a1 and a2 */
  double sig11,sig22,sig12;  /* Covariance matrix elements */
  double ymod;               /* Model y values */
  double ydiff;              /* Difference between model and observed y */
  double chisq;              /* Chi^sq */
  char line[MAXC];           /* General string for reading input */
  char infile[MAXC];         /* Filename for input position file */
  char outfile[MAXC];        /* Filename for distcalc-like output file */
  Datastruct *data=NULL;     /* Observed (x,y,err_y) */
  Datastruct *dptr;          /* Pointer to navigate data */
  FILE *ofp=NULL;            /* Output file pointer */

  /*
   * Check the command line invocation
   */

  if(argc < 3) {
    fprintf(stderr,"\nUsage: fitline [infile] [outfile]\n\n");
    fprintf(stderr,"This program fits a line (i.e., y = a1 + a2 x) to a\n");
    fprintf(stderr," data set consisting of (x,y) pairs that have errors \n");
    fprintf(stderr," on the y values.\n");
    fprintf(stderr," The input file should have 3 columns: x y err_y\n");
    fprintf(stderr," The output file will have 4 columns:  x y err_y y_mod\n");
    fprintf(stderr,"\n");
    return 1;
  }
  printf("\n");

  /*
   * Get the inputs from the command line.
   */

  strcpy(infile,argv[1]);
  strcpy(outfile,argv[2]);

  /*
   * Read the input file
   */

  if(!(data = read_datastruct(infile,'#',&ndata,3)))
      no_error = 0;

  /*
   * Get information about logs
   */

  printf("\nShould x be converted to log_10(x) before fitting? [n] ");
  fgets(line,MAXC,stdin);
  if(line[0]=='y' || line[0]=='Y') {
    printf("Converting x to log_10(x)\n");
    logx = 1;
  }

  printf("\nShould y be converted to log_10(y) before fitting? [n] ");
  fgets(line,MAXC,stdin);
  if(line[0]=='y' || line[0]=='Y') {
    printf("Converting y to log_10(y)\n");
    logy = 1;
  }

  /*
   * Initialize
   */

  S = 0.0;
  Sx = 0.0;
  Sxx = 0.0;
  Sy = 0.0;
  Sxy = 0.0;

  /*
   * Loop through the input file to generate the temporary variables
   */

  if(no_error) {
    for(i=0,dptr=data; i<ndata; i++,dptr++) {
      if(logx) {
	x = log10(dptr->x);
      } else {
	x = dptr->x;
      }
      if(logy) {
	y = log10(dptr->y);
	sig = dptr->z / (log(10.0) * dptr->y);
      } else {
	y = dptr->y;
	sig = dptr->z;
      }
      sig2 = sig * sig;
      S += 1.0 / sig2;
      Sx += x / sig2;
      Sxx += x * x / sig2;
      Sy += y / sig2;
      Sxy += x * y / sig2;
    }
    Delta = S * Sxx - Sx * Sx;
  }

  /*
   * Convert temporary variables into fitted parameters and
   *  covariance matrix
   */

  if(no_error) {
    a[0] = (Sy * Sxx - Sx * Sxy)/Delta;
    a[1] = (Sxy * S - Sx * Sy)/Delta;
    sig11 = Sxx / Delta;  
    sig22 = S / Delta;    
    sig12 = -Sx / Delta;  
  }


  /*
   * Give results of fitting
   */

  if(no_error) {
    printf("\nResults of fitting, expressed as y = a1 + a2 * x\n");
    printf("  a1 = %f\n",a[0]);
    printf("  a2 = %f\n",a[1]);
    printf("\nCovariance matrix\n");
    printf("   %7g    %7g\n",sig11,sig12);
    printf("   %7g    %7g\n",sig12,sig22);
  }

  /*
   * Open output file
   */

  if(no_error)
    if(!(ofp = open_writefile(outfile)))
      no_error = 0;

  /*
   * Write to output file and, at the same time, calculate chi^2
   */

  if(no_error) {
    chisq = 0.0;
    ndof = ndata - 2;
    for(i=0,dptr=data; i<ndata; i++,dptr++) {
      if(logx) {
	x = log10(dptr->x);
      } else {
	x = dptr->x;
      }
      if(logy) {
	y = log10(dptr->y);
	sig = dptr->z / (log(10.0) * dptr->y);
      } else {
	y = dptr->y;
	sig = dptr->z;
      }
      ymod = a[0] + a[1] * x;
      ydiff = y - ymod;
      chisq += ydiff * ydiff / (sig * sig);
      fprintf(ofp,"%f %f %f %f\n",x,y,sig,ymod);
    }
    printf("\nFinal chi^2 = %f\n",chisq);
    printf("Degrees of freedom = %d\n",ndof);
    printf("Reduced chi^2 = %f\n",chisq/ndof);
  }

  /*
   * Clean up and exit
   */

  data = del_datastruct(data);
  /* ymod = del_doubarray(ymod); */
  if(ofp)
    fclose(ofp);

  if(no_error) {
    printf("\nProgram fitline finished.\n\n");
    return 0;
  }
  else {
    fprintf(stderr,"\nERROR.  Exiting fitline.\n\n");
    return 1;
  }
}

