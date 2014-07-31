/*
 * fitquad.c
 *
 * Usage: fitquad [infile] [outfile]
 *
 * This program fits a quadratic (i.e., y = a1 + a2 x + a3 x^2) to a data set
 *  consisting of (x,y) pairs that have errors on the y values.
 * The input file should have 3 columns: x y err_y
 * The output file will have 4 columns:  x y err_y y_mod
 *
 * 08Jun2006 CDF
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
 * Function declarations
 *
 */

double determinant3x3(double a11, double a12, double a13, double a21,
		      double a22, double a23, double a31, double a32,
		      double a33);
double determinant2x2(double a11, double a12, double a21, double a22);
double *invert3x3(double a11, double a12, double a13, double a21,
		  double a22, double a23, double a31, double a32,
		  double a33);


/*.......................................................................
 *
 * Main program
 *
 */

int main(int argc, char *argv[])
{
  int i,j;                   /* Looping variables */
  int no_error=1;            /* Flag set to 0 on error */
  int logx=0;                /* Flag set to 1 if x --> logx before fitting */
  int logy=0;                /* Flag set to 1 if y --> logy before fitting */
  int ndata;                 /* Number of data points */
  int ndof;                  /* Number of degrees of freedom */
  double x,y,sig;            /* Temporary values of input data */
  double sig2;               /* Square of current value of err_y */
  double S,Sx,Sxx;           /* Temporary values for the fitting */
  double Sx3,Sx4;            /* Temporary values for the fitting */
  double alph11,alph12;      /* Temporary values for the fitting */
  double alph13,alph22;      /* Temporary values for the fitting */
  double alph23,alph33;      /* Temporary values for the fitting */
  double a[3];               /* Fitted parameters a1 and a2 */
  double beta[3];            /* Beta vector  */
  double *covar=NULL;        /* Covariance matrix elements */
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
    fprintf(stderr,"\nUsage: fitquad [infile] [outfile]\n\n");
    fprintf(stderr,"This program fits a quadratic ");
    fprintf(stderr,"(i.e., y = a1 + a2 x + a3 x^2) to a\n");
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
  Sx3 = 0.0;
  Sx4 = 0.0;
  beta[0] = 0.0;
  beta[1] = 0.0;
  beta[2] = 0.0;

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
      Sx3 += x * Sxx;
      Sx4 += x * Sx3;
      beta[0] += y / sig2;
      beta[1] += x * y / sig2;
      beta[2] += x * x * y / sig2;
    }
    for(i=0; i<3; i++) 
      printf("beta[%d] = %g\n",i,beta[i]);
    printf("\n");
  }

  /*
   * Put into alpha array
   */

  alph11 = S;
  alph12 = Sx;
  alph13 = Sxx;
  alph22 = Sxx;
  alph23 = Sx3;
  alph33 = Sx4;
  printf("%g %g %g\n",alph11,alph12,alph13);
  printf("%g %g %g\n",alph12,alph22,alph23);
  printf("%g %g %g\n",alph13,alph23,alph33);

  /*
   * Invert alpha matrix (which is symmetric).  The result will
   *  be the covariance matrix
   */

  if(!(covar = invert3x3(alph11,alph12,alph13,alph12,alph22,alph23,
			 alph13,alph23,alph33)))
    no_error = 0;

  /*
   * Use the covariance matrix to find the values of the parameters
   */

  if(no_error) {
    for(i=0; i<3; i++) {
      a[i] = 0.0;
      printf("For a[%d] multiplying: ",i+1);
      for(j=0; j<3; j++) {
	a[i] += *(covar + 3*i + j) * beta[j];
	printf("covar[%d] * beta[%d], ",3*i+j,j);
      }
      printf("\n  a[%d] = %g\n",i+1,a[i]);
      printf("\n");
    }
  }

  /*
   * Give results of fitting
   */

  if(no_error) {
    printf("\nResults of fitting, expressed as y = a1 + a2 * x + a3 * x^2\n");
    for(i=0; i<3; i++)
      printf("  a%d = %g\n",i+1,a[i]);
    printf("\nCovariance matrix\n");
    printf("   %7g    %7g     %7g\n",covar[0],covar[1],covar[2]);
    printf("   %7g    %7g     %7g\n",covar[3],covar[4],covar[5]);
    printf("   %7g    %7g     %7g\n",covar[6],covar[7],covar[8]);
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
    ndof = ndata - 3;
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
      ymod = a[0] + a[1] * x + a[2] * x * x;
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
  covar = del_doubarray(covar);
  if(ofp)
    fclose(ofp);

  if(no_error) {
    printf("\nProgram fitquad finished.\n\n");
    return 0;
  }
  else {
    fprintf(stderr,"\nERROR.  Exiting fitquad.\n\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function invert3x3
 *
 * Inverts a 3x3 matrix
 *
 * Inputs: double a11, etc.    Input matrix elements
 *
 * Output: double *invmatx     Inverse matrix, where the nine elements
 *                              correspond to:
 *                                m11 = invmatx[0]
 *                                m12 = invmatx[1]
 *                                m13 = invmatx[2]
 *                                m21 = invmatx[3]
 *                                etc.
 */

double *invert3x3(double a11, double a12, double a13, double a21,
		  double a22, double a23, double a31, double a32,
		  double a33)
{
  double deta;             /* Determinant of input matrix */
  double *invmatx=NULL;    /* Inverse matrix */

  /*
   * First calculate matrix determinant
   */

  deta = determinant3x3(a11,a12,a13,a21,a22,a23,a31,a32,a33);
  if(deta == 0.0) {
    fprintf(stderr,"ERROR: invert3x3.  Input matrix is singular.\n");
    return NULL;
  }
  printf("invert3x3: determinant = %g\n",deta);

  /*
   * Allocate memory for output matrix
   */

  if(!(invmatx = new_doubarray(9))) {
    fprintf(stderr,"ERROR: invert3x3\n");
    return NULL;
  }

  /*
   * Do ugly inversion
   */

  invmatx[0] = determinant2x2(a22,a23,a32,a33) / deta;
  invmatx[1] = determinant2x2(a13,a12,a33,a32) / deta;
  invmatx[2] = determinant2x2(a12,a13,a22,a23) / deta;
  invmatx[3] = determinant2x2(a23,a21,a33,a31) / deta;
  invmatx[4] = determinant2x2(a11,a13,a31,a33) / deta;
  invmatx[5] = determinant2x2(a13,a11,a23,a21) / deta;
  invmatx[6] = determinant2x2(a21,a22,a31,a32) / deta;
  invmatx[7] = determinant2x2(a12,a11,a32,a31) / deta;
  invmatx[8] = determinant2x2(a11,a12,a21,a22) / deta;

  /*
   * Return
   */

  return invmatx;
}

/*.......................................................................
 *
 * Function determinant2x2
 *
 * Determines the determinant of a 2x2 matrix
 *
 * Inputs: double a11          The elements of the matrix
 *         double a12
 *         double a21
 *         double a22
 *
 * Output: double determ       Determinant of the matrix
 *
 */

double determinant2x2(double a11, double a12, double a21, double a22)
{
  double determ=0.0;   /* Determinant */

  determ = a11 * a22 - a12 * a21;

  return determ;
}

/*.......................................................................
 *
 * Function determinant3x3
 *
 * Calculates the determinant of a 3x3 matrix
 *
 * Inputs: double a11, etc.    Members of the input matrix
 *
 * Output: double determ       Determinant of the matrix
 *
 */

double determinant3x3(double a11, double a12, double a13, double a21,
		      double a22, double a23, double a31, double a32,
		      double a33)
{
  double determ=0.0;    /* Determinant of the matrix */

  determ =  a11*a22*a33 - a11*a23*a32;
  determ += a12*a23*a31 - a12*a21*a33;
  determ += a13*a21*a32 - a13*a22*a31;

  return determ;
}

