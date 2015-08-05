/* fitting.c
 *
 * A library of multi-dimensional fitting routines.
 *
 */

#include <stdio.h>
#include <math.h>
#include "structdef.h"
#include "list_tools.h"
#include "fitting.h"

#define NMAX 10000

/*.......................................................................
 *
 * Function simplex_setup
 *
 * Creates the arrays necessary for doing a downhill simplex fit.  Note that
 *  the initial nvar members of the vars array should have been set before
 *  calling this function, ie. those values define the initial point in
 *  parameter space that serves as one vertex of the simplex.  The function
 *  will define the other vertices of the simplex by using the dvars
 *  displacement array.  The value of chisq associated with each vertex of
 *  the simplex is calculated using the passed function chifunc, and stored
 *  in the chisqlist array, which is returned.
 *
 * Inputs: float *vars         array containing initial point in parameter
 *                              space (nvar coordinates) plus spaces for
 *                              the additional simplex points (modified
 *                              by this function)
 *         int nvar            number of varying points
 *         float *dvars        array of parameter displacements
 *         float (*chifunc)    function used to calculate chisq
 *         Modinfo *modinfo    lens model info
 *         int verbose         flag set to 1 for verbose output
 *
 * Output: float *chisqlist    array of chisq values associated with simplex
 */

float *simplex_setup(float *vars, int nvar, float *dvars, 
		     float (*chifunc)(float *, void *), void *info,
		     int verbose)
{
  int i,j;                /* Looping variables */
  float *fptr;            /* Pointer used to navigate vars */
  float *fptr2;           /* Pointer used to navigate dvars */
  float *chisqlist=NULL;  /* Chisq values associated with simplex */

  if(verbose)
    printf("\nSetting up the simplex.\n\n");

  /*
   * Go through the vars array, assigning values to each point based on
   *  the initial point in parameter space and the displacement array.
   *  Set fptr to vars+nvar to start with the second point of the simplex.
   */

  fptr = vars + nvar;

  for(i=1; i<=nvar; i++) {
    fptr2 = dvars;
    for(j=0; j<nvar; j++,fptr2++) {
      if(j == i-1)
	*fptr++ = *(vars+j) + *fptr2;
      else
	*fptr++ = *(vars+j);
    }
  }

  /*
   * Now that the simplex is set up, calculate the value of chisq that
   *  is associated with each vertex of the simplex.
   *
   * First allocate the chisq array
   */

  if(!(chisqlist = new_array(nvar+1,1))) {
    fprintf(stderr,"Exiting program\n");
    return NULL;
  }

  /*
   * Now do the calculation, using the function chifunc
   */

  fptr = chisqlist;
  for(i=0; i<=nvar; i++,fptr++) {
    fptr2 = vars + i*nvar;
    if((*fptr = (*chifunc)(fptr2,info)) < 0.0) {
      fprintf(stderr,"ERROR: simplex_setup.  Bad value of chisq\n");
      return NULL;
    }
    if(verbose)
      printf("For simplex point %d chisq = %f\n",i,*fptr);
  }

  return chisqlist;
}

/*.......................................................................
 *
 * Function amoeba
 *
 * A slight modification of the downhill simplex minimization 
 *  routine from Numerical Recipes.
 *
 */

int amoeba(float *p, float y[], int ndim, float ftol,
	   float (*funk)(float [], void *), int *nfunk, void *info,
	   int verbose)
{
  int i,j;
  int ilo,ihi,inhi;
  float ysave,ytry;
  float rtol;
  float sum,swap;
  float *psum;
  float besty;
  
  if(!(psum = new_array(1,ndim)))
    return 1;
  *nfunk=0;
  for (j=0;j<ndim;j++) {
    for (sum=0.0,i=0;i<=ndim;i++) sum+= p[i*ndim+j];
    psum[j]=sum;
  }
  besty = y[0];

  while(1) {
    ilo = 0;
    ihi = y[0]>y[1] ? (inhi=1,0) : (inhi=0,1);
    
    for (i=0;i<=ndim;i++) {
      if (y[i] <= y[ilo])
	ilo = i;
      if (y[i] > y[ihi]) {
	inhi = ihi;
	ihi = i;
      }
      else if (y[i] > y[inhi] && i != inhi)
	inhi = i;
    }
    if(*nfunk == 0)
      if(verbose)
	printf("Iteration %d: Chisq = %f\n",*nfunk + 1,besty);
    if (y[ilo] < besty) {
      besty = y[ilo];
      if(verbose)
	printf("Iteration %d: Chisq = %f\n",*nfunk + 1,besty);
    }
    rtol = 2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo]));
    if (rtol < ftol) {
      swap = y[0];
      y[0] = y[ilo];
      y[ilo] = swap;
      for(i=0;i<ndim;i++) {
	swap = p[i];
	p[i] = p[ilo*ndim+i];
	p[ilo*ndim+i] = swap;
      }
      if(verbose)
	printf("Within %g of best fit.\n",ftol);
      return 0;
    }

    if (*nfunk>NMAX) {
      swap = y[0];
      y[0] = y[ilo];
      y[ilo] = swap;
      for(i=0;i<ndim;i++) {
	swap = p[i];
	p[i] = p[ilo*ndim+i];
	p[ilo*ndim+i] = swap;
      }
      if(verbose)
	printf("Exceeded max iterations.\n");
      return 0;
    }
    *nfunk += 2;

    ytry = amotry(p,y,psum,ndim,funk,ihi,-1.0,info);

    if(ytry <= y[ilo])
      ytry = amotry(p,y,psum,ndim,funk,ihi,2.0,info);
    else if (ytry >= y[inhi]) {
      ysave = y[ihi];
      ytry = amotry(p,y,psum,ndim,funk,ihi,0.5,info);
      if (ytry >=ysave) {
	for(i=0;i<=ndim;i++) {
	  if (i != ilo) {
	    for (j=0;j<ndim;j++)
	      p[i*ndim+j]=psum[j]=0.5*(p[i*ndim+j]+p[ilo*ndim+j]);
	    y[i] = (*funk)(psum,info);
	  }
	}
	*nfunk += ndim;
	for (j=0;j<ndim;j++) {
	  for (sum=0.0,i=0;i<=ndim;i++) sum+= p[i*ndim+j];
	  psum[j]=sum;
	}
      }
    }
    else
      --(*nfunk);
  }

  psum = del_array(psum);
}

/*
 * Function amotry
 *
 * A sub-function of amoeba.  This function extrapolates by a factor
 *   fac through the face of the simplex across from the high point,
 *   tries the function there, and replaces the high point if the
 *   new point is better.
 *
 */

float amotry(float *p, float y[], float psum[], int ndim,
	     float (*funk)(float [], void *), int ihi, float fac,
	     void *info)

{
  int j;
  float fac1,fac2;
  float ytry,*ptry;

  if(!(ptry=new_array(1,ndim)))
    return -999.9;

  fac1=(1.0 - fac)/ndim;
  fac2 = fac1 - fac;
  for(j=0;j<ndim;j++)
    ptry[j]=psum[j]*fac1 - p[ihi*ndim+j]*fac2;
  ytry = (*funk)(ptry,info);
  if (ytry<y[ihi]) {
    y[ihi] = ytry;
    for(j=0;j<ndim;j++) {
      psum[j] += ptry[j]-p[ihi*ndim+j];
      p[ihi*ndim+j]=ptry[j];
    }
  }

  ptry = del_array(ptry);

  return ytry;
}

