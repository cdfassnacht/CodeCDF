#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "structdef.h"
#include "list_tools.h"
#include "plotfuncs.h"
#include "cpgplot.h"
#include "fitting.h"
#include "lensprogs.h"

#define BERR 1e-7
#define CERR 0.01
#define MAXITER 500
#define RAD2ASEC 206265.0  /* To convert radians to arcsec */
#define H 3.23e-18         /* Because H=100 h km/sec/Mpc = 3.2e-18 h 1/sec */
#define C 3.0e10           /* Speed of light in cm/sec */
#define TINY 1.0e-18       /* To avoid division by zero */

/*......................................................................
 *
 * Function sis
 *
 * Calculates the potential, gradient of the potential and second
 *   derivatives of the potential at a list of positions using the
 *   singular isothermal sphere model.  These values are added to
 *   the previously calculated part of the potential, if any.
 *
 * Inputs:  Poten *poten       holder of the potential to be modified
 *          Posinfo *posinfo   list of positions on the lens plane
 *          Lensparams *params lens potential  parameters
 *          int n1,n2          size of the position array
 *          int fast           if 1 then only calculate the first derivs.
 *
 * Output: No output - *poten is modified as a result of this function.
 *
 */

void sis(Poten *poten, Posinfo *posinfo, Lensparams *lp, 
	 int n1, int n2, int fast)
{
  int i;
  float x,y,r,r3;
  float b;
  Posinfo *posptr;
  Poten *tmpptr;

  b = lp->b.par;
  posptr = posinfo;
  tmpptr = poten;
  for(i=0;i<n1*n2;i++,posptr++,tmpptr++) {
    x = posptr->x - lp->xl.par;
    y = posptr->y - lp->yl.par;
    r = sqrt(x*x + y*y);
    r3 = r*r*r;
    tmpptr->dphdx += b*x/r;
    tmpptr->dphdy += b*y/r;
    if(!fast) {
      tmpptr->ph += b*r;
      tmpptr->dphdxx += b/r - b*x*x/r3;
      tmpptr->dphdyy += b/r - b*y*y/r3;
      tmpptr->dphdxy += -b*x*y/r3;
    }
  }
}

/*......................................................................
 *
 * Function point
 *
 * Calculates the potential, gradient of the potential and second
 *   derivatives of the potential at a list of positions using the
 *   singular point mass model.  These values are added to
 *   the previously calculated part of the potential, if any.
 *
 * Inputs:  Poten *poten       holder of the potential to be modified
 *          Posinfo *posinfo   list of positions on the lens plane
 *          Lensparams *params lens potential  parameters
 *          int n1,n2          size of the position array
 *          int fast           if 1 then only calculate the first derivs.
 *
 * Output: No output - *poten is modified as a result of this function.
 *
 */

void point(Poten *poten, Posinfo *posinfo, Lensparams *lp, 
	   int n1, int n2, int fast)
{
  int i;
  float x,y,r,r2,r4;
  float b,b2;
  Posinfo *posptr;
  Poten *tmpptr;

  b = lp->b.par;
  b2 = b*b;
  posptr = posinfo;
  tmpptr = poten;
  for(i=0;i<n1*n2;i++,posptr++,tmpptr++) {
    x = posptr->x - lp->xl.par;
    y = posptr->y - lp->yl.par;
    r = sqrt(x*x + y*y);
    r2 = r*r;
    r4 = r2*r2;
    tmpptr->dphdx += b2*x/r2;
    tmpptr->dphdy += b2*y/r2;
    if(!fast) {
      tmpptr->ph += b2*log(r);
      tmpptr->dphdxx += b2/r2 - 2*b2*x*x/r4;
      tmpptr->dphdyy += b2/r2 - 2*b2*y*y/r4;
      tmpptr->dphdxy += -2*b2*x*y/r4;
    }
  }
}


/*......................................................................
 *
 * Function shear_int
 *
 * Calculates the potential, gradient of the potential and second
 *   derivatives of the potential at a list of positions using the
 *   internal shear model.  These values are added to
 *   the previously calculated part of the potential, if any.
 *
 * Inputs:  Poten *poten       holder of the potential to be modified
 *          Posinfo *posinfo   list of positions on the lens plane
 *          Lensparams *params lens potential  parameters
 *          int n1,n2          size of the position array
 *          int fast           if 1 then only calculate the first derivs.
 *
 * Output: No output - *poten is modified as a result of this function.
 *
 */

void shear_int(Poten *poten, Posinfo *posinfo, Lensparams *lp, 
	       int n1, int n2, int fast)
{
  int i;
  float x,y,r,r4,r6,r8;
  float sqdif,cross,angterm;
  float b,b4,gamma;
  float A,B;
  Posinfo *posptr;
  Poten *tmpptr;

  b = lp->b.par;
  b4 = b*b*b*b;
  gamma = lp->gamma.par;
  A = cos(2*lp->theta_g.par);
  B = sin(2*lp->theta_g.par);
  posptr = posinfo;
  tmpptr = poten;
  for(i=0;i<n1*n2;i++,posptr++,tmpptr++) {
    x = posptr->x - lp->xl.par;
    y = posptr->y - lp->yl.par;
    r = sqrt(x*x + y*y);
    r4 = r*r*r*r;
    r6 = r4*r*r;
    r8 = r4*r4;
    sqdif = A*(x*x-y*y);
    cross = 2*B*x*y;
    angterm = sqdif+cross;
    tmpptr->dphdx += 2*gamma*b4*x*angterm/r6 - gamma*b4*(A*x+B*y)/r4;
    tmpptr->dphdy += 2*gamma*b4*y*angterm/r6 + gamma*b4*(A*y-B*x)/r4;
    if(!fast) {
      tmpptr->ph += gamma*b4/(2*r4)*angterm;
      tmpptr->dphdxx += -2*gamma*b4*(4*A*x*x+sqdif+3*cross)/r6 - gamma*A*b4/r4
	- 12*gamma*b4*x*x*angterm/r8;
      tmpptr->dphdyy += -2*gamma*b4*(sqdif-4*A*y*y+3*cross)/r6 - gamma*A*b4/r4
	- 12*gamma*b4*y*y*angterm/r8;
      tmpptr->dphdxy += 3*gamma*B*b4/r4 - 12*gamma*b4*x*y*angterm/r8;
    }
  }
}

/*......................................................................
 *
 * Function shear_mixed
 *
 * Calculates the potential, gradient of the potential and second
 *   derivatives of the potential at a list of positions using the
 *   mixed shear model.  These values are added to
 *   the previously calculated part of the potential, if any.
 *
 * Inputs:  Poten *poten       holder of the potential to be modified
 *          Posinfo *posinfo   list of positions on the lens plane
 *          Lensparams *params lens potential  parameters
 *          int n1,n2          size of the position array
 *          int fast           if 1 then only calculate the first derivs.
 *
 * Output: No output - *poten is modified as a result of this function.
 *
 */

void shear_mixed(Poten *poten, Posinfo *posinfo, Lensparams *lp, 
		 int n1, int n2, int fast)
{
  int i;
  float x,y,r,r3,r5;
  float sqdif,cross,angterm;
  float b,b4,gamma;
  float A,B;
  Posinfo *posptr;
  Poten *tmpptr;

  b = lp->b.par;
  b4 = b*b*b*b;
  gamma = lp->gamma.par;
  A = cos(2*lp->theta_g.par);
  B = sin(2*lp->theta_g.par);
  posptr = posinfo;
  tmpptr = poten;
  for(i=0;i<n1*n2;i++,posptr++,tmpptr++) {
    x = posptr->x - lp->xl.par;
    y = posptr->y - lp->yl.par;
    r = sqrt(x*x + y*y);
    r3 = r*r*r;
    r5 = r3*r*r;
    sqdif = A*(x*x-y*y);
    cross = 2*B*x*y;
    angterm = sqdif+cross;
    tmpptr->dphdx += -gamma*b*x*angterm/r3 + 2*gamma*b*(A*x+B*y)/r;
    tmpptr->dphdy += -gamma*b*y*angterm/r3 - 2*gamma*b*(A*y-B*x)/r;
    if(!fast) {
      tmpptr->ph += gamma*b*angterm/r;
      tmpptr->dphdxx += 2*gamma*b*A/r - gamma*b*(4*A*x*x+sqdif+3*cross)/r3
	+ 3*gamma*b*x*x*angterm/r5;
      tmpptr->dphdyy += -2*gamma*b*A/r - gamma*b*(sqdif-4*A*y*y+3*cross)/r3
	+ 3*gamma*b*y*y*angterm/r5;
      tmpptr->dphdxy += 3*gamma*b*x*y*angterm/r5;
    }
  }
}

/*......................................................................
 *
 * Function shear_ext
 *
 * Calculates the potential, gradient of the potential and second
 *   derivatives of the potential at a list of positions using the
 *   external shear model.  These values are added to
 *   the previously calculated part of the potential, if any.
 *
 * Inputs:  Poten *poten       holder of the potential to be modified
 *          Posinfo *posinfo   list of positions on the lens plane
 *          Lensparams *params lens potential  parameters
 *          int n1,n2          size of the position array
 *          int fast           if 1 then only calculate the first derivs.
 *
 * Output: No output - *poten is modified as a result of this function.
 *
 */

void shear_ext(Poten *poten, Posinfo *posinfo, Lensparams *lp, 
	       int n1, int n2, int fast)
{
  int i;
  float x,y,r;
  float sqdif,cross,angterm;
  float b,b4,gamma;
  float A,B;
  Posinfo *posptr;
  Poten *tmpptr;

  b = lp->b.par;
  b4 = b*b*b*b;
  gamma = lp->gamma.par;
  A = cos(2*lp->theta_g.par);
  B = sin(2*lp->theta_g.par);
  posptr = posinfo;
  tmpptr = poten;
  for(i=0;i<n1*n2;i++,posptr++,tmpptr++) {
    x = posptr->x - lp->xl.par;
    y = posptr->y - lp->yl.par;
    r = sqrt(x*x + y*y);
    sqdif = A*(x*x-y*y);
    cross = 2*B*x*y;
    angterm = sqdif+cross;
    tmpptr->dphdx += gamma*(A*x+B*y);
    tmpptr->dphdy += -gamma*(A*y-B*x);
    if(!fast) {
      tmpptr->ph += gamma*angterm/2;
      tmpptr->dphdxx += gamma*A;
      tmpptr->dphdyy += -gamma*A;
      tmpptr->dphdxy += gamma*B;
    }
  }
}

