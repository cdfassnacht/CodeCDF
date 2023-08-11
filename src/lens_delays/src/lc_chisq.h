#ifndef lc_chisq_h
#define lc_chisq_h

#include "structdef.h"
#include "lc_setup.h"
#include "lc_funcs.h"

/*.......................................................................
 *
 * Structure definitions
 *
 */

typedef struct {
  float tau;       /* Delay */
  float mu;        /* Magnification */
  float chisq;     /* Chisq for this (delay,magnification) gridpoint */
  int ndof;        /* Number of degrees of freedom */
  float rchisq;    /* Reduced chisq */
  int arraypos;    /* Grid or array position -- only occasionally used */
} LCchisq;         /* Chisq information produced by delays.c */

/*.......................................................................
 *
 * Function declarations
 *
 */

LCchisq *new_lcchisq(int size);
LCchisq *del_lcchisq(LCchisq *lcchisq);
LCchisq *read_lcchisq(char *filename, Axisinfo *tau, Axisinfo *mu);

int do_chi(Fluxrec *raw[], Fluxrec *interp[], int nraw, int nint, 
	   Setup *setup, FILE *logfp, int doprint, int speed);
int shift_chi(Fluxrec *control, Fluxrec *comp, int nctrl, int ncomp, 
	      float tau0, float mu0, Setup *setup, int speed, 
	      LCchisq *bestlag, char *outname, int doprint);
int step_delays(Fluxrec *control, Fluxrec *comp, int nctrl, int ncomp, 
		float ratio, Prange tau, LCchisq *chisq, 
		Setup *setup, int speed, FILE *ofp, int doprint);
LCchisq *find_min_chisq(LCchisq *chisq, int npoints, int fitparabola);
int print_chisq_slice(LCchisq *chisq, int size, char *outname);


#endif
