#ifndef noninterp_fns_h
#define noninterp_fns_h

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
  float disp;      /* Dispersion for this (delay,magnification) gridpoint */
  int arraypos;    /* Grid or array position -- only occasionally used */
} LCdisp;          /* Dispersion information produced by delays.c */

/*.......................................................................
 *
 * Function declarations
 *
 */

LCdisp *new_lcdisp(int size);
LCdisp *del_lcdisp(LCdisp *lcdisp);
LCdisp *read_lcdisp(char *filename, Axisinfo *tau, Axisinfo *mu);

int call_disp(Fluxrec *flux[], int npoints, Setup *setup,
	      int doprint);
int disp_setup(Fluxrec *flux[], int ncurves, int *npoints, 
	       int *index, Setup *setup, LCdisp *bestdisp, char *outname, 
	       int doprint);
int two_curve_disp(Fluxrec *flux[], int *npoints, int *index,
		       Prange *tau0, Prange *mu0, Setup *setup, 
		       LCdisp *bestdisp, char *outname, int doprint);
int four_curve_disp(Fluxrec *flux[], int *npoints, int *index,
		    Prange *tau0, Prange *mu0, Setup *setup, 
		    LCdisp *bestdisp, char *outname, int doprint);
float disp_d1(Fluxrec *compos, int nccompos, float delta);
float disp_d2(Fluxrec *compos, int nccompos, float delta);
float disp_lovell(Fluxrec *a, Fluxrec *b, int npoints, float tau, float mu,
		  float delta);
int print_disp_slice(LCdisp *disp, int size, char *outname);
int call_dcf(Fluxrec *flux[], int size, char *filename, FILE *logfp, 
	     int doprint);
Fluxrec *discrete_corr(Fluxrec *flux1, Fluxrec *flux2, int npoints, 
		       float binsize, float maxlag, int ndcf);

#endif
