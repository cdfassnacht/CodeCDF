#ifndef lc_interp_h
#define lc_interp_h

#include "lc_funcs.h"
#include "lc_setup.h"

#define DAYSTEP 1.0    /* Stepsize for creating interpolated light curves */

/*.......................................................................
 *
 * Function declarations
 *
 */

Fluxrec *interp_switch(Fluxrec *raw, int nraw, Setup *setup, int nbad);
Fluxrec *lin_interp(Fluxrec *data, int ndata, int ngrid, float dx, int nbad);
Fluxrec *csmooth(Fluxrec *raw, int nraw, Setup *setup);
Fluxrec *smooth_in_place(Fluxrec *raw, int nraw, Setup *setup);
int boxcar(Fluxrec *raw, int nraw, Fluxrec *smooth, int nsmooth, float width, 
	   float step);
int medsmooth(Fluxrec *raw, int nraw, Fluxrec *smooth, int nsmooth, 
	      float width, float step);
int triangle(Fluxrec *raw, int nraw, Fluxrec *smooth, int nsmooth, float width, 
	     float step);
int gaussian(Fluxrec *raw, int nraw, Fluxrec *smooth, int nsmooth, float sigma, 
	     float step);
int varbox(Fluxrec *raw, int nraw, Fluxrec *smooth, int nsmooth, int nvar, 
	   float step);
int vartri(Fluxrec *raw, int nraw, Fluxrec *smooth, int nsmooth, int nvar, 
	   float step);
int var_wmean(Fluxrec *flux, float *pweight, int npoints, float *sig_wmean);

#endif
