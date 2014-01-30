#ifndef plotfuncs_h
#define plotfuncs_h

#include "structdef.h"

void open_plot_window();
int plot_open(float cheight);
int plot_open_gen(float cheight, char *plotdev);
int plot_open_land();
int plot_open_tall(float cheight, char *plotdev);
int plot_open_spec(float cheight, char *plotdev);
int plot_close();
Pos return_position_from_click(char *prompt, int verbose, int *status);
int plot_xyerr(float *x, float *y, float *err, int npoints,char *xlab,
	       char *ylab, char *title, int flipx, int flipy);
int plot_lcurve(Fluxrec *flux, int npoints, char *xlab, char *ylab, 
		char *title);
int plot_gray(float *array, int size);
int plot_cont(float *array, int size, float stepsize, float zeropt, 
	      int ncont, int color);
int plot_pts(Posinfo *pts, int size, float pixscale, int npts, 
	     int color, int ptype);
int plot_obsim(Posinfo *pi, int size, Modinfo *modinfo);
int plot_labs(char *name, int size, float pixscale);
int plot_spec(Pos *spectrum, int nlines, float minflux, float maxflux, 
	      float minlambda, float maxlambda, int ltype, int dobox,
	      float cheight);
int plot_specpts(Pos *spectrum, int nlines, float minflux, float maxflux, 
		 float minlambda, float maxlambda, int dobox,
		 float cheight);
void plot_specbox(float minlambda, float maxlambda, float minflux, 
		  float maxflux, int dobox, float cheight);

#endif
