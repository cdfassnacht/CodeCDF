/*
 * plotfuncs.c
 *
 * A collection of functions used to make plots by calling PGPLOT routines.
 *
 * Revision history:
 *  1998Apr27, CDF, Modified plot_pts so that the action of plot_best 
 *                   could be reproduced.  Deleted plot_best.
 *                  Renamed plot_boxes to plot_obsim and changed boxes to
 *                   circles within that function.
 *  1998Aug15, CDF, Added plot_xyerr to plot x, y, and errors on y.
 *  2001Oct15, CDF, Added plot_specbox to plot box surrounding a spectrum
 *                   with several labeling options.
 *  2008Jul16, CDF, Moved open_plot_window function (a generic PGPLOT window
 *   opening function) from match_astrom.c into this library.
 *  2010Dec21, CDF, Added functionality to plot_xyerr
 *                  Added a new function, plot_xy
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cpgplot.h"
#include "structdef.h"
#include "plotfuncs.h"

/*.......................................................................
 *
 * Function open_plot_window
 *
 * Opens a generic pgplot window
 *
 * Inputs: (none)
 *
 * Output: (none)
 *
 */

void open_plot_window()
{
  char *plotdev=NULL;     /* PGPLOT display device */

  plotdev = getenv("PGPLOT_DEV");
  if(!plotdev) {
    printf("\nEnvironment variable PGPLOT_DEV is not set.\n");
    while(cpgopen("?") <= 0)
      fprintf(stderr,"ERROR: Not a supported device.  Try again");
  }
  else {
    printf("\nplotdev = %s\n",plotdev);
    if(cpgopen(plotdev) <= 0) {
      while(cpgopen("?") <= 0)
	fprintf(stderr,"ERROR: Not a supported device.  Try again");
    }
  }
}

/*.......................................................................
 *
 * Function plot_open
 *
 * Opens a PGPLOT plotting device and sets a square viewport, with the
 *  appropriate amount of room for plot labels with a scale of cheight.
 *  The character font is selected to be Roman.
 *
 * Input:  float cheight       character scale height for labels
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int plot_open(float cheight)
{
  if(cpgbeg(0, "?", 1, 1) != 1)
    return 1;
  cpgpap(0.0,1.0);
  cpgpage();
  cpgsch(cheight);
  cpgvstd();
  cpgscf(2);

  return 0;
}

/*.......................................................................
 */

int plot_open_gen(float cheight, char *plotdev)
{
  if(!plotdev) {
    while(cpgbeg(0, "?", 1, 1) != 1)
      fprintf(stderr,"ERROR: Not a supported device.  Try again");
  }
  else {
    if(cpgbeg(0, plotdev, 1, 1) != 1)
      while(cpgbeg(0,"?",1,1) != 1)
	fprintf(stderr,"ERROR: Not a supported device.  Try again");
  }
  cpgpage();
  cpgsch(cheight);
  cpgscf(2);
  cpgvstd();

  return 0;
}

/*.......................................................................
 *         char *plotdev       name of PGPLOT display device
 */

int plot_open_tall(float cheight, char *plotdev)
{
  plotdev = getenv("PGPLOT_DEV");
  if(!plotdev) {
    while(cpgbeg(0, "?", 1, 1) != 1)
      fprintf(stderr,"ERROR: Not a supported device.  Try again");
  }
  else {
    if(cpgbeg(0, plotdev, 1, 1) != 1)
      while(cpgbeg(0,"?",1,1) != 1)
	fprintf(stderr,"ERROR: Not a supported device.  Try again");
  }
  cpgpap(0.0,1.618);
  cpgpage();
  cpgsch(cheight);
  cpgvstd();

  return 0;
}

/*.......................................................................
 *         char *plotdev       name of PGPLOT display device
 */

int plot_open_spec(float cheight, char *plotdev)
{
  if(!plotdev) {
    while(cpgbeg(0, "?", 1, 1) != 1)
      fprintf(stderr,"ERROR: Not a supported device.  Try again");
  }
  else {
    if(cpgbeg(0, plotdev, 1, 1) != 1)
      while(cpgbeg(0,"?",1,1) != 1)
	fprintf(stderr,"ERROR: Not a supported device.  Try again");
  }
  cpgpap(0.0,0.618);
  cpgpage();
  cpgsch(cheight);
  cpgvstd();

  return 0;
}

/*.......................................................................
 */

int plot_open_land()
{
  if(cpgbeg(0, "?", 1, 1) != 1)
    return 1;
  cpgpage();
  cpgvstd();

  return 0;
}

/*.......................................................................
 */

int plot_close()
{
  cpgend();

  return 0;
}

/*.......................................................................
 *
 * Function return_location_from_click
 *
 * Returns the position of the cursor when the mouse is clicked within
 *  the display window
 *
 * Inputs: char *prompt        prompt for mouse click
 *         int verbose         flag set to 1 for verbose output
 *         int *status         flag set to ERROR on error
 *
 * Output: Pos cursloc         cursor location
 */

Pos return_position_from_click(char *prompt, int verbose, int *status)
{
  float cx,cy;            /* Cursor position */
  char cchar='0';         /* Character returned by cursor command */
  Pos cpos;               /* Cursor position */

  printf("%s\n",prompt);
  if(cpgcurs(&cx,&cy,&cchar)>0) {
    cpos.x = cx;
    cpos.y = cy;
    sprintf(cpos.label,"%c",cchar);
    if(verbose)
      printf("Mouse clicked at: %7.2f %7.2f\n",cx,cy);
    *status = SUCCESS;
  }
  else {
    fprintf(stderr,"*** Invalid cursor input ***\n");
    *status = ERROR;
  }

  return cpos;
}

/*.......................................................................
 *
 * Function plot_xy
 *
 * Given input arrays of x and y, draws a box and then plots
 *  points (no errorbars).
 *
 * Inputs: float *x            x array
 *         float *y            y array
 *         int npoints         number of points in curve
 *         char *xlab          label for x axis
 *         char *ylab          label for y axis
 *         char *title         title
 *         int flipx           set to 1 if x-axis should increase to the left
 *         int flipy           set to 1 if y-axis should increase to the left
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int plot_xy(float *x, float *y, int npoints,char *xlab,
	    char *ylab, char *title, int flipx, int flipy)
{
  int i;             /* Looping variable */
  float tmp;         /* Temporary storage variable */
  float xmin,xmax;   /* Min and max x values */
  float ymin,ymax;   /* Min and max y values */
  float size;        /* Temporary difference between min and max values */
  float *xptr;       /* Pointer for navigating x */
  float *yptr;       /* Pointer for navigating y */
  float *eptr;       /* Pointer for navigating err */

  /*
   * Initialize
   */

  xmin = xmax = *x;
  ymin = ymax = *y;

  /*
   * Find min and max values
   */

  for(i=0,xptr=x,yptr=y; i<npoints; i++,xptr++,yptr++) {
    if(*xptr > xmax)
      xmax = *xptr;
    if(*xptr < xmin)
      xmin = *xptr;
    if(*yptr > ymax)
      ymax = *yptr;
    if(*yptr < ymin)
      ymin = *yptr;
  }

  /*
   * Increase size of box a bit
   */

  size = xmax - xmin;
  xmin -= 0.15*size;
  xmax += 0.15*size;
  size = ymax - ymin;
  ymin -= 0.15*size;
  ymax += 0.15*size;

  /*
   * Flip x and/or y axis if requested
   */

  if(flipx) {
    tmp = xmax;
    xmax = xmin;
    xmin = tmp;
  }
  if(flipy) {
    tmp = ymax;
    ymax = ymin;
    ymin = tmp;
  }

  /*
   * Draw box and label axes
   */

  cpgenv(xmin,xmax,ymin,ymax,0,0);
  cpglab(xlab,ylab,title);

  /*
   * Put in the points
   */

  cpgpt(npoints,x,y,16);

  return 0;
}

/*.......................................................................
 *
 * Function plot_xyerr
 *
 * Given input arrays of x, y, and y_err, draws a box and then plots
 *  points with errorbars.
 *
 * Inputs: float *x            x array
 *         float *y            y array
 *         float *err          error array
 *         int npoints         number of points in curve
 *         char *xlab          label for x axis
 *         char *ylab          label for y axis
 *         char *title         title
 *         int flipx           set to 1 if x-axis should increase to the left
 *         int flipy           set to 1 if y-axis should increase to the left
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int plot_xyerr(float *x, float *y, float *err, int npoints,char *xlab,
	       char *ylab, char *title, int flipx, int flipy)
{
  int i;             /* Looping variable */
  float tmp;         /* Temporary storage variable */
  float xmin,xmax;   /* Min and max x values */
  float ymin,ymax;   /* Min and max y values */
  float size;        /* Temporary difference between min and max values */
  float *xptr;       /* Pointer for navigating x */
  float *yptr;       /* Pointer for navigating y */
  float *eptr;       /* Pointer for navigating err */

  /*
   * Initialize
   */

  xmin = xmax = *x;
  ymin = ymax = *y;

  /*
   * Find min and max values
   */

  for(i=0,xptr=x,yptr=y; i<npoints; i++,xptr++,yptr++) {
    if(*xptr > xmax)
      xmax = *xptr;
    if(*xptr < xmin)
      xmin = *xptr;
    if(*yptr > ymax)
      ymax = *yptr;
    if(*yptr < ymin)
      ymin = *yptr;
  }

  /*
   * Increase size of box a bit
   */

  size = xmax - xmin;
  xmin -= 0.15*size;
  xmax += 0.15*size;
  size = ymax - ymin;
  ymin -= 0.15*size;
  ymax += 0.15*size;

  /*
   * Flip x and/or y axis if requested
   */

  if(flipx) {
    tmp = xmax;
    xmax = xmin;
    xmin = tmp;
  }
  if(flipy) {
    tmp = ymax;
    ymax = ymin;
    ymin = tmp;
  }

  /*
   * Draw box and label axes
   */

  cpgenv(xmin,xmax,ymin,ymax,0,0);
  cpglab(xlab,ylab,title);

  /*
   * Put in the points
   */

  cpgpt(npoints,x,y,16);
  for(i=0,xptr=x,yptr=y,eptr=err; i<npoints; i++,xptr++,yptr++,eptr++) {
    ymin = *yptr - *eptr;
    ymax = *yptr + *eptr;
    cpgerry(1,xptr,&ymax,&ymin,1.0);
  }

  return 0;
}

/*.......................................................................
 *
 * Function plot_lcurve
 *
 * Given input lightcurve, draws a box and then plots points with errorbars.
 *
 * Inputs: Fluxrec *flux       input lightcurve
 *         int npoints         number of points in curve
 *         char *xlab          label for x axis
 *         char *ylab          label for y axis
 *         char *title         title
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int plot_lcurve(Fluxrec *flux, int npoints, char *xlab, char *ylab, 
		char *title)
{
  int i;             /* Looping variable */
  float xmin,xmax;   /* Min and max x values */
  float ymin,ymax;   /* Min and max y values */
  float size;        /* Temporary difference between min and max values */
  Fluxrec *fptr;     /* Pointer for navigating flux */

  /*
   * Initialize
   */

  xmin = xmax = flux->day;
  ymin = ymax = flux->flux;

  /*
   * Find min and max values
   */

  for(i=0,fptr=flux; i<npoints; i++,fptr++) {
    if(fptr->day > xmax)
      xmax = fptr->day;
    if(fptr->day < xmin)
      xmin = fptr->day;
    if(fptr->flux > ymax)
      ymax = fptr->flux;
    if(fptr->flux < ymin)
      ymin = fptr->flux;
  }

  /*
   * Increase size of box a bit
   */

  size = xmax - xmin;
  xmin -= 0.15*size;
  xmax += 0.15*size;
  size = ymax - ymin;
  ymin -= 0.15*size;
  ymax += 0.15*size;

  /*
   * Draw box and label axes
   */

  cpgenv(xmin,xmax,ymin,ymax,0,0);
  cpglab(xlab,ylab,title);

  /*
   * Put in the points
   */

  for(i=0,fptr=flux; i<npoints; i++,fptr++) {
    ymin = fptr->flux - fptr->err;
    ymax = fptr->flux + fptr->err;
    cpgpt(1,&fptr->day,&fptr->flux,16);
    cpgerry(1,&fptr->day,&ymax,&ymin,1.0);
  }

  return 0;
}

/*.......................................................................
 */

int plot_gray(float *array,int size)
{
  static float tr[6] = {0,1,0,0,0,1};  
  int i;
  float hi,lo;
  float *arrayptr;

  arrayptr = array;
  hi=lo=*arrayptr;
  for(i = 0;i < size*size; i++,arrayptr++){
    if(hi < *arrayptr)
      hi = *arrayptr;
    if(lo > *arrayptr)
      lo = *arrayptr;
  }

  /* 
   * Plot the greyscale
   */

  cpgwnad(1.0*size,1.0,1.0,1.0*size);
  cpggray(array,size,size,1,size,1,size,hi,lo,tr);
  
  return 0;
}

/*.......................................................................
 *
 * Function plot_cont
 *
 * Plots contours levels for a given array of data values.
 *
 * Inputs: float *array        data array
 *         int size            if array is (n x n) then size = n
 *         float stepsize      steps between contour levels
 *         float zeropt        lowest contouring level
 *         int ncont           number of contours
 *         int color           color of contour lines
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int plot_cont(float *array, int size, float stepsize, float zeropt,
	      int ncont, int color)
{
  static float tr[6] = {0,1,0,0,0,1};  
  int i;
  float hi,lo;
  float *cont;
  float *arrayptr;

  arrayptr = array;
  hi=lo=*arrayptr;
  for(i = 0;i < size*size; i++){
    if(hi < *arrayptr)
      hi = *arrayptr;
    if(lo > *arrayptr)
      lo = *arrayptr;
    arrayptr++;
  }

  if (!(cont = new_array(ncont,1)))
    return 1;
  arrayptr = cont;
  if(zeropt != -999.9)
    lo = zeropt;
  if(stepsize == 0.0)
    for(i=0; i<ncont; i++)
      *arrayptr++ = lo + i*0.005*(hi-lo);
  else
    for(i=0; i<ncont; i++)
      *arrayptr++ = lo + stepsize*(1+i);

  cpgsci(color);
  cpgwnad(1.0*size,1.0,1.0,1.0*size);
  cpgsls(2);
  cpgslw(2);
  cpgcont(array,size,size,1,size,1,size,cont,-ncont,tr);
  cpgsci(1);
  cpgsls(1);
  cpgslw(1);

  cont = del_array(cont);
  return 0;
}

/*.......................................................................
 *
 * Function plot_pts
 *
 * Plots a series of points, which is passed to the function in a Posinfo
 *  array.  The color and type of the points desired are passed to the 
 *  function.
 *
 * Inputs: Posinfo *pts        array containing points
 *         int size            size of surface on which points are
 *                              plotted
 *         float pixscale      pixel scale
 *         int npts            size of array of points
 *         int color           color of points
 *         int ptype           type of point to be plotted (PGPLOT)
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int plot_pts(Posinfo *pts, int size, float pixscale, int npts,
	     int color, int ptype)
{
  int i;
  Posinfo *ptsptr;

  ptsptr = pts;
  cpgsci(color);
  cpgwnad(size*pixscale/2.0,-size*pixscale/2.0,
	  -size*pixscale/2.0,size*pixscale/2.0);
  for(i=0;i<npts;i++,ptsptr++)
    cpgpt(1,&ptsptr->x,&ptsptr->y,ptype);
  cpgsci(1);

  return 0;
}

/*.......................................................................
 *
 * Function plot_obsim
 *
 * Plots circles  at the positions of the observed images, as well as
 *  marking the position of the lens and source(s).
 *
 * Inputs: Posinfo *pi         image positions
 *         int size            size of grid on which boxes are plotted
 *         Modinfo *modinfo    other model, source and image info
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int plot_obsim(Posinfo *pi, int size, Modinfo *modinfo)
{
  int i;            /* Looping variable */
  Srcinfo *pptr;    /* Pointer to navigate modinfo->fsource */
  Posinfo *ptsptr;  /* Pointer to navigate pi */

  /*
   * Set window size in angular units
   */

  cpgwnad(size*modinfo->pixscale/2.0,-size*modinfo->pixscale/2.0,
	  -size*modinfo->pixscale/2.0,size*modinfo->pixscale/2.0);

  /*
   * Plot all source positions
   */

  cpgsfs(2);
  cpgslw(2);
  for(i=0,pptr=modinfo->fsource; i<modinfo->nfsrc; i++,pptr++)
    cpgpt(1,&(pptr->x.par),&(pptr->y.par),-4);
  cpgslw(1);

  /*
   * Plot the lens position
   */

  cpgsci(11);
  cpgpt(1,&(modinfo->lp->xl.par),&(modinfo->lp->yl.par),8);
  cpgsci(1);

  /*
   * Plot the observed image positions
   */

  cpgslw(1);
  cpgsci(5);
  for(i=0,ptsptr=pi; i<modinfo->nobs; i++,ptsptr++)
    cpgcirc(ptsptr->x,ptsptr->y,10*modinfo->pixscale);
  cpgsci(1);

  return 0;
}

/*.......................................................................
 *
 * Function plot_labs
 *
 * Draws the border for the plot and then labels it.
 *
 * Inputs: char *name          title of plot
 *         int size            size of plot in pixels
 *         float pixscale      size of pixels in angular units
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int plot_labs(char *name, int size, float pixscale)
{

  /*
   * Set font, character height, and line width
   */

  cpgscf(2);
  cpgsch(1.5);
  cpgslw(2);

  /*
   * Set size of box in angular units
   */

  cpgwnad(size*pixscale/2.0,-size*pixscale/2.0,
	  -size*pixscale/2.0,size*pixscale/2.0);

  /*
   * Draw box and label it
   */

  cpgbox("BCNST", 0, 0, "BCNST", 0, 0);
  cpglab("\\gD x (arcsec)","\\gD y (arcsec)",name);

  /*
   * Reset values
   */

  cpgwnad(1.0*size,1.0,1.0,1.0*size);
  cpgscf(1);
  cpgsch(1.0);
  cpgslw(1);

  return 0;
}

/*.......................................................................
 *
 * Function plot_spec
 *
 * Plots a spectrum.  The various inputs are used to control the plotting
 *  of the box surrounding the spectrum and the line type used in the
 *  plotting.
 *
 * Inputs: Pos *spectrum       array of wavelengths and fluxes
 *         int nlines          number of lines in array
 *         float minflux       min value of y axis
 *         float maxflux       max value of y axis
 *         float minlambda     min value of x axis
 *         float maxlambda     max value of x axis
 *         int ltype           line type used for plotting
 *         int dobox           flag that determines box labeling style
 *                              0 ==> no box
 *                              1 ==> fully-labeled box
 *                              2 ==> y labels only
 *                              3 ==> x labels only
 *                              4 ==> no bottom axis, y labels only
 *                              5 ==> no top axis, y labels only
 *                              6 ==> no top axis, x labels only
 *         float cheight       character height
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error.
 *
 * v12Oct01 CDF, Moved splitting of spectrum (a Pos array) into x and y
 *                components into a call to the new function pos2xy (in the
 *                structdef library).
 * v15Oct01 CDF, Moved drawing of box into new plot_specbox function below.
 */

int plot_spec(Pos *spectrum, int nlines, float minflux, float maxflux, 
	      float minlambda, float maxlambda, int ltype, int dobox,
	      float cheight) 
{
  int i;               /* Looping variable */
  int ngood=0;         /* Number of good points in spectrum */
  int no_error=1;      /* Flag set to 0 on error */
  float **xy={NULL};   /* Temp container for lambda and flux arrays */

  /*
   * Convert spectrum to flux and lambda arrays
   */

  if(!(xy = pos2xy(spectrum,nlines,&ngood))) {
    fprintf(stderr,"ERROR: plot_spec\n");
    return 1;
  }


  if(ngood == 0) {
    fprintf(stderr,"ERROR: plot_spec.  No valid points.\n");
    no_error = 0;
  }

  /*
   * Do the plotting
   */

  if(no_error) {
    plot_specbox(minlambda,maxlambda,minflux,maxflux,dobox,cheight);
    cpgsls(ltype);
    cpgbin(nlines,xy[0],xy[1],1);
    cpgsls(1);
  }

  /*
   * Clean up and exit
   */

  xy = del_float2(xy,2);

  if(no_error)
    return 0;
  else {
    fprintf(stderr,"ERROR: plot_spec\n");
    return 1;
  }
}
  

/*.......................................................................
 *
 * Function plot_specpts
 *
 * Plots a spectrum as a series of discrete points.  The various inputs 
 *  are used to control the plotting of the box surrounding the spectrum 
 *  and the line type used in the plotting.
 *
 * Inputs: Pos *spectrum       array of wavelengths and fluxes
 *         int nlines          number of lines in array
 *         float minflux       min value of y axis
 *         float maxflux       max value of y axis
 *         float minlambda     min value of x axis
 *         float maxlambda     max value of x axis
 *         int dobox           flag that determines box labeling style
 *                              0 ==> no box
 *                              1 ==> fully-labeled box
 *                              2 ==> y labels only
 *                              3 ==> x labels only
 *                              4 ==> no bottom axis, y labels only
 *                              5 ==> no top axis, y labels only
 *                              6 ==> no top axis, x labels only
 *         float cheight       character height
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error.
 *
 * v12Oct01 CDF, Moved splitting of spectrum (a Pos array) into x and y
 *                components into a call to the new function pos2xy (in the
 *                structdef library).
 * v15Oct01 CDF, Moved drawing of box into new plot_specbox function below.
 */

int plot_specpts(Pos *spectrum, int nlines, float minflux, float maxflux, 
		 float minlambda, float maxlambda, int dobox,
		 float cheight) 
{
  int i;               /* Looping variable */
  int ngood=0;         /* Number of good points in the array */
  int no_error=1;      /* Flag set to 0 on error */
  float **xy={NULL};   /* Temp container for lambda and flux arrays */

  /*
   * Convert spectrum to flux and lambda arrays
   */

  if(!(xy = pos2xy(spectrum,nlines,&ngood))) {
    fprintf(stderr,"ERROR: plot_spec\n");
    return 1;
  }


  if(ngood == 0) {
    fprintf(stderr,"ERROR: plot_specpts.  No valid points.\n");
    no_error = 0;
  }

  /*
   * Do the plotting
   */

  if(no_error) {
    plot_specbox(minlambda,maxlambda,minflux,maxflux,dobox,cheight);
    cpgpt(ngood,xy[0],xy[1],1);
    cpgsls(1);
  }

  /*
   * Clean up and exit
   */

  xy = del_float2(xy,2);

  if(no_error)
    return 0;
  else {
    fprintf(stderr,"ERROR: plot_specpts\n");
    return 1;
  }
}
  
/*.......................................................................
 *
 * Function plot_specbox
 *
 * Plots the box inside which the spectrum will be plotted.  There are
 *  several labeling options, which are selected via the passed dobox
 *  variable.
 *
 * Inputs: float minflux       min value of y axis
 *         float maxflux       max value of y axis
 *         float minlambda     min value of x axis
 *         float maxlambda     max value of x axis
 *         int dobox           flag that determines box labeling style
 *                              0 ==> no box
 *                              1 ==> fully-labeled box
 *                              2 ==> y labels only
 *                              3 ==> x labels only
 *                              4 ==> no bottom axis, y labels only
 *                              5 ==> no top axis, y labels only
 *                              6 ==> no top axis, x labels only
 *         float cheight       character height
 *
 * Output: (none)
 */

void plot_specbox(float minlambda, float maxlambda, float minflux, 
		  float maxflux, int dobox, float cheight)
{
    cpgslw(2);
    cpgsch(cheight);
    cpgswin(minlambda,maxlambda,minflux,maxflux);
    switch(dobox) {
    case 1:
      cpgbox("BCNST", 0, 0, "BCNSTV", 0, 0);
      break;
    case 2:
      cpgbox("BCST", 0, 0, "BCNSTV", 0, 0);
      break;
    case 3:
      cpgbox("BCNST", 0, 0, "BCST", 0, 0);
      break;
    case 4:
      cpgbox("CST", 0, 0, "BCNSTV", 0, 0);
      break;
    case 5:
      cpgbox("BST", 0, 0, "BCNSTV", 0, 0);
      break;
    case 6:
      cpgbox("BNST", 0, 0, "BC", 0, 0);
      break;
    default:
      cpgbox("BCST", 0, 0, "BCST", 0, 0);
      break;
    }
    cpgslw(1);
}
