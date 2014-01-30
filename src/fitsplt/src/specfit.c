/* specfit.c
 *
 * Usage: specfit filename
 *
 * Description:  Takes a two-dimension optical spectrum, stored in fits
 *                format, and fits a double Gaussian to cuts along the spatial
 *                dimension.
 *
 * To compile use the appropriate makefile in this directory.
 *
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "libfits.h"
#include "cpgplot.h"
#include "structdef.h"
#include "plotfuncs.h"
#include "fitsim.h"

/*.......................................................................
 *
 * Function declarations
 *
 */

Pos *get_trace(Image *image);
int clip_trace(Pos *trace, int ntrace);
int find_wmean(Pos *slice, int npoints, float m, float b, int ysmax,
	       float *ymean);

/*.......................................................................
 *
 * Main program
 *
 */

int main(int argc, char *argv[])
{
  int i;
  int no_error=1;       /* Flag set to 0 on error */
  float slmax;          /* Max value along a slice */
  float b1lo,b1hi;      /* Definition for background region 1 */
  float b2lo,b2hi;      /* Definition for background region 2 */
  float m,b;            /* Parameters for linear fit to bkgd */
  float *x=NULL;
  float *y=NULL;
  float ymean;
  float smin,smax;      /* Minimum and maximum along a slice */
  double ysmax;         /* Y value associated with slmax */
  char *file_name;      /* The name of the FITS file to load */
  char *setupname;      /* Name of the optional setup file */
  Image *image=NULL;    /* The container of the FITS image array */
  Setup *setup=NULL;    /* Container for the image display info */
  Pos **xslice={NULL};  /* 1D vertical slices */
  Pos **yslice={NULL};  /* 1D horizontal slices */
  Pos *slice=NULL;      /* A 1D slice through the image */
  Pos *trace=NULL;      /* Y position of the trace */
  Pos *pptr;            /* Pointer to navigate Pos arrays */

  /*
   * Check the command line
   */

  if(argc < 2 || argc > 3) {
    setup_help();
    return 1;
  }
  else {
    fitsplt_copyright_notice();
  }

  /*
   * Get the name of the file.
   */

  file_name = argv[1];

  /*
   * Read the data array.
   */

  image = new_Image(file_name);
  if(!image) {
    fprintf(stderr,"ERROR.  Exiting program\n\n");
    return 1;
  }

  /*
   * Print image statistics
   */

  if(imstat(image))
    no_error = 0;

  /*
   *  Allocate memory for the slice arrays
   */

  xslice = (Pos **) malloc(sizeof(Pos *) * image->nx);
  if(!xslice) {
    fprintf(stderr,"ERROR:  Insufficient memory for xslice array.\n");
    no_error = 0;
  }
  yslice = (Pos **) malloc(sizeof(Pos *) * image->ny);
  if(!yslice) {
    fprintf(stderr,"ERROR:  Insufficient memory for yslice array.\n");
    no_error = 0;
  }

  if(!(x = new_array(image->ny,1)))
    no_error = 0;
  if(!(y = new_array(image->ny,1)))
    no_error = 0;

  /*
   * Get the trace
   */

  if(no_error) {
    if(!(trace = get_trace(image)))
      no_error = 0;
  }

  /*
   * Plot initial trace
   */

  plot_open_spec(1.5,"/xs");
  if(no_error)
    if(plot_vslice(trace,0,image->nx-1,smin,smax,1))
       no_error = 0;

  /*
   * Clip the trace
   */

  if(no_error)
    if(clip_trace(trace,image->nx) < 0)
      no_error = 0;

  /*
   * Plot the clipped trace
   */

  if(no_error) {
    cpgpage();
    if(plot_vslice(trace,0,image->nx-1,smin,smax,1))
       no_error = 0;
  }

  /*
   * Plot a vertical cut through the image at the central column
   *  and find its maximum value.
   */

  if(no_error) {
    if(!(slice =  make_vslice(image,image->nx/2.0,0,image->ny-1,
			      &smin,&smax,&ysmax)))
      no_error = 0;
    else {
      cpgpage();
      if(plot_vslice(slice,0,image->ny-1,smin,smax,0))
	no_error = 0;
    }
  }

  /*
   * Replot the central column, expanded around the max
   */

  if(no_error) {
    cpgpage();
    if(plot_vslice(slice,ysmax-40,ysmax+40,smin,smax,0))
       no_error = 0;
  }

  /*
   * Define the background region
   */

  if(no_error)
    if(get_bkgd(&b1lo,&b1hi,&b2lo,&b2hi))
      no_error = 0;

  /*
   * Fit a linear continuum
   */

  if(no_error) {
    m = b = 0.0;
    if(fit_continuum(slice,image->ny,ysmax-20,ysmax-10,ysmax+15,ysmax+25,
		     &m,&b))
      no_error = 0;
    else {
      cpgsci(2);
      for(i=0; i<image->ny; i++) {
	x[i] = i;
	y[i] = b + m * i;
      }
      cpgline(image->ny,x,y);
      cpgsci(1);
    }
    find_wmean(slice,image->ny,m,b,ysmax,&ymean);
  }
#if 0
  /*
   * Slice up the image
   */

  if(no_error)
    if(imslice(image))
      no_error = 0;
#endif

  /*
   * Clean up and exit
   */

  image = del_Image(image);
  setup = del_setup(setup);
#if 0
  for(i=0; i<image->nx; i++) {
    xslice[i] = del_pos(xslice[i]);
  }
#endif
  x = del_array(x);
  y = del_array(y);
  trace = del_pos(trace);
  plot_close();
  printf("\nPlot closed.\n");
  if(xslice)
    free(xslice);
  if(yslice)
    free(yslice);

  if(no_error) {
    printf("\nFinished with program specfit.\n\n");
    return 0;
  }
  else {
    fprintf(stderr,"ERROR. Exiting specfit\n\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function get_trace
 *
 * Gets the "trace", that is the y position of the spectrum of the
 *  source, as a function of the x position.
 *
 * Inputs: Image *image        2D spectrum and associated info.
 *
 * Output: Pos *trace          trace of spectrum
 *
 */

Pos *get_trace(Image *image)
{
  int i;                /* Looping variable */
  int no_error=1;       /* Flag set to 0 on error */
  float smin,smax;      /* Minimum and maximum along a slice */
  Pos *slice=NULL;      /* A 1D slice through the image */
  Pos *trace=NULL;      /* Y position of the trace */
  Pos *pptr;            /* Pointer to navigate Pos arrays */

  /*
   * Allocate memory for trace
   */

  if(!(trace = new_pos(image->nx,1))) {
    fprintf(stderr,"ERROR: get_trace.\n");
    return NULL;
  }

  printf("\nget_trace: Finding trace.\n\n");

  /*
   * Slice up the image vertically.
   */

  if(no_error) {
    for(i=0,pptr=trace; i<image->nx; i++,pptr++) {
      pptr->x = i;
      if(!(slice = make_vslice(image,i,0,image->ny-1,&smin,&smax,
			       &pptr->y)))
	no_error = 0;
      slice = del_pos(slice);
    }
  }

  if(no_error)
    return trace;
  else {
    fprintf(stderr,"ERROR: get_trace.\n");
    return del_pos(trace);
  }
}

/*.......................................................................
 *
 * Function clip_trace
 *
 * Does a sigma-clipping on the trace.
 *
 * Inputs: Pos *trace          spectral trace
 *         int ntrace          number of points in trace
 *
 * Output: int ngood           number of good points after clipping.
 *                             Set to -1 on error.
 *
 */

int clip_trace(Pos *trace, int ntrace)
{
  int i;            /* Looping variable */
  int no_error=1;   /* Flag set to 0 on error */
  int ngood=0;      /* Number of good points in trace */
  float *y=NULL;    /* Temporary container for trace position */
  float *yptr;      /* Pointer to navigate y */
  double rms;       /* RMS of trace */
  double mean;      /* Mean of trace */
  Pos *pptr;        /* Pointer to navigate trace */

  /*
   * Allocate memory for temporary container
   */

  if(!(y = new_array(ntrace,1))) {
    fprintf(stderr,"ERROR: clip_trace.\n");
    return -1;
  }

  /*
   * Load trace into temporary container and calculate mean.
   */

  mean = 0.0;
  for(i=0,yptr=y,pptr=trace; i<ntrace; i++,yptr++,pptr++) {
    *yptr = pptr->y;
    mean += *yptr;
  }

  mean /= ntrace;

  /*
   * Calculate rms
   */

  if((rms = im_rms(y,ntrace,mean)) < 0.0)
    no_error = 0;
  printf("clip_trace: Initial mean=%6.1f, rms=%6.2f\n",mean,rms);

  /*
   * Do the clipping at 3-sigma level
   */

  if(no_error)
    if(sigclip(y,ntrace,&mean,&rms,3.0,0.1) < 0)
      no_error = 0;

  /*
   * Flag the clipped points
   */

  if(no_error) {
    for(i=0,pptr=trace; i<ntrace; i++,pptr++) {
      if(fabs(pptr->y - mean) > 3.0 * rms)
	pptr->flag = 1;
    }
  }

  return 0;
}

/*.......................................................................
 *
 * Function get_bkgd
 *
 * Gets background regions.
 *
 * Inputs: float *b1lo         low limit of left background region
 *         float *b1hi         high limit of left background region
 *         float *b2lo         low limit of right background region
 *         float *b2hi         high limit of right background region
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int get_bkgd(float *b1lo, float *b1hi, float *b2lo, float *b2hi)
{
  char line[MAXC];   /* General string for getting input */

  return 0;
}

/*.......................................................................
 *
 * Function find_wmean
 *
 * Calculates the weighted mean of the background-subtracted profile.
 *
 * Inputs: Pos *slice          uncorrected profile
 *         int npoints         number of points in the profile
 *         float m,b           parameters describing background
 *         int ysmax           location of profile maximum
 *         float *ymean        weighted mean (set by this function)
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int find_wmean(Pos *slice, int npoints, float m, float b, int ysmax,
	       float *ymean)
{
  int i;           /* Looping variable */
  float count=0.0; /* Running count */
  float sum=0.0;   /* Running sum */
  Pos *sptr;   /* Pointer for navigating slice */

  for(i=0,sptr=slice; i<npoints; i++,sptr++) {
    sum += ((sptr->y - (m*sptr->x + b)) * sptr->x);
    count += (sptr->y - (m*sptr->x + b));
  }

  *ymean = sum / count;
  printf("ymean = %f, ysmax = %d\n",*ymean,ysmax);

  return 0;
}


int MAIN_(void)
{
  return 0;
}
