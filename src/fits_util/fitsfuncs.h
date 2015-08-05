#ifndef fitsfuncs_h
#define fitsfuncs_h

#include "structdef.h"
#include "coords.h"
#include "libfits.h"
#include "fitswcs.h"

/*.......................................................................
 *
 * Structure definitions
 *
 */

typedef struct {
  Fits *fits;          /* The FITS file context descriptor */
  Imaxis *imaxis;      /* The data axes */
  char filename[80];   /* The filename of the fits file */
  int badaxis;         /* Flag set to 1 if problems reading an axis */
  char telescope[80];  /* Telescope used for the observation */
  char instrument[80]; /* Instrument used for the observation */
  int nx,ny;           /* The dimensions of data[] */
  int ntotal;          /* nx * ny */
  float *data;         /* 2D data array read from FITS file (FORTRAN order) */
  float datamin;       /* The minimum value in the data array */
  float datamax;       /* The maximum value in the data array */
  double pixmean;      /* Mean pixel value */
  double pixmedian;    /* Median pixel value */
  double pixmode;      /* Estimate of pixel mode */
  double pixrms;       /* Estimate of pixel rms */
  WCSinfo wcs;         /* WCS information (CD matrix, CRVAL, etc) */
  int haswcsinfo;      /* Flag set to 1 if the file has valid WCS info */
  int cdmatx;          /* Flag set to 1 if the file has a CD matrix */
  float pixscale;      /* Pixel scale in the image */
  double rotation;     /* WCS axis rotation (N->E) */
  char fitsra[80];     /* RA of crval in hh:mm:ss.ss format */
  char fitsdec[80];    /* Dec of crval in dd:mm:ss.ss format */
} Image;

/*.......................................................................
 *
 * Function declarations
 *
 */

Image *new_Image(char *file_name);
Image *del_Image(Image *image);
int load_image_data(Image *image);
int image_info(Fits *fits);
int print_image_info(Image *image);
int imstat(Image *image);
double im_rms(float *data, int ndata, double mean);
int im_median(float *data, int ndata, double *median);
int sigclip(float *data, int ndata, double *mean, double *rms, float nsig,
	    float ftol);
int floatcmp(const void *v1, const void *v2);

#endif
