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
#include "libfits.h"
#include "cpgplot.h"
#include "structdef.h"
#include "plotfuncs.h"
#include "fitsim.h"

/*.......................................................................
 *
 * Main program
 *
 */

int main(int argc, char *argv[])
{
  int i;
  int no_error=1;     /* Flag set to 0 on error */
  int ymax;           /* Y value associated with slmax */
  float immean;       /* Mean value of image */
  float imrms;        /* RMS of image */
  float slmax;        /* Max value along a slice */
  float sky=0.0;      /* Temporary variable */
  float m,b;          /* Parameters for linear fit to bkgd */
  float *x=NULL;
  float *y=NULL;
  char *file_name;    /* The name of the FITS file to load */
  char *setupname;    /* Name of the optional setup file */
  Image *image=NULL;  /* The container of the FITS image array */
  Setup *setup=NULL;  /* Container for the image display info */
  Pos *slice=NULL;    /* A 1D slice through the image */

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
   * Calculate image statistics
   */

  if(calc_imstat(image,&immean,&imrms))
    no_error = 0;

  /*
   *  Allocate memory for the slice array
   */

  if(!(slice = new_pos(image->ny,1)))
    no_error = 0;
  if(!(x = new_array(image->ny,1)))
    no_error = 0;
  if(!(y = new_array(image->ny,1)))
    no_error = 0;

  /*
   * Plot a vertical cut through the image at the central column
   *  and find its maximum value.
   */

  plot_open_spec(1.5);
  if(no_error) {
    if(plot_vslice(image,image->nx/2,0,image->ny-1))
      no_error = 0;
    else
      find_vslice_max(image,image->nx/2,0,image->ny-1,&slmax,&ymax);
  }

  /*
   * Replot the central column, expanded around the max
   */

  if(no_error) {
    cpgpage();
    if(plot_vslice(image,image->nx/2,ymax-50,ymax+50))
       no_error = 0;
  }

  /*
   * Fit a linear continuum
   */

  if(no_error) {
    m = b = 0.0;
    if(fit_continuum(image,image->nx/2,ymax,-30,-15,15,30,slice,
		     &m,&b,&sky))
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
  }

  /*
   * Slice up the image
   */

  if(no_error)
    if(imslice(image))
      no_error = 0;

  /*
   * Clean up and exit
   */

  image = del_Image(image);
  setup = del_setup(setup);
  slice = del_pos(slice);
  x = del_array(x);
  y = del_array(y);
  plot_close();

  if(no_error)
    return 0;
  else {
    fprintf(stderr,"ERROR. Exiting program\n\n");
    return 1;
  }
}


int MAIN_(void)
{
  return 0;
}
