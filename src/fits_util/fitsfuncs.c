/*.......................................................................
 *
 * fitsfuncs.c
 *
 * A series of generic functions dealing with fits files.  These heavily
 *  use Martin Shepherd's fits libraries.
 *
 * Functions included in this library:
 * -----------------------------------
 *  new_Image
 *  del_Image
 *  load_image_data
 *  image_info
 *  print_image_info
 *  imstat
 *  im_rms
 *  im_median
 *  sigclip
 *  floatcmp
 *
 *
 *-----------------------------------------------------------------------
 * Revision history:
 * -----------------
 *  2008Jul10 - Chris Fassnacht (CDF) - First split off from fitsim.c. in
 *    the fitsplt directory
 *   Split new_Image into an informational part (still called new_Image) and
 *    the part that actually loads the data (called load_image_data)
 *  2008Jul14 - CDF - Moved the imstat, im_rms, im_median, sigclip, and 
 *    floatcmp functions from fitsim.c to this library.
 *  2008Jul15 - CDF - Moved the calculation of the image mean to load_image_data
 *    function from the imstat function.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "structdef.h"
#include "coords.h"
#include "libfits.h"
#include "fitsfuncs.h"

/*.......................................................................
 *
 * Function new_Image
 *
 * Allocate memory for and make a copy of the 2D image contained in a
 *  FITS file primary HDU.
 *
 * Input:  char *file_name     name of the fits file

 * Output: Image *image        image data, or NULL on error
 *
 */

Image *new_Image(char *file_name)
{
  Image *image;     /* The dynamically allocated output container */
  Phdu *phdu;       /* Primary header-data-unit context descriptor */

  /*
   * Allocate the return container.
   */

  image = (Image *) malloc(sizeof(Image));
  if(!image) {
    fprintf(stderr, "new_Image: Insufficient memory.\n");
    return NULL;
  };

  /*
   * Before attempting any operation that might fail, initialize
   * the members of the container at least up to the point at which
   * it can be safely passed to del_Image().
   */

  image->fits = NULL;
  image->imaxis = NULL;
  sprintf(image->telescope,"N/A");
  sprintf(image->instrument,"N/A");
  image->badaxis = 0;
  image->nx = 0;
  image->ny = 0;
  image->data = NULL;
  image->datamin = image->datamax = 0.0f;
  image->cdmatx = 0;
  image->haswcsinfo = 0;
  image->pixscale = 0.0;
  image->rotation = 0.0;
  image->pixmean = 0.0;
  image->pixmedian = 0.0;
  image->pixmode = 0.0;
  image->pixrms = 0.0;

  /*
   * Attempt to open the file.
   */

  printf("\nReading image from %s.....\n\n",file_name);
  image->fits = new_Fits(file_name, 1, 1, 1, 0);
  if(!image->fits) {
    fprintf(stderr,"\n\nERROR: new_Image.  Could not read fits file %s\n\n",
	    file_name);
    return del_Image(image);
  }
  sprintf(image->filename,"%s",file_name);

  /*
   * Get a pointer to the primary header-data-unit context descriptor.
   */

  phdu = (Phdu *) image->fits->hdu;

  /*
   * We can only cope with 2D arrays.
   */

  if(phdu->naxis < 2) {
    fprintf(stderr, "new_Image: Can't handle %dD image arrays.\n",
	    phdu->naxis);
    return del_Image(image);
  }

  /*
   * Record the dimensions of the array.
   */

  image->nx = phdu->dims[0];
  image->ny = phdu->dims[1];
  image->ntotal = image->nx * image->ny;

  /*
   * Sanity check the dimensions.
   */

  if(image->nx < 1 || image->ny < 1) {
    fprintf(stderr, "new_Image: Array size <= 0.\n");
    return del_Image(image);
  };

  /*
   * Print out some basic image information
   */

  if(print_image_info(image) == ERROR)
    return del_Image(image);
  else
    return image;
}

/*.......................................................................
 *
 * Function del_Image
 *
 * Delete an Image context descriptor.
 *
 * Input:  Image *image        image to be freed
 *
 * Output: NULL
 *
 */

Image *del_Image(Image *image)
{
  if(image) {
    image->fits = del_Fits(image->fits);
    image->imaxis = del_Imaxis(image->imaxis);
    if(image->data)
      free(image->data);
  };

  return NULL;
}

/*.......................................................................
 *
 * Function load_image_data
 *
 * Loads the data from a fits file into a previously-defined Image container
 *
 * Inputs:
 *  Image *image               container for the data
 *
 * Output:
 *  int (SUCCESS or ERROR)
 *
 */

int load_image_data(Image *image)
{
  int i;        /* Looping variable */
  Phdu *phdu;   /* Pointer to primary header-data unit */

  /*
   * Allocate the image array.
   */

  if(!(image->data = new_array(image->nx,image->ny))) {
    fprintf(stderr,
	    "load_image_data: Insufficient memory for %dx%d image array.\n",
	    image->nx, image->ny);
    return ERROR;
  };

  /*
   * Read data into the image array.
   */

  phdu = (Phdu *) image->fits->hdu;
  if(rimage(image->fits, phdu, 0, 0, image->ntotal, DAT_FLT, 1,
	    NULL, image->data) < image->ntotal) {
    fprintf(stderr, "load_image_data: Unable to read Fits image array.\n");
    return ERROR;
  };

  /*
   * Determine the range of data values in and mean of the image array.
   */

  image->datamin = image->datamax = image->data[0];
  image->pixmean = 0.0;
  for(i=0; i<image->ntotal; i++) {
    float value = image->data[i];
    image->pixmean += value;
    if(value < image->datamin)
      image->datamin = value;
    else if(value > image->datamax)
      image->datamax = value;
  };
  image->pixmean /= image->ntotal;

  return SUCCESS;
}

/*....................................................................... 
 *
 * Function image_info
 *
 * Gets general image information from image header and prints it out.
 *
 *
 * Inputs: 
 *  Fits *fits          image container
 *
 * Output: 
 *  int (SUCCESS or ERROR)
 *
 * Revision history
 * ----------------
 *  2003Aug26 - CDF - Moved acquisition of axis information into the new
 *   axis_info function.
 *  2008Jul10 - CDF - Moved this function from fitsim.c to fitsfuncs.c
 *   and changed passed parameter from an Image structure to a Fits structure
 */

int image_info(Fits *fits) 
{
#if 0
  int i;            /* Looping variable */
  Phdu *phdu;       /* Primary header-data-unit context descriptor */
  Imaxis *axptr;    /* Pointer to navigate image->imaxis */
  
  /*
   * Get a pointer to the primary header-data-unit context descriptor.
   */

  phdu = (Phdu *) fits->hdu;

  /*
   * Load telescope and instrument info
   */

  if(phdu->telescop)
    strcpy(image->telescope,phdu->telescop);

  if(phdu->instrume)
    strcpy(image->instrument,phdu->instrume);

  /*
   * Get pixel scale from image header
   */

  if(read_file_pixscale(image,3600.0)) {
    fprintf(stderr,"WARNING: image_info. Could not get pixel scale.***\n");
  }

  /*
   * Get axis information
   */

  if(axis_info(image)) {
    fprintf(stderr,"ERROR: image_info\n");
    return 1;
  }

  /*
   * Print out information about the image
   */

  print_image_info(image);
#endif
  return 0;
}

/*.......................................................................
 *
 * Function print_image_info
 *
 * Prints out information obtained from running image_info.
 *
 * Inputs: Image *image          fits container
 *
 * Output:
 *  int (SUCCESS or ERROR)
 *
 */

int print_image_info(Image *image)
{
  int i;               /* Looping variable */
  Phdu *phdu;          /* Primary header-data-unit context descriptor */
  WCSinfo wcs;         /* Container for WCS information */
  Imaxis *axes=NULL;   /* Container for axis information */
  Imaxis *axptr;       /* Pointer to navigate image->imaxis */

  /*
   * Get a pointer to the primary header-data-unit context descriptor.
   */

  phdu = (Phdu *) image->fits->hdu;

  /*
   * Get image pixel scale
   */

  if(read_wcs_info(image->fits,&image->wcs)) {
    fprintf(stderr,"WARNING: print_fits_info: Could not get pixel scale.\n");
    /* return ERROR; */
  }
  else if(image->wcs.validwcs) {
    print_wcs_info(image->wcs);
    printf("\n");
    print_rscale(image->wcs);
  }

  /*
   * Print out the info.
   */

  printf("\n------------------------------------------------------------\n\n");
  printf("Image Information:\n");
  if(phdu->telescop) {
    printf("  Telescope:     %s\n",phdu->telescop);
  }
  if(phdu->instrume)
    printf("  Instrument:    %s\n",phdu->instrume);
  if(phdu->object)
    printf("  Object:        %s\n",phdu->object);
  if(phdu->date_obs)
    printf("  Date Observed: %s\n",phdu->date_obs);
  printf("  Image size:    %d x %d\n",phdu->dims[0],phdu->dims[1]);
  printf("  Pixel scale: %8.4f x %8.4f arcsec/pix\n",
	 image->wcs.cdelt[0]*-3600.0,image->wcs.cdelt[1]*3600);
  printf("  Data axes:\n");
  for(i=0,axptr=image->wcs.axisinfo; i<phdu->naxis; i++,axptr++) {
    if(axptr->ctype) {
      if(strncmp(axptr->ctype,"RA",2) == 0) {
	printf("    %d = %-8s:  CRPIX = %g, CRVAL =  %s, ",
	       i+1,axptr->ctype,axptr->crpix,image->wcs.fitsra);
	printf("CDELT = %g\n",axptr->cdelt);
      }
      else if(strncmp(axptr->ctype,"DEC",3) == 0) {
	printf("    %d = %-8s:  CRPIX = %g, CRVAL = %s, ",
	       i+1,axptr->ctype,axptr->crpix,image->wcs.fitsdec);
	printf("CDELT = %g\n",axptr->cdelt);
      }
      else {
	if(axptr->crpix) {
	  printf("    %d = %-8s:  CRPIX = %g, CRVAL = %g, CDELT = %g\n",
		 i+1,axptr->ctype,axptr->crpix,axptr->crval,axptr->cdelt);
	}
	else {
	  printf("No axis information on axis %d\n",i+1);
	}
      }
    }
    else
      printf("   NULL\n");
  }

  if(phdu->naxis > 2) {
    fprintf(stderr,"    ***** WARNING! Using only first two axes for plot.");
    fprintf(stderr," *****\n");
  }
  printf("\n------------------------------------------------------------\n\n");

  return SUCCESS;
}

/*.......................................................................
 *
 * Function imstat
 *
 * Calculates image statistics, such as the median pixel value, the 
 *  RMS, and estimate of the mode.
 * NB: image->pixmean is now calculated in load_image_data
 *
 * Inputs: Image *image        image container
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int imstat(Image *image)
{
  int i;                /* Looping variable */
  int no_error=1;       /* Flag set to 0 on error */

  printf("Calculating image statistics.....\n");
  printf("---------------------------------\n");
  printf("  Minimum:    %9.3f\n",image->datamin);
  printf("  Maximum:    %9.3f\n",image->datamax);
  printf("  Mean:       %9.3f\n",image->pixmean);

  /*
   * Check that there are valid pixels in the image.
   */

  if(image->ntotal == 0)  {
    fprintf(stderr,"---------------------------------\n");
    fprintf(stderr,"ERROR: imstat.\n");
    return 1;
  }

  /*
   * Estimate the RMS of the data.
   */

  if((image->pixrms = 
      im_rms(image->data,image->ntotal,image->pixmean)) < 0.0)
    no_error = 0;
  else
    printf("  RMS:        %9.3f\n",image->pixrms);

  /*
   * Find the image median
   */
#if 0
  if(no_error)
    if(im_median(image->data,image->ntotal,&image->pixmedian))
      no_error = 0;
    else
      printf("  Median:     %9.3f\n",image->pixmedian);

  /*
   * Estimate the mode following the Kendall and Stuart (1977)
   *
   *   mode = 3 x median - 2 x mean
   *
   * NB: If we ever do iterative 3-sigma clipping, then we should
   *     consider using the SExtractor formulation for estimating the
   *     mode in a clipped histogram:
   *
   *     mode = 2.5 x median - 1.5 x mean
   */

  image->pixmode = 3 * image->pixmedian - 2 * image->pixmean;
  printf("  Mode:       %9.3f\n",image->pixmode);
#endif
  /*
   * Clean up and return
   */

  printf("---------------------------------\n");
  if(no_error)
    return 0;
  else {
    fprintf(stderr,"ERROR: imstat\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function im_rms
 *
 * Estimates the rms scatter about the mean pixel value.  The function
 *  uses the unbiassed estimator of the variance:
 *
 *
 *             1
 *   s^2 = --------- SUM [(x_i - mu_x)^2]
 *           N - 1
 *
 * Inputs: float *data         array of pixel values
 *         int ndata           number of data points
 *         double mean         mean value of array
 *
 * Output: double rms          estimated rms about the mean.  NB: set to
 *                              -1.0 on error.
 *
 */

double im_rms(float *data, int ndata, double mean)
{
  int i;                /* Looping variable */
  float *dptr;          /* Pointer to navigate data */
  double sum=0.0;       /* Sum of pixel values */
  double diff;          /* (x_i - mu_x) */
  double rms;           /* Estimated rms of the distribution */

  sum = 0.0;
  for(i=0,dptr=data; i < ndata; i++,dptr++) {
    diff = *dptr - mean;
    sum += (diff * diff);
  }
  rms = sqrt(sum / (ndata - 1));

  return rms;
}

/*.......................................................................
 *
 * Function im_median
 *
 * Calculates the image median by sorting the data array (with qsort)
 *  and then taking the middle value.  Note that qsort is slow so this
 *  takes a while with large arrays.
 *
 * Inputs: float *data         array of pixel values
 *         int ndata           number of data points
 *         double *median      array median (set by this function)
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error.
 *
 */

int im_median(float *data, int ndata, double *median)
{
  int i;                /* Looping variable */
  int medindex;         /* Array position of median value */
  float *tmpdata=NULL;  /* Temporary twin to data array */
  float *dptr;          /* Pointer to navigate data */
  float *tptr;          /* Pointer to navigate tmpdata */

  /*
   * Allocate memory for temporary twin of data array.
   */

  if(!(tmpdata = new_array(ndata,1))) {
    fprintf(stderr,"ERROR: im_median.\n");
    return 1;
  }

  /*
   * Transfer data to temporary array
   */

  for(i=0,dptr=data,tptr=tmpdata; i < ndata; i++,dptr++,tptr++) {
    *tptr = *dptr;
  }

  /*
   * Sort the temporary array.
   */

  qsort(tmpdata,ndata,sizeof(tmpdata[0]),floatcmp);

  /*
   * Find median by taking middle value of array.  The definition of
   *  "middle value" depends on whether the number of elements in
   *  the array (image->ntotal) is even or odd.
   */

  if((ndata/2.0) - (int) (ndata/2.0) == 0) {
    medindex = ndata / 2;
    *median = (tmpdata[medindex-1] + tmpdata[medindex])/2.0;
  }
  else {
    medindex = (ndata - 1) / 2;
    *median = tmpdata[medindex];
  }

  /*
   * Clean up and return
   */

  tmpdata = del_array(tmpdata);

  return 0;
}

/*.......................................................................
 *
 * Function sigclip
 *
 * Function to clip the data in an array at the [nsig]-sigma level
 *  until convergence is achieved.  nsig is one of the passed parameters.
 *  Once convergence is achieved, return the clipped mean and rms values.
 *
 * Inputs: float *data         array of pixel values
 *         int ndata           number of data points
 *         double *mean        mean of full data set (re-set by this function)
 *         double *rms         rms of full data set (re-set by this function)
 *         float nsig          sigma-clipping level
 *         float ftol          minimum fractional change in rms to 
 *                              continue clipping process.
 *
 * Output: int ngood           number of good points after clipping.
 *                              Set to -1 on error.
 */

int sigclip(float *data, int ndata, double *mean, double *rms, float nsig,
	    float ftol)
{
  int i;                /* Looping variable */
  int ngood=0;          /* Number of good (unclipped) points in array */
  float *tmpdata=NULL;  /* Temporary storage for clipped data array */
  float *dptr,*tptr;    /* Pointers to navigate image->data and tmpdata */
  double cliprms;       /* RMS about mean in clipped array */
  double lastrms;       /* Holds previous version of clipped rms */
  double clipsum;       /* Sum of pixel values in clipped array */
  double clipmean;      /* Mean of clipped array */

  /*
   * Allocate memory for tmpdata
   */

  if(!(tmpdata = new_array(ndata,1))) {
    fprintf(stderr,"ERROR: sigclip\n");
    return -1;
  }

  /*
   * Initialize variables for clipping loop.
   */

  clipmean = *mean;
  cliprms = *rms;
  lastrms = 1000*cliprms;

  /*
   * Loop on clipping
   */

  printf("\nsigclip: Clipping at %5.2f sigma.\n",nsig);
  while((fabs(lastrms - cliprms)/lastrms) > ftol) {
    lastrms = cliprms;
    ngood = 0;
    clipsum = 0.0;
    tptr=tmpdata;
    for(i=0,dptr=data; i<ndata; i++,dptr++) {
      if(fabs(*dptr - clipmean) < nsig * cliprms) {
	*tptr = *dptr;
	clipsum += *tptr;
	tptr++;
	ngood++;
      }
    }
    if(ngood > 0) {
      clipmean = clipsum / ngood;
      cliprms = im_rms(tmpdata,ngood,clipmean);
      printf("sigclip: Mean = %9.3f, RMS = %9.3f ngood = %d\n",
	     clipmean,cliprms,ngood);
    }
    else {
      fprintf(stderr,"ERROR: sigclip.  No valid points in clipped array.\n");
      ngood = -1;
      break;
    }
  }

  /*
   * Set the array mean and rms to the clipped values.
   */

  *mean = clipmean;
  *rms = cliprms;

  return ngood;
}

/*.......................................................................
 *
 * Function floatcmp
 *
 * Compares the two float values and returns 1 if the first is greater 
 *  than the second, 0 if they're equal, and -1 if the first is less 
 *  than the second.  This function is called by qsort.
 *
 * Inputs: void *v1            first structure (cast to void)
 *         void *v2            second structure (cast to void)
 *
 * Output: int (-1,0,1)        as described above
 *
 */

int floatcmp(const void *v1, const void *v2)
{
  float *f1 = (float *) v1;  /* float casting of v1 */
  float *f2 = (float *) v2;  /* float casting of v2 */

  /*
   * Do the comparison
   */

  if(*f1 > *f2)
    return 1;
  else if(*f1 == *f2)
    return 0;
  else
    return -1;
}
