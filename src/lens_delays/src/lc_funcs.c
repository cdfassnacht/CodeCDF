/*
 * lc_funcs.c
 *
 * A library of functions to perform operations on 1608 light curves.
 *
 * 18Jun98 CDF,  Moved functions over from lcurve.c.  For history prior
 *                to this date, see lcurve.c
 *               Added functions read_fluxrec and flag_bad.  Modified
 *                many other functions to deal with flagged points on curves
 * v19Jun98 CDF, Make smooth_1608 return size of smoothed curves.
 *               Add array initialization to new_fluxrec.
 *               Add logfile option to do_chi and call_xcorr.
 *               Move optional interpolation from shift_chi to do_chi and
 *                make shift_chi return the best-fit value.
 * v22Jun98 CDF, Add proper normalization to cross_corr_fft.
 *               Added rand_curves to create light curves with the same
 *                flux distributions as the observed light curves, but with
 *                the time series randomized.
 *               Added the triangle function to do a triangle smooth on the
 *                light curves.
 * v23Jun98 CDF, Added setup functions setup_file, setup_interact, and
 *                read_setup_line, as well as new_setup and del_setup.
 *                These function allow input from setup files similar to
 *                those used in fitsplt.c/fitsim.c
 * v24Jun98 CDF, Added cross_corr_nr, which calls the Numerical Recipes
 *                version of cross-correlation.  This can provide a
 *                check for the Fourier transform method in cross_corr_fft
 *               Added norm_zero_mean, which normalizes a light curve and
 *                then subtracts 1 to create a zero-mean light curve.
 * v25Jun98 CDF, Add calculation of rms in cross-correlation curves to assess
 *                the significance of the peaks.
 *               Slight modifications to calc_mean and do_chi.
 *               Added varbox function to do variable-width boxcar smoothing.
 * v28Jun98 CDF, Moved setup_summary from lcurve.
 * v29Jun98 CDF, Added gaussian to do gaussian-weighted smoothing.  
 *               Modified setup_file and setup_interact to deal with 
 *                gaussian weighting.
 *               Changed uncertainties on boxcar-smoothed values to be
 *                max of weighted sum of variances and rms scatter about
 *                mean.
 * v04Jul98 CDF, Modification of output in call_xcorr
 *               Added make_monte to create fake light curves
 *                with sparse sampling from an idealized, regularly-
 *                sampled light curve.
 * v08Jul98 CDF, Changed output flux ratios in do_chi from B/[ACD] to
 *                [ACD]/B.
 *               Added the gaussian random displacements to the idealized
 *                light curve in make_monte.  Also delete bad days from
 *                output.
 *               Added interactive file request to read_fluxrec.
 * v09Jul98 CDF, Moved choice of linear interpolation out of do_chi and up
 *                into the calling function.  As a consequence, added a
 *                dosmooth flag to the Setup structure.
 *               Split logfile option into two -- a logfile for the chisq
 *                portion and a logfile for the cross-correlation portion.
 * v12Jul98 CDF, Added functions medsmooth and vartri to do median smoothing
 *                and variable-width triangle smoothing, respectively.
 *                Modified the setup container functions appropriately.
 * v13Jul98 CDF, Fixed a bug in setup_interact in which defaults weren't being
 *                handled properly.
 *               Added cross_corr_time to do a "by hand" cross-correlation
 *                without using FFTs.
 *               Added better output to do_chi.
 * v14Jul98 CDF, Add variance-weighting to all smoothing functions by adding
 *                the var_wmean function and having all the smoothing functions
 *                call it.
 * v16Jul98 CDF, Changed zero_pad to produce float arrays rather than
 *                Fluxrec arrays.  Modified do_corr, cross_corr_fft, 
 *                cross_corr_nr, and cross_corr_time accordingly.
 *               Added wt_hanning function to apply Hanning weighting to the
 *                zero-mean curves that are input to the cross_corr_fft function.
 * v19Jul98 CDF, Took out Hanning weighting in do_corr.
 *               Fixed bug in cross_corr_time.
 *               Made steps in flux density ratio into fractions of the initial
 *                guess rather than absolute steps (shift_chi).
 * v20Jul98 CDF, Added flexibility in the number of zeros used in zero_pad
 *                by putting in the ZEROFAC factor.
 * v21Jul98 CDF, Fixed bug in smooth_1608 in which the errors on the
 *                flat-fielded curves were not being calculated correctly.
 *                The function now calls flat_field, which correctly incorporates
 *                the fractional rms contribution to the flux density errors.
 *               Better logfile output from do_chi.
 * v22Jul98 CDF, Modified smooth_1608 to make printing to the output file
 *                optional. 
 * v04Aug98 CDF, Commented out "reverse" case in do_chi
 * v13Aug98 CDF, Added call_dcf and discrete_corr to calculate the
 *                discrete correlation function for the light curves.
 * v15Aug98 CDF, Added functions fit_poly and poly to fit a polynomial
 *                to a light curve.
 * v16Aug98 CDF, Added function chisq_fit to calculate a reduced chisq between
 *                a light curve and the polynomial fit.
 * v20Aug98 CDF, Added print_log to eliminate lines needed for adding to
 *                log files.
 *               Added a choice of smoothing onto a coarse grid and then
 *                interpolating onto a finer one.
 * v21Aug98 CDF, Modifications of flag_bad, interp_curve and smooth_1608.
 * v22Aug98 CDF, Completely re-wrote cross_corr_time, and as a result,
 *                split call_xcorr into call_xcorr_fft and call_xcorr_time.
 *               Added fit_parab function from modfuncs.c
 * v23Aug98 CDF, Modified shift_chi to fit a parabola to the lag-chisq curve
 *                at the location of the minimum value of chisq found on the
 *                gridded curve.  This function, performed in find_min_chisq,
 *                allows the determination of the "true" value of the minimum
 *                chisq.
 * v25Aug98 CDF, Changed rand_curves to only output 1608 curves since
 *                mkrand.c now flat-fields the curves before passing them
 *                to rand_curves.
 * v30Aug98 CDF, Re-wrote do_chi and shift_chi to handle input curves of
 *                different lengths and spacings.
 *               Moved flat-fielding from interp_curve to calling programs.
 *               Moved initialization of smoothed curves out of the smoothing
 *                functions (boxcar, triangle, etc.) and into csmooth, which
 *                calls them.  This is to prepare for being able to call the
 *                smoothing functions from a "smooth-in-place" function
 *                which will not output a regularly sampled grid.
 * v31Aug98 CDF, Modified smooth_1608 to handle in-place smoothing.
 * v01Sep98 CDF, Modified flag_bad to have an option to flag individual curves
 *                in addition to flagging all fluxes for a given day.
 * v07Sep98 CDF, Added corr_bevington to calculate correlation coefficients
 *                using the method described in Bevington and bevington_prob
 *                to calculate the probability that the value of a correlation
 *                coefficient could arise from two uncorrelated curves.
 * v11Sep98 CDF, Added gauss_noise to generate a time series of zero-mean
 *                Gaussian-distributed random noise.
 * v12Sep98 CDF, Modified cross_corr_time to estimate the number of 
 *                independent points in the overlap region rather than
 *                just taking all the points in the overlap region to
 *                calculate the Bevington probability.
 * v24Sep98 CDF, Modified shift_chi in a similar way to change the
 *                calculation of the reduced chisq.
 * v25Sep98 CDF, Split off smoothing/interpolation functions to lc_interp.c
 *               Split off correlation functions to correlate.c
 * v01Oct98 CDF, Split off Monte Carlo related functions to monte.c
 * v02Dec98 CDF, Added make_compos function to make a composite light
 *                curve from two to four input light curves.
 * v03Dec98 CDF, Modified setup functions to include choice of analysis
 *                method in the setup container.
 * v25Mar99 CDF, Changed fractional rms returned by ratio_err to be a factor
 *                of sqrt(2) less, after discussions with LVK.
 * v13Jun99 CDF, Moved Setup structure handing into new lc_setup.c library.
 *               Added choice of secondary flux curve calibration to make_flat.
 * v27Jul00 CDF, Moved all chisq calculation functions into new lc_chisq.c
 *                library.
 * v30Jul00 CDF, Added the set_mu0 function to set the initial guesses for
 *                the component relative magnifications.
 * v04Oct00 CDF, Added a new write_1608 function that writes out 1608
 *                light curves to a file in a format that can be read by the
 *                read_1608 function.
 * v09Oct00 CDF, Added new_prange and del_prange to dynamically allocate
 *                memory for Prange arrays.
 * v03Apr02 CDF, Added a new write_fluxrec function.
 * v02Sep02 CDF, Improved file I/O in read_fluxrec function.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "structdef.h"
#include "dataio.h"
#include "nrutil.h"
#include "nr.h"
#include "lc_setup.h"
#include "lc_funcs.h"

/*.......................................................................
 *
 * Function new_prange
 *
 * Allocates dynamic memory for a pointer array of Prange structures
 *
 * Input:  int size            size of the array
 *
 * Output: Prange *newinfo     pointer to the new array.  NULL if error
 *
 */

Prange *new_prange(int size)
{
  int i;
  Prange *newinfo;
  Prange *fptr;
 
  newinfo = (Prange *) malloc(sizeof(Prange) * size);
  if(!newinfo) {
    fprintf(stderr,"Insufficient memory for data array.\n");
    return NULL;
  }

  /*
   * Initialize values in the array
   */

  for(i=0,fptr=newinfo; i<size; i++,fptr++) {
    fptr->val0 = fptr->dval = 0.0;
    fptr->nval = fptr->minstep = fptr->maxstep = 0;
  }

  return newinfo;
}

/*.......................................................................
 *
 * Function del_prange
 *
 * Frees up memory allocated to Prange array
 *
 * Input:  Prange *prange      array to be freed
 *
 * Output: NULL
 */

Prange *del_prange(Prange *prange)
{
  if(prange)
    free(prange);
 
  return NULL;
}

/*.......................................................................
 *
 * Function new_fluxrec
 *
 * Allocates dynamic memory for a pointer array of Fluxrec structures
 *
 * Input:  int size            size of the array
 *
 * Output: Fluxrec *newinfo    pointer to the new array.  NULL if error
 *
 */

Fluxrec *new_fluxrec(int size)
{
  int i;
  Fluxrec *newinfo;
  Fluxrec *fptr;
 
  newinfo = (Fluxrec *) malloc(sizeof(Fluxrec) * size);
  if(!newinfo) {
    fprintf(stderr,"Insufficient memory for data array.\n");
    return NULL;
  }

  /*
   * Initialize values in the array
   */

  for(i=0,fptr=newinfo; i<size; i++,fptr++) {
    fptr->day = fptr->flux = fptr->err = 0.0;
    fptr->match = 0;
  }

  return newinfo;
}

/*.......................................................................
 *
 * Function del_fluxrec
 *
 * Frees up memory allocated to Fluxrec array
 *
 * Input:  Fluxrec *fluxrec    array to be freed
 *
 * Output: NULL
 */

Fluxrec *del_fluxrec(Fluxrec *fluxrec)
{
  if(fluxrec)
    free(fluxrec);
 
  return NULL;
}

/*.......................................................................
 * 
 * Function print_log
 *
 * Prints a string into the chisq and cross-correlation logfiles, if
 *  they exist.
 *
 * Inputs: FILE *chifp         chisq file pointer
 *         FILE *xcfp          cross-correlation file pointer
 *         char *logstring     string to be printed
 *
 * Output: none
 *
 */

void print_log(FILE *chifp, FILE *xcfp, char *logstring)
{
  if(chifp)
    fprintf(chifp,logstring);

  if(xcfp)
    fprintf(xcfp,logstring);
}

/*.......................................................................
 *
 * Function read_data
 *
 * Reads data from input file and fills Fluxrec arrays.
 *
 * Inputs: Fluxrec *fl34       1634 lightcurve
 *         Fluxrec *fl35       1635 lightcurve
 *         Fluxrec *fl08[]     1608 lightcurves
 *         int nlines          number of lines in the array
 *         FILE *lfp           lens input file pointer
 *         FILE *cfp           comp input file pointer
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> failure
 *
 * v30Aug98 CDF, Added line to set the match member of the new structures to 1.
 */

int read_data(Fluxrec *fl34, Fluxrec *fl35, Fluxrec *fl08[], 
	      int nlines, FILE *lfp, FILE *cfp)
{
  int count=0;        /* Counts number of input lines read */
  float junk;         /* Variable for reading data which won't be used */
  char line[MAXC];    /* General char string for reading input */
  Fluxrec *f34p;      /* Pointer for fl34 */
  Fluxrec *f35p;      /* Pointer for fl35 */
  Fluxrec *fap;       /* Pointer for fla */
  Fluxrec *fbp;       /* Pointer for flb */
  Fluxrec *fcp;       /* Pointer for flc */
  Fluxrec *fdp;       /* Pointer for fld */

  /*
   * Initialize pointers
   */

  f34p = fl34;
  f35p = fl35;
  fap = fl08[0];
  fbp = fl08[1];
  fcp = fl08[2];
  fdp = fl08[3];

  /*
   * Read in lens data
   */

  while(fgets(line,MAXC,lfp) != NULL) {
    if(line[0] != '#' && count < nlines) {
      if(sscanf(line,"%f %f %f %f %f %f %f %f %f %f",&fap->day,
		&fap->flux,&fbp->flux,&fcp->flux,&fdp->flux,
		&junk,&junk,&junk,&junk,&fap->err) != 10) {
	fprintf(stderr,"ERROR: read_data.  Input file has incorrect format\n");
	return 1;
      }
      fbp->day = fcp->day = fdp->day = fap->day;
      fbp->err = fcp->err = fdp->err = fap->err;
      fap->match = fbp->match = fcp->match = fdp->match = 1;
      count++;
      fap++;
      fbp++;
      fcp++;
      fdp++;
    }
  }

  /*
   * Read in comparison source data
   */

  count = 0;
  while(fgets(line,MAXC,cfp) != NULL) {
    if(line[0] != '#' && count < nlines) {
      if(sscanf(line,"%f %f %f %f %f %f",&f34p->day,
		&f34p->flux,&f34p->err,&junk,&f35p->flux,&f35p->err) != 6) {
	fprintf(stderr,"ERROR: read_data.  Input file has incorrect format\n");
	return 1;
      }
      f35p->day = f34p->day;
      count++;
      f34p++;
      f35p++;
    }
  }

  /*
   * Error checking and return
   */

  if(count != nlines) {
    fprintf(stderr,"ERROR: read_data. Read %d lines -- expected %d\n",
	    count,nlines);
    return 1;
  }

  return 0;
}

/*.......................................................................,
 *
 * Function read_1608
 *
 * Reads a 1608 light curve input file and puts results into a previously
 *  allocated Fluxrec array.  The input file contains the light
 *  curves and their errors in the following format:
 *
 *    day S_A S_B S_C S_D err_A err_B err_C err_D
 *
 * Inputs: Fluxrec *fl08[]     1608 lightcurves
 *         int nlines          number of lines in the array
 *         FILE *lfp           1608 input file pointer
 *         char *filename      name of input file
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> failure
 *
 */

int read_1608(Fluxrec *fl08[], int nlines, FILE *lfp, char *filename)
{
  int count=0;        /* Counts number of input lines read */
  char line[MAXC];    /* General char string for reading input */
  Fluxrec *fap;       /* Pointer for fla */
  Fluxrec *fbp;       /* Pointer for flb */
  Fluxrec *fcp;       /* Pointer for flc */
  Fluxrec *fdp;       /* Pointer for fld */

  /*
   * Initialize pointers
   */

  fap = fl08[0];
  fbp = fl08[1];
  fcp = fl08[2];
  fdp = fl08[3];

  /*
   * Read in data
   */

  while(fgets(line,MAXC,lfp) != NULL) {
    if(line[0] != '#' && count < nlines) {
      if(sscanf(line,"%f %f %f %f %f %f %f %f %f",&fap->day,
		&fap->flux,&fbp->flux,&fcp->flux,&fdp->flux,
		&fap->err,&fbp->err,&fcp->err,&fdp->err) != 9) {
	fprintf(stderr,"ERROR: read_1608.  Input file has incorrect format\n");
	return 1;
      }
      fbp->day = fcp->day = fdp->day = fap->day;
      fap->match = fbp->match = fcp->match = fdp->match = 1;
      count++;
      fap++;
      fbp++;
      fcp++;
      fdp++;
    }
  }

  /*
   * Error checking and return
   */

  if(count != nlines) {
    fprintf(stderr,"ERROR: read_1608. Read %d lines -- expected %d\n",
	    count,nlines);
    return 1;
  }
  else {
    printf("read_1608: Read %d lines from %s.\n",count,filename);
  }

  return 0;
}

/*.......................................................................
 *
 * Function write_1608
 *
 * Writes out the day, four flux densities, and four errors for a 
 *  four-image lens light curve in a format that can be read by
 *  read_1608.
 *
 * Inputs: Fluxrec *flux[]     lens lightcurves
 *         int npoints         number of points in the curve
 *         char *filename      name of output file
 *         int verbose         flag set to 1 for verbose output
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> failure
 *
 */

int write_1608(Fluxrec *flux[], int npoints, char *filename, int verbose)
{
  int i;                            /* Looping variable */
  Fluxrec *ptr1,*ptr2,*ptr3,*ptr4;  /* Pointers to navigate flux */
  FILE *ofp=NULL;                   /* Output file pointer */

  /*
   * Open output file
   */

  if(!(ofp = open_writefile(filename))) {
    fprintf(stderr,"ERROR: write_1608\n");
    return 1;
  }

  /*
   * Print out light curves
   */

  for(i=0,ptr1=flux[0],ptr2=flux[1],ptr3=flux[2],ptr4=flux[3]; 
      i<npoints; i++,ptr1++,ptr2++,ptr3++,ptr4++) {
    fprintf(ofp,"%7.2f %8.5f %8.5f %8.5f %8.5f %7.5f %7.5f %7.5f %7.5f\n",
	    ptr1->day,ptr1->flux,ptr2->flux,ptr3->flux,ptr4->flux,
	    ptr1->err,ptr2->err,ptr3->err,ptr4->err);
  }

  /*
   * Clean up
   */

  if(ofp)
    fclose(ofp);

  return 0;
}


/*.......................................................................,
 *
 * Function read_fluxrec
 *
 * Reads a 2 or 3 column input file and puts results into a Fluxrec array,
 *  which has its memory allocation performed in the function.
 *
 * Inputs: char *inname        name of input file
 *         char comment        comment character
 *         int *nlines         number of lines in input file -- set by
 *                              this function.
 *         
 * Output: Fluxrec *newflux    filled array
 */

Fluxrec *read_fluxrec(char *inname, char comment, int *nlines)
{
  int no_error=1;        /* Flag set to 0 on error */
  int ncols;             /* Number of columns in input file */
  char line[MAXC];       /* General input string */
  Fluxrec *newflux=NULL; /* Filled fluxrec array */
  Fluxrec *fptr;         /* Pointer to navigate fluxrec */
  FILE *ifp=NULL;        /* Input file pointer */

  /*
   * Open output file
   */

  if(!(ifp = open_readfile(inname))) {
    fprintf(stderr,"ERROR: read_fluxrec.  Cannot open %s.\n",inname);
    return NULL;
  }

  /*
   * Get number of lines in input file
   */

  if((*nlines = n_lines(ifp,comment)) == 0) {
    fprintf(stderr,"ERROR: read_fluxrec.  No valid data in input file.\n");
    no_error = 0;
  }
  else
    rewind(ifp);

  /*
   * Allocate memory for fluxrec
   */

  if(no_error)
    if(!(newflux = new_fluxrec(*nlines)))
      no_error = 0;

  /*
   * Read in data
   */

  if(no_error) {
    fptr = newflux;
    while(fgets(line,MAXC,ifp) != NULL) {
      if(line[0] != comment && no_error) {
	if(sscanf(line,"%f %f %f %d",&fptr->day,&fptr->flux,&fptr->err,
		       &fptr->match) == 4) {
	  ncols = 4;
	  fptr++;
	}
	else if(sscanf(line,"%f %f %f",&fptr->day,&fptr->flux,&fptr->err) 
		== 3) {
	  ncols = 3;
	  fptr++;
	}
	else if(sscanf(line,"%f %f",&fptr->day,&fptr->flux) == 2) {
	  ncols = 2;
	  fptr++;
	}
	else {
	  fprintf(stderr,"ERROR: read_fluxrec. Bad input file format.\n");
	  no_error = 0;
	}
      }
    }
  }

  /*
   * Clean up and exit
   */

  if(ifp)
    fclose(ifp);

  if(no_error) {
    printf("\nread_fluxrec: %s has %d columns and %d lines\n\n",inname,
	   ncols,*nlines);
    return newflux;
  }
  else {
    fprintf(stderr,"ERROR: read_fluxrec.\n");
    return del_fluxrec(newflux);
  }
}

/*.......................................................................
 *
 * Function write_fluxrec
 *
 * Writes out a fluxrec array.  If requested, the function will check for
 *  duplicate days and only print out unique days.
 *
 * Inputs: Fluxrec *fluxrec    fluxrec array
 *         int npoints         number of points in array
 *         char *outfile       output file name
 *         int no_dup_days     flag set to 1 to reject duplicate days
 *         float dtol          minimum separation between days for those
 *                              days to be considered unique.
 *
 * Output: int 0 or 1          0 ==> success, 1 ==> error
 *
 */

int write_fluxrec(Fluxrec *fluxrec, int npoints, char *outfile,
		  int no_dup_days, float dtol)
{
  int i;           /* Looping variable */
  int ndelete=0;   /* Number of duplicate days deleted */
  Fluxrec *fptr;   /* Pointer to navigate fluxrec */
  Fluxrec *fhold;  /* Another pointer to navigate fluxrec */
  FILE *ofp=NULL;  /* Output file pointer */

  /*
   * Open output file
   */


  if(!(ofp = open_writefile(outfile))) {
    fprintf(stderr,"ERROR: write_fluxrec.\n");
    return 1;
  }

  /*
   * Initialize
   */

  fhold = fluxrec;

  /*
   * If the no_dup_days flag is set then only print out unique days.  
   */

  if(no_dup_days) {
    for(i=0,fptr=fluxrec; i<npoints; i++,fptr++) {

      /*
       * Test whether or not the day at the current point in the array is 
       *  the same as the previous day by comparing fptr->day to fhold->day.
       *
       * Of course if this is the first day (i==0) then no checking is done.
       */

      if(i>0 && (fptr->day - fhold->day) < dtol) {
	ndelete++;
      }
      else {
	fprintf(ofp,"%7.2f %10.4f %8.5f %d\n",fptr->day,fptr->flux,fptr->err,
		fptr->match);
	fhold = fptr;
      }
    }
    printf("write_fluxrec: %d duplicate days were discarded from %s.\n",
	   ndelete,outfile);
  }

  /*
   * If no_dup_days is not set, just print out all the days.
   */

  else {
    for(i=0,fptr=fluxrec; i<npoints; i++,fptr++)
      fprintf(ofp,"%7.2f %10.4f %8.5f %d\n",fptr->day,fptr->flux,fptr->err,
	      fptr->match);
  }

  /*
   * Clean up and exit
   */

  if(ofp)
    fclose(ofp);

  return 0;
}

/*.......................................................................
 *
 * Function read_bad_ext
 *
 * Reads a bad day list that is supposed to contain flags for the
 *  individual light curves.
 *
 * Inputs: char *inname        name of input file
 *         char comment        comment character
 *         Fluxrec *bada       flags for A curve (set by this function)
 *         Fluxrec *badb       flags for B curve (set by this function)
 *         Fluxrec *badc       flags for C curve (set by this function)
 *         Fluxrec *badd       flags for D curve (set by this function)
 *         int *nlines         number of lines in input file -- set by
 *                              this function.
 *
 * Output: int 0 or 1          0 ==> success, 1 ==> error
 *
 */

int read_bad_ext(char *inname, char comment, Fluxrec *bada, Fluxrec *badb,
		 Fluxrec *badc, Fluxrec *badd)
{
  int no_error=1;        /* Flag set to 0 on error */
  float junk;            /* Holds unwanted input */
  char line[MAXC];       /* General input string */
  Fluxrec *aptr;         /* Pointer to navigate bada */
  Fluxrec *bptr;         /* Pointer to navigate badb */
  Fluxrec *cptr;         /* Pointer to navigate badc */
  Fluxrec *dptr;         /* Pointer to navigate badd */
  FILE *ifp=NULL;        /* Input file pointer */

  /*
   * Open input file
   */

  if(!(ifp = open_readfile(inname))) {
    fprintf(stderr,"ERROR: read_bad_ext\n");
    return 1;
  }

  /*
   * Read in data
   */

  if(no_error) {
    aptr = bada;
    bptr = badb;
    cptr = badc;
    dptr = badd;
    while(fgets(line,MAXC,ifp) != NULL && no_error) {
      if(line[0] != comment) {
	if(sscanf(line,"%f %f %f %f %f %f",&junk,&junk,
		  &aptr->flux,&bptr->flux,&cptr->flux,&dptr->flux) != 6) {
	  fprintf(stderr,"ERROR: read_bad_ext. Bad input format.\n");
	  no_error = 0;
	  break;
	}
	else {
	  aptr++;
	  bptr++;
	  cptr++;
	  dptr++;
	}
      }
    }
  }

  /*
   * Clean up and exit
   */

  if(ifp)
    fclose(ifp);

  if(no_error) {
    printf("\nread_bad_ext: Read individual flags from %s.\n",inname);
    return 0;
  }
  else {
    fprintf(stderr,"ERROR: read_bad_ext.\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function flag_bad
 *
 * Flags the bad days in all the Fluxrec arrays by changing their
 *  match variables to -1.  There are two different approaches used in
 *  this function.  If setup->flagbad == 1, then just use the second
 *  column in badfilename to flag days as a whole.  If setup->flagbad == 2, 
 *  then use the third through sixth columns in badfilename to flag
 *  the 1608 curves individually.  In this case, only flag a 1634/1635
 *  point if all four 1608 points for that day are bad.
 *
 * Inputs: Fluxrec *fl34       1634 light curve
 *         Fluxrec *fl35       1635 light curve
 *         Fluxrec *fl08[]     1608 light curves
 *         int nlines          number of points in the light curves
 *         int *flagbad        type of flagging to do
 *         char *badfilename   file containing list of bad days
 *         int *nbad           number of bad days in each curve (set by
 *                              this function)
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 * v21Aug98 CDF, Move function that reads in bad day list from lcurve.c
 *                into this function.
 * v01Sep98 CDF, Modify to allow different flaggings for each 1608 curve.
 */

int flag_bad(Fluxrec *fl34, Fluxrec *fl35, Fluxrec *fl08[], int nlines,
	     int *flagbad, char *badfilename, int *nbad)
{
  int i;                   /* Looping variable */
  int no_error=1;          /* Flag set to 0 on error */
  int nlbad=0;             /* Number of data lines in bad day file */
  int indbad;              /* Number of bad individual point on a given day */
  int na,nb,nc,nd;         /* Temporary containers for nbad */
  Fluxrec *baddays=NULL;   /* List of bad days */
  Fluxrec *bada=NULL;      /* List of bad days for component A */
  Fluxrec *badb=NULL;      /* List of bad days for component B */
  Fluxrec *badc=NULL;      /* List of bad days for component C */
  Fluxrec *badd=NULL;      /* List of bad days for component D */
  Fluxrec *p34,*p35;       /* Pointers for navigating fl34 and fl35 */
  Fluxrec *pa,*pb,*pc,*pd; /* Pointers for navigating 1608 curves */
  Fluxrec *pbad;           /* Pointer for navigating bad day list */
  Fluxrec *ba,*bb,*bc,*bd; /* Pointers for navigating ind. bad day lists */

  /*
   * Read in list of bad days.  First read in the total day flagging
   *  list whether or not we're doing individual curve flagging.  Note
   *  that this also sets nlbad.
   */

  if(!(baddays = read_fluxrec(badfilename,'#',&nlbad))) {
    fprintf(stderr,"ERROR: flag_bad\n");
    return 1;
  }
  
  /*
   * Error check on number of lines in bad day list
   */

  if(nlbad != nlines) {
    fprintf(stderr,"ERROR: flag_bad.  ");
    fprintf(stderr,"Bad day list size does not match data file size.\n");
    no_error = 0;
  }

  /*
   * Now, if flagging individual curves, allocate memory for the
   *  curves and call read_bad_ext to fill the lists.
   */

  if(*flagbad == 2 && no_error) {

    /*
     * Allocate memory
     */

    if(!(bada = new_fluxrec(nlbad)))
      no_error = 0;
    if(!(badb = new_fluxrec(nlbad)))
      no_error = 0;
    if(!(badc = new_fluxrec(nlbad)))
      no_error = 0;
    if(!(badd = new_fluxrec(nlbad)))
      no_error = 0;

    /*
     * Now fill the lists
     */

    if(no_error) {
      if(read_bad_ext(badfilename,'#',bada,badb,badc,badd)) {
	fprintf(stderr,"ERROR: flag_bad. Cannot flag individual curves.\n");
	fprintf(stderr,"  *** Will try total day flagging. ***\n");
	*flagbad = 1;
      }
    }
  }

  /*
   * Flag the bad days in the list
   */

  if(no_error) {
    switch(*flagbad) {
    case 2:

      /*
       * If *flagbad == 2 then do individual curve flagging
       */

      na = nb = nc = nd = 0;
      for(i=0,p34=fl34,p35=fl35,pa=fl08[0],pb=fl08[1],pc=fl08[2],pd=fl08[3],
	    ba=bada,bb=badb,bc=badc,bd=badd; 
	  i<nlines; i++,p34++,p35++,pa++,pb++,pc++,pd++,ba++,bb++,bc++,bd++) {
	indbad = 0;
	if(ba->flux > 0) {
	  pa->match = -1;
	  na++;
	  indbad++;
	}
	if(bb->flux > 0) {
	  pb->match = -1;
	  nb++;
	  indbad++;
	}
	if(bc->flux > 0) {
	  pc->match = -1;
	  nc++;
	  indbad++;
	}
	if(bd->flux > 0) {
	  pd->match = -1;
	  nd++;
	  indbad++;
	}
	if(indbad == 4) {
	  p34->match = -1;
	  p35->match = -1;
	}
      }
      break;

    default:
      
      /*
       * The default case is to flag days as a whole
       */
      
      na = 0;
      for(i=0,p34=fl34,p35=fl35,pa=fl08[0],pb=fl08[1],pc=fl08[2],pd=fl08[3],
	    pbad=baddays; 
	  i<nlines; i++,p34++,p35++,pa++,pb++,pc++,pd++,pbad++) {
	if(pbad->flux > 0) {
	  p34->match = -1;
	  p35->match = -1;
	  pa->match = -1;
	  pb->match = -1;
	  pc->match = -1;
	  pd->match = -1;
	  na++;
	}
      }
      nb = nc = nd = na;
    }
  }

  /*
   * Clean up and return
   */

  baddays = del_fluxrec(baddays);
  bada = del_fluxrec(bada);
  badb = del_fluxrec(badb);
  badc = del_fluxrec(badc);
  badd = del_fluxrec(badd);

  if(no_error) { 
    nbad[0] = na;
    nbad[1] = nb;
    nbad[2] = nc;
    nbad[3] = nd;
    return 0;
  }
  else {
    fprintf(stderr,"ERROR: flag_bad.\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function make_flat
 *
 * Constructs a light curve "flat field" from the two supposedly
 *  non-varying sources, 1634 and 1635.  It does this by normalizing each
 *  of the light curves, then taking the average of the two normalized curves 
 *  for each day.  The normalization is carried out in one of three ways,
 *  which are selected for using the passed meanchoice variable:
 *
 *     meanchoice == 0   --  Use the total flux density from the source
 *                            model.
 *     meanchoice == 1   --  Determine the mean from all points in the
 *                            light curve.
 *     meanchoice == 2   --  Determine separate means for the A, BnA and B 
 *                            configurations, using the norm_config function.
 *
 * Called by lcurve.c, mk_monte_lcurve.c, and mkrand.c
 *
 * Inputs: Fluxrec *fl34       1634 light curve
 *         Fluxrec *fl35       1635 light curve
 *         int nlines          size of array
 *         int meanchoice      normalization method
 *
 * Output: float *flat         flat field -- NULL on error
 *
 * v30Aug98 CDF, Set curve(s) used to make the flat-field curve based on
 *                the FLATCHOICE value.
 *                1 ==> 1634 only, 2 ==> 1635 only, 3 ==> both.
 * v13Jun99 CDF, Added the choice of normalization methods via a new passed
 *                parameter called meanchoice
 */

float *make_flat(Fluxrec *fl34, Fluxrec *fl35, int nlines, int meanchoice)
{
  int i;                 /* Looping variable */
  int no_error=1;        /* Flag set to 0 on error */
  float *flat=NULL;      /* Flat field curve */
  float *flptr;          /* Pointer for navigating flat array */
  Fluxrec *n34=NULL;     /* Normalized 1634 light curve */
  Fluxrec *n35=NULL;     /* Normalized 1635 light curve */
  Fluxrec *nptr1,*nptr2; /* Pointers for navigating float arrays */

  /*
   * Allocate arrays
   */

  if(!(flat = new_array(nlines,1))) {
    fprintf(stderr,"ERROR: make_flat\n");
    return NULL;
  }

  /*
   * Normalize the curves, with the method chosen via the meanchoice
   *  parameter.
   */

  switch(meanchoice) {
  case 0:
    if(!(n34 = norm_constant(fl34,nlines,F34_0)))
      no_error = 0;
    if(!(n35 = norm_constant(fl35,nlines,F35_0)))
      no_error = 0;
    break;
  case 1:
    if(!(n34 = norm_curve(fl34,nlines,"1634",1)))
      no_error = 0;
    if(!(n35 = norm_curve(fl35,nlines,"1635",1)))
      no_error = 0;
    break;
  case 2:
    if(!(n34 = norm_config(fl34,nlines,"1634",1)))
      no_error = 0;
    if(!(n35 = norm_config(fl35,nlines,"1635",1)))
      no_error = 0;
    break;
  default:
    fprintf(stderr,"\n*** WARNING: make_flat - bad value for meanchoice.\n");
    no_error = 0;
  }

  /*
   * Make the flat
   */

  if(no_error)
    for(i=0,nptr1=n34,nptr2=n35,flptr=flat; i<nlines; 
	i++,nptr1++,nptr2++,flptr++) {
      switch(FLATCHOICE) {
      case 1:
	*flptr = nptr1->flux;
	break;
      case 2:
	*flptr = nptr2->flux;
	break;
      case 3:
	*flptr = (nptr1->flux + nptr2->flux)/2.0;
	break;
      default:
	*flptr = (nptr1->flux + nptr2->flux)/2.0;
      }
    }

  /*
   * Clean up
   */

  n34 = del_fluxrec(n34);
  n35 = del_fluxrec(n35);

  if(no_error)
    return flat;
  else {
    fprintf(stderr,"ERROR: make_flat\n");
    return NULL;
  }
}

/*.......................................................................
 *
 * Function flat_field
 *
 * Flat-fields the light curve, using the flat created by calling make_flat.
 *
 * Inputs: Fluxrec *raw        data to be flat-fielded
 *         float *flat         flat field
 *         int nlines          number of data points
 *         float fracrms       fractional rms scatter in 1634/1635 ratio
 *                             -- contributes to error bars.
 *         char *outname       name of output file
 *         int doprint         flag set to 1 for printing to output file
 *
 * Output   Fluxrec *ff        flat-fielded light curve -- NULL on error
 *
 * v02Oct98 CDF, Added output option.
 */

Fluxrec *flat_field(Fluxrec *raw, float *flat, int nlines, float fracrms,
		    char *outname, int doprint)
{
  int i;             /* Looping variable */
  float *fptr;       /* Pointer to navigate flat */
  Fluxrec *ff=NULL;  /* Flat-fielded data */
  Fluxrec *rptr;     /* Pointer to navigate raw */
  Fluxrec *ffptr;    /* Pointer to navigate ff */
  FILE *ofp=NULL;    /* Output file pointer */

  /*
   * Open the output file
   */
  if(doprint) {
    if(!(ofp = fopen(outname,"w"))) {
      fprintf(stderr,"ERROR: flat_field\n");
      return NULL;
    }
  }

  /*
   * Allocate memory for the flat-fielded curves
   */

  if(!(ff = new_fluxrec(nlines))) {
    fprintf(stderr,"ERROR: flat_field\n");
    if(ofp)
      fclose(ofp);
    return NULL;
  }

  /*
   * Do the flat-fielding for the light curves.
   */

  for(i=0,rptr=raw,fptr=flat,ffptr=ff; i<nlines; 
      i++,rptr++,fptr++,ffptr++) {
    *ffptr = *rptr;
    ffptr->flux = (rptr->flux / *fptr);
#if 1
    ffptr->err = sqrt(rptr->err * rptr->err + 
		      fracrms * fracrms * rptr->flux * rptr->flux);
#else
    ffptr->err = rptr->err;
#endif
    if(doprint)
      fprintf(ofp,"%7.2f %6.4f %7.5f %d\n",ffptr->day,ffptr->flux,ffptr->err,
	      ffptr->match);
  }

  /*
   * Clean up and exit
   */

  if(ofp)
    fclose(ofp);
  return ff;
}



/*.......................................................................
 *
 * Function norm_constant
 *
 * Normalizes a curve by dividing each point by a constant value
 *
 * Inputs: Fluxrec *raw        raw input curve
 *         int nlines          number of points in curve
 *         float constant      constant normalization
 *
 * Output: Fluxrec *norm       normalized curve -- NULL on error
 *
 */

Fluxrec *norm_constant(Fluxrec *raw, int nlines, float constant)
{
  int i;                /* Looping variable */
  int n=0;              /* Number of valid points in curve */
  Fluxrec *norm=NULL;   /* Normalized curve */
  Fluxrec *nptr;        /* Pointer for navigating norm */
  Fluxrec *flptr;       /* Pointer for navigating raw light curve */

  if(!(norm = new_fluxrec(nlines))) {
    fprintf(stderr,"ERROR: norm_constant\n");
    return NULL;
  }

  /*
   * Normalize
   */

  for(i=0,flptr=raw,nptr=norm; i<nlines; i++,flptr++,nptr++) {
    *nptr = *flptr;
    nptr->flux /= constant;
    nptr->err /= constant;
  }

  return norm;
}

/*.......................................................................
 *
 * Function norm_curve
 *
 * Normalizes a curve by dividing each point by the mean value of the curve
 *
 * Inputs: Fluxrec *raw        raw input curve
 *         int nlines          number of points in curve
 *         char *source        name of source
 *         int doprint         if == 1 then print out mean value
 *
 * Output: Fluxrec *norm       normalized curve -- NULL on error
 *
 * v12Feb97 CDF, Added print control parameter "doprint"
 * v18Jun98 CDF, Don't include bad days (match == -1) in calculation of mean
 */

Fluxrec *norm_curve(Fluxrec *raw, int nlines, char *source, int doprint)
{
  int i;                /* Looping variable */
  int n=0;              /* Number of valid points in curve */
  float mean=0.0;       /* Mean value */
  Fluxrec *norm=NULL;   /* Normalized curve */
  Fluxrec *nptr;        /* Pointer for navigating norm */
  Fluxrec *flptr;       /* Pointer for navigating raw light curve */

  if(!(norm = new_fluxrec(nlines))) {
    fprintf(stderr,"ERROR: norm_curve\n");
    return NULL;
  }

  /*
   * Calculate mean
   */

  for(i=0,flptr=raw; i<nlines; i++,flptr++) {
    if(flptr->match != -1) {
      mean += flptr->flux;
      n++;
    }
  }

  mean /= n;

  if(doprint)
    printf("%s has mean flux of %7.2f mJy\n",source,mean);

  /*
   * Normalize
   */

  for(i=0,flptr=raw,nptr=norm; i<nlines; i++,flptr++,nptr++) {
    *nptr = *flptr;
    nptr->flux /= mean;
    nptr->err /= mean;
  }

  return norm;
}

/*.......................................................................
 *
 * Function norm_zero_mean
 *
 * Creates a zero-mean curve by first normalizing the curve by dividing
 *  by its mean value and then subtracting 1.  The reason for doing this
 *  slightly roundabout procedure is to bring the fractional variations
 *  in the light curves to the same level.
 *
 * Inputs: Fluxrec *raw        raw input curve
 *         int nlines          number of points in curve
 *
 * Output: Fluxrec *zmean      zero-mean curve -- NULL on error
 *
 */

Fluxrec *norm_zero_mean(Fluxrec *raw, int nlines)
{
  int i;                /* Looping variable */
  int n=0;              /* Number of valid points in curve */
  float mean=0.0;       /* Mean value */
  Fluxrec *zmean=NULL;  /* Zero-mean curve */
  Fluxrec *zptr;        /* Pointer for navigating zmean */
  Fluxrec *flptr;       /* Pointer for navigating raw light curve */

  /*
   * Allocate memory for zero-mean array
   */

  if(!(zmean = new_fluxrec(nlines))) {
    fprintf(stderr,"ERROR: norm_zero_mean\n");
    return NULL;
  }

  /*
   * Calculate mean
   */

  for(i=0,flptr=raw; i<nlines; i++,flptr++) {
    if(flptr->match != -1) {
      mean += flptr->flux;
      n++;
    }
  }

  mean /= n;

  /*
   * Normalize and then subtract 1 to get a zero-mean curve
   */

  for(i=0,flptr=raw,zptr=zmean; i<nlines; i++,flptr++,zptr++) {
    *zptr = *flptr;
    zptr->flux = (flptr->flux / mean) - 1.0;
    zptr->err /= mean;
  }

  return zmean;
}

/*.......................................................................
 *
 * Function norm_config
 *
 * Normalizes a curve by dividing each point by the mean value of the curve
 *  in the appropriate configuration, i.e. the curve is normalized
 *  with separate values for the A, BnA and B configuration sections.
 *
 * Inputs: Fluxrec *raw        raw input curve
 *         int nlines          number of points in curve
 *         char *source        name of source
 *         int doprint         if == 1 then print out mean value
 *
 * Output: Fluxrec *norm       normalized curve -- NULL on error
 *
 * v12Feb97 CDF, added print control parameter "doprint"
 */

Fluxrec *norm_config(Fluxrec *raw, int nlines, char *source, int doprint)
{
  int i;             /* Looping variable */
  int acount=0;      /* Number of points in A config */
  int bcount=0;      /* Number of points in B config */
  int bnacount=0;    /* Number of points in BnA config */
  float amean=0.0;   /* Mean value for A config */
  float bmean=0.0;   /* Mean value for B config */
  float bnamean=0.0; /* Mean value for BnA config */
  Fluxrec *norm;     /* Normalized curve */
  Fluxrec *nptr;     /* Pointer for navigating norm */
  Fluxrec *flptr;    /* Pointer for navigating raw light curve */

  if(!(norm = new_fluxrec(nlines))) {
    fprintf(stderr,"ERROR: norm_config\n");
    return NULL;
  }

  /*
   * Calculate mean values for A, BnA and B configurations
   */

  for(i=0,flptr=raw; i<nlines; i++,flptr++) {
    if(flptr->day < AEND && flptr->day > ASTART && flptr->match != -1) {
      acount++;
      amean += flptr->flux;
    }
    else if(flptr->day > AEND && flptr->day < BSTART && flptr->match != -1) {
      bnacount++;
      bnamean += flptr->flux;
    }
    else if(flptr->day > BSTART && flptr->match != -1) {
      bcount++;
      bmean += flptr->flux;
    }
  }

  amean /= acount;
  bmean /= bcount;
  bnamean /= bnacount;

  if(doprint) {
    printf("\n%s has mean flux density of %7.2f mJy in A config (%d points)\n",
	   source,amean,acount);
    printf("%s has mean flux density of %7.2f mJy in BnA config (%d points)\n",
	   source,bnamean,bnacount);
    printf("\%s has mean flux density of %7.2f mJy in B config (%d points)\n",
	   source,bmean,bcount);
  }

  /*
   * Normalize
   */

  for(i=0,flptr=raw,nptr=norm; i<nlines; i++,flptr++,nptr++) {
    *nptr = *flptr;
    if(flptr->day < AEND) {
      nptr->flux = flptr->flux / amean;
      nptr->err = flptr->err / amean;
    }
    else if(flptr->day > AEND && flptr->day < BSTART) {
      nptr->flux = flptr->flux / bnamean;
      nptr->err = flptr->err / bnamean;
    }
    else {
      nptr->flux = flptr->flux / bmean;
      nptr->err = flptr->err / bmean;
    }
  }

  return norm;
}

/*.......................................................................
 *
 * Function calc_flrat
 *
 * Calculates flux ratios between two flux curves
 *
 * Inputs: Fluxrec *flux1      curve 1 (numerator)
 *         Fluxrec *flux2      curve 2 (denominator)
 *         int size            size of arrays
 *
 * Output: Fluxrec *flrat      flux ratio curve
 *
 */

Fluxrec *calc_flrat(Fluxrec *flux1, Fluxrec *flux2, int size)
{
  int i;                 /* Looping variable */
  int divzero=0;         /* Flag set to 1 if attempted div. by 0 */
  Fluxrec *flrat=NULL;   /* Flux ratio curve to be returned */
  Fluxrec *fl1,*fl2;     /* Pointers to navigate flux1 and flux2 */
  Fluxrec *flr;          /* Pointer to navigate flrat */

  if(!(flrat = new_fluxrec(size))) {
    fprintf(stderr,"ERROR: calc_flrat\n");
    return NULL;
  }

  for(i=0,fl1=flux1,fl2=flux2,flr=flrat; i<size; i++,fl1++,fl2++,flr++) {
    if(fl2->flux == 0.0) {
      divzero = 1;
      flr->flux = 0.0;
      flr->err = 0.0;
      flr->day = fl1->day;
    }
    else {
      flr->flux = log10(fl1->flux / fl2->flux);
      flr->err = (fl1->flux / fl2->flux) * 
	sqrt((fl1->err*fl1->err/(fl1->flux*fl1->flux)) +
	     (fl2->err*fl2->err/(fl2->flux*fl2->flux)));
      flr->day = fl1->day;
    }
  }

  if(divzero)
    fprintf(stderr,"WARNING: Division by zero attempted in calc_flrat\n");

  return flrat;
}

/*.......................................................................
 *
 * Function set_mu0
 *
 * Calculates the mean values in the overlap regions between the 1608
 *  light curves after they have been shifted by a initial guess for the
 *  time delays.  If no initial guess is indicated, then the program uses
 *  the innermost 50% of the points in each light curve to calculate the
 *  mean values.  After the mean values are calculated, the ratios of
 *  curves A, C, and D to curve B are placed in the setup container as
 *  initial guesses for the curve-fitting routines.  These ratios are
 *  stored in setup->mu0.
 *
 * Inputs: Fluxrec *fl08[]     1608 light curves
 *         int nlines          number of points in light curve
 *         Setup *setup        setup information.  Note that 
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int set_mu0(Fluxrec *fl08[], int nlines, Setup *setup)
{
  int i;                /* Looping variable */
  int no_error=1;       /* Flag set to 0 on error */
  int startindex;       /* Array index corresponding to startday */
  int endindex;         /* Array index corresponding to startday */
  float startday;       /* Starting day of overlap region */
  float endday;         /* Starting day of overlap region */
  float amean;          /* Mean of curve A in the overlap region */
  float bmean;          /* Mean of curve B in the overlap region */
  float cmean;          /* Mean of curve C in the overlap region */
  float dmean;          /* Mean of curve D in the overlap region */
  float rms;            /* RMS in light-curve section */
  Fluxrec *fptr;        /* Pointer to navigate the fl08 arrays */

  /*
   * First check if the setup->nmu value has been set by the user.
   *  If it hasn't, set it to the default NFLUXSTEP.
   */

  if(setup->nmu == 0)
    setup->nmu = NFLUXSTEP;
  printf("\nset_mu0: Number of magnification steps set to %d.\n",
	 setup->nmu);

  /*
   * Now check to see if the setup->mu0 values have already been
   *  set by the user.  If they have, don't bother continuing.
   */

  if(setup->mu0[0] > 0.0) {
    printf("set_mu0:  Initial guesses for magnifications have ");
    printf("already been set:\n");
    printf("   mu_AB = %7.4f\n",setup->mu0[0]);
    printf("   mu_CB = %7.4f\n",setup->mu0[2]);
    printf("   mu_DB = %7.4f\n",setup->mu0[3]);
    return 0;
  }

  /*
   * Calculate the means using the inner 50% of points from each light 
   *  curve.
   */

  printf("set_mu0: ------------------------------------------------------\n");
  printf("set_mu0: Calculating flux ratios for inner 50%% of light curves.\n");
  fptr = fl08[0];
  startindex = (int) floor(nlines/4);
  endindex = (int) floor(3*nlines/4);
  startday = (fptr+startindex-1)->day;
  endday = (fptr+endindex-1)->day;
  if(calc_mean_dt(fl08[0],nlines,&amean,&rms,startday,endday,0.0))
    no_error = 0;
  if(calc_mean_dt(fl08[1],nlines,&bmean,&rms,startday,endday,0.0))
    no_error = 0;
  if(calc_mean_dt(fl08[2],nlines,&cmean,&rms,startday,endday,0.0))
    no_error = 0;
  if(calc_mean_dt(fl08[3],nlines,&dmean,&rms,startday,endday,0.0))
    no_error = 0;
  printf("set_mu0: Means are: %6.3f %6.3f %6.3f %6.3f.\n",amean,bmean,cmean,
	 dmean);
  printf("set_mu0: Flux ratios are: %6.4f %6.4f %6.4f.\n",amean/bmean,
	 cmean/bmean,dmean/bmean);
  printf("set_mu0: ------------------------------------------------------\n");

  /*
   * If setup->dooverlap is set to NO (the default), stop here.
   */
  

  /*
   * If, however, setup->dooverlap is set to YES, then shift the day
   *  curves and clip all points that fall outside the overlap regions.
   */

  fptr = fl08[1];
  startday = fptr->day;
  endday = (fptr+nlines-1)->day - DLAG;
  printf("set_mu0: Calculating flux ratios in overlap region.\n");
  printf("set_mu0: Start and end of overlap region are %6.1f %6.1f.\n",
	 startday,endday);
  if(calc_mean_dt(fl08[0],nlines,&amean,&rms,startday,endday,ALAG))
    no_error = 0;
  if(calc_mean_dt(fl08[1],nlines,&bmean,&rms,startday,endday,0.0))
    no_error = 0;
  if(calc_mean_dt(fl08[2],nlines,&cmean,&rms,startday,endday,CLAG))
    no_error = 0;
  if(calc_mean_dt(fl08[3],nlines,&dmean,&rms,startday,endday,DLAG))
    no_error = 0;
  printf("set_mu0: Means are: %6.3f %6.3f %6.3f %6.3f.\n",amean,bmean,cmean,
	 dmean);
  printf("set_mu0: Flux ratios are: %6.4f %6.4f %6.4f.\n",amean/bmean,
	 cmean/bmean,dmean/bmean);
  printf("set_mu0: ------------------------------------------------------\n");

  /*
   * Fill the setup->mu0 containers
   */

  setup->mu0[0] = amean/bmean;
  setup->mu0[1] = bmean/bmean;
  setup->mu0[2] = cmean/bmean;
  setup->mu0[3] = dmean/bmean;

  return 0;
}

/*.......................................................................
 *
 * Function set_mu_grid
 *
 * Calculates the default values of the flux ratios between the input
 *  curves.  The function uses the innermost 50% of the points in each 
 *  light curve to calculate the mean values.  After the mean values are 
 *  calculated, the ratios of the mean values are placed in the setup 
 *  container as initial guesses for the curve-fitting routines.  
 *  These ratios are stored in setup->mu0.
 *
 * Inputs: Fluxrec *lc[]       input light curves
 *         int *npoints        number of points in each light curves
 *         Setup *setup        setup information.  Note that 
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int set_mu_grid(Fluxrec *lc[], int *npoints, Setup *setup)
{
  int i;                /* Looping variable */
  int no_error=1;       /* Flag set to 0 on error */
  int startindex;       /* Array index corresponding to startday */
  int endindex;         /* Array index corresponding to startday */
  float startday;       /* Starting day of overlap region */
  float endday;         /* Starting day of overlap region */
  float *mean;          /* Array to hold means of curves */
  float rms;            /* RMS in light-curve section */
  Fluxrec *fptr;        /* Pointer to navigate the lc arrays */

  /*
   * First check if the setup->nmu value has been set by the user.
   *  If it hasn't, set it to the default NFLUXSTEP.
   */

  if(setup->nmu == 0)
    setup->nmu = NFLUXSTEP;
  printf("\nset_mu_grid: Number of magnification steps set to %d.\n",
	 setup->nmu);

  /*
   * Now check to see if the setup->mu0 values have already been
   *  set by the user.  If they have, don't bother continuing.
   */

  if(setup->mu0[0] > 0.0) {
    printf("set_mu_grid:  Initial guesses for magnifications have ");
    printf("already been set:\n");
    for(i=0; i<setup->ncurves; i++) {
      printf("   mu0_[%d] = %7.4f\n",i,setup->mu0[0]);
    }
    return 0;
  }

  /*
   * Allocate memory for mean array
   */

  if(!(mean = new_array(setup->ncurves,1))) {
    fprintf(stderr,"ERROR: set_mu_grid.\n");
    return 1;
  }

  /*
   * Calculate the means using the inner 50% of points from each light 
   *  curve.
   */

  printf("set_mu_grid:  --------------------------------------------------\n");
  printf("set_mu_grid: Calculating means for inner 50%% of light curves.\n");
  for(i=0; i<setup->ncurves; i++) {
    fptr = lc[i];
    startindex = (int) floor(npoints[i] / 4);
    endindex = (int) floor(3 * npoints[i] / 4);
    startday = (fptr+startindex-1)->day;
    endday = (fptr+endindex-1)->day;
    if(calc_mean_dt(lc[i],npoints[i],mean+i,&rms,startday,endday,0.0))
      no_error = 0;
  }

  /*
   * Print out results and fill the setup->mu0 containers
   */

  printf("set_mu_grid: Means are: ");
  for(i=0; i<setup->ncurves; i++)
    printf("%6.3f ",mean[i]);
  printf("\nset_mu_grid: Flux ratios are: ");
  for(i=0; i<setup->ncurves; i++) {
    setup->mu0[i] = mean[i]/mean[0];
    printf("%6.4f ",setup->mu0[i]);
  }
  printf("\n\n");

  return 0;
}

/*.......................................................................
 *
 * Function set_tau_grid
 *
 * Sets the default values for the delays to be used in the grid search
 *  for the best-fit delays.  The values are only set if not already set
 *  by the input setup file.  The values are stored in setup->tu0.
 *
 * Inputs: Fluxrec *lc[]       input light curves
 *         int *npoints        number of points in each light curves
 *         int *index          array showing which curves are being compared
 *         Setup *setup        setup information.  Note that 
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int set_tau_grid(Fluxrec *lc[], int *npoints, int *index, Setup *setup)
{
  int i;                /* Looping variable */
  int no_error=1;       /* Flag set to 0 on error */
  float startday;       /* Starting day of light curve */
  float endday;         /* Starting day of light curve */
  float ttotal;         /* Total length of run */
  float *mean;          /* Array to hold means of curves */
  float rms;            /* RMS in light-curve section */
  char line[MAXC];      /* General string for reading input */
  Fluxrec *fptr;        /* Pointer to navigate the lc arrays */

  /*
   * Print out information about the input curves
   */

  printf("\nset_tau_grid: Curve  Start    End    Midpt   Length  <dt> \n");
  printf("set_tau_grid: -----  ------  ------  ------  ------  -----\n");
  for(i=0; i<setup->ncurves; i++) {
    startday = lc[index[i]]->day;
    endday = (lc[index[i]]+npoints[index[i]]-1)->day;
    ttotal = endday - startday;
    printf("set_tau_grid: %5d  %6.1f  %6.1f  %6.1f  %5.1f   %4.1f\n",
	   i+1,startday,endday,(startday+endday)/2.0,ttotal,
	   ttotal/npoints[index[i]]);
  }

  /*
   * See if the values of tau0 have been set in the input file.  If they
   *  have been, just echo the values.
   */

  if(setup->tauset == YES) {
    printf("set_tau_grid: Using values for tau0 set in input file.\n");
    printf("set_tau_grid: ");
    for(i=0; i<setup->ncurves; i++)
      printf("%6.1f  ",setup->tau0[i]);
  }

  /*
   * If tau0 values have not been set, query the user for them
   */

  else {
    printf("\nset_tau_grid: Values of tau0 (delays) have not been set.\n");
    for(i=0; i<setup->ncurves; i++) {
      printf("set_tau_grid: Enter value of tau0 to use for curve %d: [%6.1f] ",
	     i+1,setup->tau0[i]);
      fgets(line,MAXC,stdin);
      if(line[0] != '\n') {
	while(sscanf(line,"%f",&setup->tau0[i]) != 1) {
	  fprintf(stderr," ERROR. Bad input.  Enter value again:  ");
	  fgets(line,MAXC,stdin);
	}
      }
    }
  }

  /*
   * Get the stepsize for the tau grid if it has not been set already
   */

  if(setup->dtau == 0) {
   printf("\nset_tau_grid: Stepsize to be used in delay grid (dtau) has ");
   printf("not been set.\n");
   setup->dtau = 1.0;
   printf("set_tau_grid: Enter value of dtau: [%5.2f] ",setup->dtau);
   fgets(line,MAXC,stdin);
   if(line[0] != '\n') {
     while(sscanf(line,"%lf",&setup->dtau) != 1 || setup->dtau <= 0.0) {
       fprintf(stderr," ERROR. Bad input.  Enter value again:  ");
       fgets(line,MAXC,stdin);
     }
   }
  }

  /*
   * Calculate the number of time delays to consider if this value has
   *  not already been set in the Setup function.  The default value is
   *  that which will include 1/4 of the total length of the observations
   *  on each side of tau_0, i.e., from 1/4 to 3/4 of the total length.
   *  Use the length of the unshifted curve, if the curves are of different
   *  lengths.
   */

  if(setup->ntau == 0 && no_error) {
    setup->ntau = floor(ttotal / (4.0 * setup->dtau));
    printf("set_tau_grid: Enter number of time delay steps to take on\n");
    printf("set_tau_grid:   either side of tau0: [%d] ",setup->ntau);
    fgets(line,MAXC,stdin);
    if(line[0] != '\n') {
      while(sscanf(line,"%d",&setup->ntau) != 1 || setup->ntau <= 0) {
	fprintf(stderr,"ERROR: Invalid input for number of steps.  ");
	fprintf(stderr,"Enter value again.  ");
	fgets(line,MAXC,stdin);
      }
    }
  }

  printf("\n");
#if 0
  /*
   * Put the tau information into the Prange container.  This
   *  information is stored in the Setup container, and the assignment
   *  is controlled by the index array.  Note that index[0] will always
   *  contain the index of the UNSHIFTED/UNMAGNIFIED curve.  This means
   *  that we only load the values for index[1] -> index[ncurves-1].
   */

  if(no_error) {
    for(i=0,pptr=tau0; i<ncurves-1; i++,pptr++) {
      pptr->val0 = setup->tau0[index[i+1]];
      pptr->nval = setup->ntau;
      pptr->dval = DAYSTEP;
      pptr->minstep = ((int) (pptr->val0/pptr->dval)) - pptr->nval;
      pptr->maxstep = ((int) (pptr->val0/pptr->dval)) + pptr->nval;
      printf(" disp_setup: tau0=%6.1f, tau_min=%6.1f, tau_max=%6.1f, ",
	     pptr->val0,pptr->minstep * pptr->dval,
	     pptr->maxstep * pptr->dval);
      printf("dtau=%5.2f, ntau=%d\n",pptr->dval,2 * pptr->nval + 1);
    }
  }

  /*
   * First check if the setup->nmu value has been set by the user.
   *  If it hasn't, set it to the default NFLUXSTEP.
   */

  if(setup->nmu == 0)
    setup->nmu = NFLUXSTEP;
  printf("\nset_mu_grid: Number of magnification steps set to %d.\n",
	 setup->nmu);

  /*
   * Now check to see if the setup->mu0 values have already been
   *  set by the user.  If they have, don't bother continuing.
   */

  if(setup->mu0[0] > 0.0) {
    printf("set_mu_grid:  Initial guesses for magnifications have ");
    printf("already been set:\n");
    for(i=0; i<setup->ncurves; i++) {
      printf("   mu0_[%d] = %7.4f\n",i,setup->mu0[0]);
    }
    return 0;
  }

  /*
   * Allocate memory for mean array
   */

  if(!(mean = new_array(setup->ncurves,1))) {
    fprintf(stderr,"ERROR: set_mu_grid.\n");
    return 1;
  }

  /*
   * Calculate the means using the inner 50% of points from each light 
   *  curve.
   */

  printf("set_mu_grid:  --------------------------------------------------\n");
  printf("set_mu_grid: Calculating means for inner 50%% of light curves.\n");
  for(i=0; i<setup->ncurves; i++) {
    fptr = lc[i];
    startindex = (int) floor(npoints[i] / 4);
    endindex = (int) floor(3 * npoints[i] / 4);
    startday = (fptr+startindex-1)->day;
    endday = (fptr+endindex-1)->day;
    if(calc_mean_dt(lc[i],npoints[i],mean+i,&rms,startday,endday,0.0))
      no_error = 0;
  }

  /*
   * Print out results and fill the setup->mu0 containers
   */

  printf("set_mu_grid: Means are: ");
  for(i=0; i<setup->ncurves; i++)
    printf("%6.3f ",mean[i]);
  printf("\nset_mu_grid: Flux ratios are: ");
  for(i=0; i<setup->ncurves; i++) {
    setup->mu0[i] = mean[i]/mean[0];
    printf("%6.4f ",setup->mu0[i]);
  }
  printf("\n");
#endif
  return 0;
}


/*.......................................................................
 *
 * Function ratio_err
 *
 * Calculates the rms scatter in the ratio of two fluxes.
 *
 * Inputs: Fluxrec *flux1      first flux
 *         Fluxrec *flux2      second flux
 *         int nlines          number of lines in arrays
 *         float *fracrms      fractional rms in ratio (set by this function)
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> failure
 *
 * v18Jun99 CDF, Correct the input lightcurves for slightly different
 *                mean flux densities in each array configuration BEFORE
 *                taking ratio.  This is necessary since the 1635 mean
 *                flux density changes by several percent between A and
 *                B configurations due to extended structure.  This
 *                correction is accomplished by a call to corr_config.
 */

int ratio_err(Fluxrec *flux1, Fluxrec *flux2, int nlines, float *fracrms)
{
  int i;                /* Looping variable */
  int acount=0;         /* Number of days in A configuration */
  int bcount=0;         /* Number of days in B configuration */
  int bnacount=0;       /* Number of days between A and B configurations */
  int no_error=1;       /* Flag set to 0 on error */
  float amean,bmean;    /* Mean of ratios */
  float arms,brms;      /* RMS scatter in ratios */
  float rmean;          /* Mean of ratio light curve */
  Fluxrec *fptr1;       /* Pointer to navigate first source arrays */
  Fluxrec *fptr2;       /* Pointer to navigate second source arrays */
  Fluxrec *rptr;        /* Pointer to navigate ratio array */
  Fluxrec *corr1=NULL;  /* Corrected light curve for first source */
  Fluxrec *corr2=NULL;  /* Corrected light curve for second source */
  Fluxrec *ratio=NULL;  /* Ratio between the two corrected light curves */
  Fluxrec *aratio=NULL; /* Ratio of flux1 to flux2 for A configuration data */
  Fluxrec *bratio=NULL; /* Ratio of flux1 to flux2 for B configuration data */
  Fluxrec *aptr,*bptr;  /* Pointers to navigate ratio arrays */
#if 0

  /*
   * Correct the input light curves
   */

  if(!(corr1 = corr_config(flux1,nlines)))
    no_error = 0;

  if(!(corr2 = corr_config(flux2,nlines)))
    no_error = 0;

  /*
   * Allocate memory for the ratio array
   */

  if(!(ratio = new_fluxrec(nlines)))
    no_error = 0;

  /*
   * Fill the ratio array
   */

  if(no_error) {
    for(i=0,fptr1=corr1,fptr2=corr2,rptr=ratio; i<nlines; 
	i++,fptr1++,fptr2++,rptr++) {
      *rptr = *fptr1;
      if(fptr2->flux == 0.0)
	rptr->flux = 0.0;
      else
	rptr->flux = fptr1->flux / fptr2->flux;
    }
  }

  /*
   * Calculate the mean and RMS scatter in the ratio array
   */

  if(no_error) {
    if(calc_mean(ratio,nlines,&rmean,fracrms) == 1)  {
      fprintf(stderr,"ERROR: ratio_err.  Setting fracrms = 0.0\n");
      *fracrms = 0.0;
    }
  }

  /*
   * Print out results, after reducing fracrms by a futher factor of
   *  sqrt(2) since the fractional RMS comes from a ratio of two curves and
   *  thus the individual curves should have RMS errors approximately a factor
   *  of sqrt(2) less than the RMS errors on the ratio.
   */

  *fracrms /= sqrt(2.0);

  if(no_error) {
    printf("\nratio_err: RMS scatter in 1634/1635 ratio is %-7.3f%%\n",
	   *fracrms * 100.0 * sqrt(2.0));
    printf("ratio_err: ==> Fractional error in flux densities is %-7.3f%%\n\n",
	   *fracrms * 100.0);
  }
#endif
  /*
   * Calculate number of days in A and B configurations
   */

  for(i=0,fptr1=flux1; i< nlines; i++,fptr1++) {
    if(fptr1->day < AEND && fptr1->day > ASTART && fptr1->match != -1)
      acount++;
    else if(fptr1->day > BSTART && fptr1->match != -1)
      bcount++;
    else
      bnacount++;
  }
  
  /*
   * Allocate arrays
   */

  if(!(aratio = new_fluxrec(acount))) {
    fprintf(stderr,"ERROR:  ratio_err\n");
    return 1;
  }

  if(!(bratio = new_fluxrec(bcount))) {
    fprintf(stderr,"ERROR:  ratio_err\n");
    return 1;
  }

  /*
   * Fill ratio arrays
   */

  aptr = aratio;
  bptr = bratio;
  for(i=0,fptr1=flux1,fptr2=flux2; i<nlines; 
      i++,fptr1++,fptr2++) {
    if(fptr1->day < AEND && fptr1->day > ASTART && fptr1->match != -1) {
      if(fptr2->flux == 0.0)
	aptr->flux = 0.0;
      else
	aptr->flux = fptr1->flux / fptr2->flux;
      aptr++;
    }
    else if(fptr1->day > BSTART && fptr1->match != -1) {
      if(fptr2->flux == 0.0)
	bptr->flux = 0.0;
      else
	bptr->flux = fptr1->flux / fptr2->flux;
      bptr++;
    }
  }

  /*
   * Get mean and rms scatter in ratio arrays
   */

  if(calc_mean(aratio,acount,&amean,&arms) == 1)  {
    fprintf(stderr,"ERROR: ratio_err\n");
    fprintf(stderr,"  Setting arms = 0.0\n");
    arms = 0.0;
    amean = 1.0;
  }

  if(calc_mean(bratio,bcount,&bmean,&brms) == 1)  {
    fprintf(stderr,"ERROR: ratio_err\n");
    fprintf(stderr,"  Setting brms = 0.0\n");
    brms = 0.0;
    bmean = 1.0;
  }

  /*
   * Combine A and B configuration fractional means into total fractional
   *  mean by adding them in quadrature.  Reduce this by a further factor of 
   *  sqrt(2) since the fractional RMS comes from a ratio of two curves and
   *  thus the individual curves should have RMS errors approximately a factor
   *  of sqrt(2) less than the RMS errors on the ratio.
   */

  *fracrms = sqrt((arms*arms/(amean*amean) + brms*brms/(bmean*bmean))/2.0);

  if(no_error) {
    printf("\nratio_err: RMS scatter in 1634/1635 ratio is %-7.3f%%\n",
	   *fracrms * 100.0 * sqrt(2.0));
    printf("ratio_err:    (RMS scatter for A config is %-7.3f%%)\n",
	   100.0*arms/amean);
    printf("ratio_err:    (RMS scatter for B config is %-7.3f%%)\n",
	   100.0*brms/bmean);
    printf("ratio_err: ==> Fractional error in flux densities is %-7.3f%%\n\n",
	   *fracrms * 100.0);
  }

  /*
   * KLUDGE
   */

  *fracrms = 0.008925;
  printf("**** Using kludge value of 0.008925\n");

  /*
   * Clean up
   */

  corr1 = del_fluxrec(corr1);
  corr2 = del_fluxrec(corr2);
  ratio = del_fluxrec(ratio);
  aratio = del_fluxrec(aratio);
  bratio = del_fluxrec(bratio);

  if(no_error)
    return 0;
  else
    return 1;
}

/*.......................................................................
 *
 * Function corr_config
 *
 * Takes an input lightcurve and corrects the flux densities for any change
 *  in median flux density due to array configuration changes.
 *
 * Inputs: Fluxrec *flux       input light curve
 *         int nlines          number of lines in array
 *
 * Output: Fluxrec *corr       corrected light curve
 *
 */

Fluxrec *corr_config(Fluxrec *flux, int nlines)
{
  int i;                   /* Looping variable */
  int acount;              /* Number of valid points in A configuration */
  int bnacount;            /* Number of valid points in BnA configuration */
  int bcount;              /* Number of valid points in B configuration */
  float amed,bmed,bnamed;  /* Medians for the various array configurations */
  Fluxrec *fptr;           /* Pointer to navigate flux */
  Fluxrec *corr=NULL;      /* Corrected light curve */

  /*
   * Calculate number of days in A, BnA, and B configurations
   */

  for(i=0,fptr=flux; i< nlines; i++,fptr++) {
    if(fptr->day < AEND && fptr->day > ASTART && fptr->match != -1)
      acount++;
    else if(fptr->day > BSTART && fptr->match != -1)
      bcount++;
    else if(fptr->day > AEND && fptr->day < BSTART && fptr->match != -1)
      bnacount++;
  }
  
  return corr;
}

/*.......................................................................
 *
 * Function calc_mean
 *
 * Calculates the sample mean and the rms scatter of the sample.
 * Remember that an unbiassed estimate of the variance is given by:
 *
 *               1
 *      s^2 = -------  sum((x_i - mean(x_i))^2)
 *             N - 1
 *
 * Inputs: Fluxrec *flux       array of fluxes
 *         int nlines          number of lines in array (may not equal number
 *                               of measured fluxes)
 *         float *mean         mean of fluxes (set by function)
 *         float *rms          rms of sample (set by function)
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> failure
 *
 * v25Jun98 CDF, Changed criterion for including a point in the calculation
 *                of the mean from flux->flux > 0 to flux->match != -1 
 *                because the cross-correlation curves, which will now be
 *                passed to this function, have both positive and negative
 *                values.
 */

int calc_mean(Fluxrec *flux, int nlines, float *mean, float *rms)
{
  int i;           /* Looping variable */
  int nfluxes=0;   /* Number of measured fluxes */
  float sum;       /* Running sum */
  Fluxrec *dfptr;  /* Pointer to navigate array of fluxes */

  /*
   * First calculate the mean
   */

  sum = 0.0;

  for(i=0,dfptr=flux; i<nlines; i++,dfptr++) {
    if(dfptr->match != -1) {
      sum += dfptr->flux;
      nfluxes++;
    }
  }

  /*
   * Sanity check on number of fluxes
   */

  if(nfluxes == 0) {
    *mean = 0.0;
    *rms = 0.0;
    fprintf(stderr,"ERROR: calc_mean. No input values.\n");
    return 1;
  }
  else if(nfluxes > nlines) {
    fprintf(stderr,"ERROR: calc_mean. More fluxes than expected.\n");
    return 1;
  }

  *mean = sum/nfluxes;

  /*
   * Now find the sample variance
   */

  sum = 0.0;

  if(nfluxes == 1)
    *rms = 0.0;
  else {
    for(i=0,dfptr=flux; i<nlines; i++,dfptr++)
      if(dfptr->match != -1)
	sum += (dfptr->flux - *mean) * (dfptr->flux - *mean);
    *rms = sqrt(sum/(nfluxes-1));
  }

  return 0;
}

/*.......................................................................
 *
 * Function calc_mean_dt
 *
 * Calculates the sample mean and the rms scatter of the sample, in a
 *  time range given by the passed startday and endday parameters.
 * Remember that an unbiassed estimate of the variance is given by:
 *
 *               1
 *      s^2 = -------  sum((x_i - mean(x_i))^2)
 *             N - 1
 *
 * Inputs: Fluxrec *flux       array of fluxes
 *         int nlines          number of lines in array (may not equal number
 *                               of measured fluxes)
 *         float *mean         mean of fluxes (set by function)
 *         float *rms          rms of sample (set by function)
 *         float startday      starting day for calculation
 *         float endday        ending day for calculation
 *         float dayshift      delay to be subtracted from the day member
 *                              of the input flux array.  This is needed
 *                              when calculating, for example, the mean in
 *                              the overlap region between the shifted curves,
 *                              as in the set_mu0 function.
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> failure
 *
 * v25Jun98 CDF, Changed criterion for including a point in the calculation
 *                of the mean from flux->flux > 0 to flux->match != -1 
 *                because the cross-correlation curves, which will now be
 *                passed to this function, have both positive and negative
 *                values.
 */

int calc_mean_dt(Fluxrec *flux, int nlines, float *mean, float *rms,
		 float startday, float endday, float dayshift)
{
  int i;           /* Looping variable */
  int nfluxes=0;   /* Number of measured fluxes */
  float delayday;  /* The day member delayed by dayshift */
  float sum;       /* Running sum */
  Fluxrec *dfptr;  /* Pointer to navigate array of fluxes */

  /*
   * Check on passed startday and endday values.  If they are both set
   *  to zero, then just pass the calculations to the calc_mean function
   *  and return.
   */

  if(startday == 0.0 && endday == 0.0 && dayshift == 0.0) {
    if(calc_mean(flux,nlines,mean,rms)) {
      fprintf(stderr,"ERROR: calc_mean_dt\n");
    }
    else
      return 0;
  }

  /*
   * If startday and endday are NOT set to zero, then proceed with this
   *  function.  Calculate the mean.
   */

  sum = 0.0;

  for(i=0,dfptr=flux; i<nlines; i++,dfptr++) {
    delayday = dfptr->day - dayshift;
    if(dfptr->match != -1 && delayday >= startday && delayday <= endday) {
      sum += dfptr->flux;
      nfluxes++;
    }
  }

  /*
   * Sanity check on number of fluxes
   */

  if(nfluxes == 0) {
    *mean = 0.0;
    *rms = 0.0;
    fprintf(stderr,"ERROR: calc_mean. No input values.\n");
    return 1;
  }
  else if(nfluxes > nlines) {
    fprintf(stderr,"ERROR: calc_mean. More fluxes than expected.\n");
    return 1;
  }

  *mean = sum/nfluxes;

  /*
   * Now find the sample variance
   */

  sum = 0.0;

  if(nfluxes == 1)
    *rms = 0.0;
  else {
    for(i=0,dfptr=flux; i<nlines; i++,dfptr++) {
      delayday = dfptr->day - dayshift;
      if(dfptr->match != -1 && delayday >= startday && delayday <= endday)
	sum += (dfptr->flux - *mean) * (dfptr->flux - *mean);
    }
    *rms = sqrt(sum/(nfluxes-1));
  }

  return 0;
}

/*.......................................................................
 *
 * Function calc_swmean
 *
 * Calculates a statistically weighted sample mean, i.e.
 *
 *         sum(x_i / sig_i^2)
 *  mean = ------------------
 *           sum(1 / sig_i^2)
 *
 * Inputs: Fluxrec *flux       array of fluxes
 *         int nlines          number of lines in array (may not equal number
 *                               of measured fluxes)
 *         float *mean         mean of fluxes (set by function)
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> failure
 *
 */

int calc_swmean(Fluxrec *flux, int nlines, float *mean)
{
  int i;           /* Looping variable */
  int nfluxes=0;   /* Number of measured fluxes */
  float weight;    /* Weight for each point */
  float sum;       /* Running sum of weighted fluxes */
  float wsum;      /* Running sum of weights */
  Fluxrec *fptr;   /* Pointer to navigate array of fluxes */

  /*
   * Calculate the mean
   */

  sum = 0.0;
  wsum = 0.0;
  for(i=0,fptr=flux; i<nlines; i++,fptr++) {
    if(fptr->match != -1) {
      weight = 1.0 / (fptr->err * fptr->err);
      sum += fptr->flux * weight;
      wsum += weight;
      nfluxes++;
    }
  }

  /*
   * Sanity check on number of fluxes
   */

  if(nfluxes == 0) {
    *mean = 0.0;
    return 1;
  }
  else if(nfluxes > nlines) {
    fprintf(stderr,"ERROR: calc_swmean.\n");
    return 1;
  }

  *mean = sum/wsum;
  return 0;
}

/*.......................................................................
 *
 * Function fit_1608
 *
 * Combines the four 1608 light curve into one composite curve using a
 *  number of delays and flux density ratios.  The composite curve has
 *  a function fit to it by calling fit_poly and the reduced chisq is
 *  calculated using chisq_fit.  The parameters giving the best value
 *  for the reduced chisq are stored.
 * For now, only search a small region around the percieved best-fit values.
 *
 * Inputs: Fluxrec *flux[N08]  1608 light curve
 *         int npoints         number of points in each light curve
 *         int nbad[]          number of bad points in each curve
 *         float dt            step size for delay axis
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 */

int fit_1608(Fluxrec *flux[], int npoints, int nbad[], float dt)
{
  int i,j,k,m;            /* Looping variables */
  int no_error=1;         /* Flag set to 0 on error */
  int first=1;            /* Flag for first iteration */
  int nbadtot=0;          /* Total number of bad points in curves */
  int ncompos;            /* Number of points in the composite curve */
  int ndf=2;              /* Number of steps (df) on either side of mean */
  float df=0.02;          /* Flux density scaling step (percentage of mean) */
  float dtmax;            /* Maximum allowed delay */
  float bestchi;          /* Lowest value of reduced chisq */
  float mean[N08];        /* Means of the light curves */
  float m0[N08];          /* Initial guess values for means */
  float d0[N08];          /* Initial guess values for delays */
  Fluxrec bestpars[N08];  /* Delay and flux associated with best chisq */
  Fluxrec *compos=NULL;   /* Composite light curve */
  Fluxrec *fptr;          /* Pointer to navigate light curves */
  int ndeg=10;

  /*
   * First find weighted means of each curve and set initial guess
   *  values for means to these values.
   */

  for(i=0; i<N08; i++) {
    if(calc_swmean(flux[i],npoints,&mean[i])) {
      fprintf(stderr,"ERROR: fit_1608\n");
      return 1;
    }
    m0[i] = mean[i];
  }

  /*
   * Find the maximum allowed delay
   */

  fptr = flux[1] + npoints - 1;
  dtmax = (fptr->day - flux[1]->day)/12.0;

  /*
   * Find total number of bad points in the curves
   */

  for(i=0; i<N08; i++)
    nbadtot += nbad[i];
  printf("\nfit_1608: Total number of bad points = %d\n\n",nbadtot);
  ncompos = 4 * npoints - nbadtot;

  /*
   * Allocate memory for composite curve
   */

  if(!(compos = new_fluxrec(ncompos))) {
    fprintf(stderr,"ERROR: fit_1608\n");
    return 1;
  }

  /*
   * Set the B delay to 0
   */

  d0[1] = 0.0;

  /*
   * Set other delays to their initial guess values
   */

  d0[0] = ALAG;
  d0[2] = CLAG;
  d0[3] = DLAG;

  /*
   * Step through all possible combinations of flux density scalings
   *  given df and ndf.
   * Make the outermost loop on A flux density.
   */

  for(i=-ndf; i<=ndf; i++) {
    mean[0] = m0[0] * (1.0 + i*df);

    /*
     * Loop on B flux density
     */

    for(j=-ndf; j<=ndf; j++) {
      mean[1] = m0[1] * (1.0 + j*df);

      /*
       * Loop on C flux density
       */

      for(k=-ndf; k<=ndf; k++) {
	mean[2] = m0[2] * (1.0 + k*df);

	/*
	 * Loop on D flux density
	 */

	for(m=-ndf; m<=ndf; m++) {
	  mean[2] = m0[2] * (1.0 + m*df);

	  /*
	   * Call fit_delays to loop through possible delays, to do 
	   *  the chisq calculations, and to store the best-fitting 
	   *  model parameters.
	   */

	  if(no_error) {
	    if(fit_delays(flux,npoints,mean,d0,dtmax,4.5,ndeg,compos,
			  ncompos,&first,&bestchi,bestpars))
	      no_error = 0;
	  }
	}
      }
    }
  }

  /*
   * Print out best values
   */

  if(no_error) {
    printf("\nfit_1608: Best reduced chisq = %f found at:\n",bestchi);
    printf("    Delay     Flux Density\n");
    printf("  ---------   ------------\n");
    for(i=0; i<N08; i++) {
      printf("    %5.1f       %7.4f\n",bestpars[i].day,bestpars[i].flux);
    }
  }

  /*
   * Clean up and exit
   */

  compos = del_fluxrec(compos);

  if(no_error)
    return 0;
  else {
    fprintf(stderr,"ERROR: fit_1608.\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function fit_delays
 *
 * This function is called by fit_1608.  It loops over all of the possible
 *  delays and creates a composite curve for each combination of delays.
 *  The composite curve is fit with fit_poly and the reduced chisq of the
 *  fit is compared to the best chisq found.  The values of the parameters
 *  associated with each improved fit are recorded.
 *
 * Inputs: Fluxrec *flux[]     input light curves
 *         int npoints         number of points in input light curves
 *         float mean[]        means of input light curves
 *         float d0[]          initial guesses for delays
 *         float dtmax         maximum delay allowed relative to initial guess
 *         float dt            step size for incrementing delays
 *         int ndeg            degree of function to fit
 *         Fluxrec *compos     composite curve container
 *         int ncompos         number of points in composite curve
 *         int *first          flag set to 1 for first iteration, 0 otherwise
 *         float *bestchi      current minimum reduced chisq
 *         Fluxrec bestpars[]  container for parameters (delays and fluxes) 
 *                              associated with the lowest reduced chisq
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 */

int fit_delays(Fluxrec *flux[], int npoints, float mean[], float d0[], 
	       float dtmax, float dt, int ndeg, Fluxrec *compos, int ncompos, 
	       int *first, float *bestchi, Fluxrec bestpars[])
{
  int i,j;           /* Looping variables */
  int no_error=1;    /* Flag set to 0 on error */
  float chisq;       /* Reduced chisq for fit */
  float delays[N08]; /* Current value of delays */
  float *a=NULL;     /* Container for coefficients of fitting fn. */
  Fluxrec *compptr;  /* Pointer to navigate compos */
  Fluxrec *fptr;     /* Pointer to navigate input light curves */

  /*
   * Set the B delay to 0
   */

  delays[1] = d0[1];

  /*
   * Step through all possible combinations of delays given dtmax and dt.
   * Make the outermost loop on BA delay.
   */

  delays[0] = d0[0] - dtmax;
  while(delays[0] < d0[0] + dtmax+0.5*dt && no_error) {

    /*
     * Initialize BC delay for this iteration
     */

    delays[2] = d0[2] - dtmax;

    /*
     * Loop on BC delay
     */

    while(delays[2] < d0[2]+dtmax+0.5*dt && no_error) {

      /*
       * Initialize BD delay for this iteration
       */

      delays[3] = d0[3] - dtmax;

      /*
       * Loop on BD delay
       */

      while(delays[3] < d0[3]+dtmax+0.5*dt && no_error) {

	/*
	 * Fill the composite curve for this iteration
	 */

	compptr = compos;
	for(i=0; i<N08; i++){
	  for(j=0,fptr=flux[i]; j<npoints; j++,fptr++) {
	    if(fptr->match > -1) {
	      *compptr = *fptr;
	      compptr->day -= delays[i];
	      compptr->flux /= mean[i];
	      compptr->err /= mean[i];
	      compptr++;
	    }
	  }
	}

	/*
	 * Call the fitting function
	 */

	if(!(a = fit_poly(compos,ncompos,11,0)))
	  no_error = 0;

	/*
	 * Now calculate the reduced chisq
	 */

	if(no_error)
	  if((chisq = chisq_fit(compos,ncompos,a,ndeg+1,7,0)) < 0.0)
	    no_error = 0;

	if((chisq < *bestchi || *first) && no_error) {
	  *first = 0;
	  *bestchi = chisq;
	  printf("Reduced chisq = %f ",*bestchi);
	  for(i=0; i<N08; i++) {
	    bestpars[i].day = delays[i];
	    printf("%5.1f ",delays[i]);
	  }
	  for(i=0; i<N08; i++) {
	    bestpars[i].flux = mean[i];
	    printf("%6.3f ",mean[i]);
	  }
	  printf("\n");
	}
	delays[3] += dt;
      }
      delays[2] += dt;
    }
    delays[0] += dt;
  }

  printf("fit_delays: finished loop.\n");

  if(no_error)
    return 0;
  else{
    fprintf(stderr,"ERROR: fit_delays\n");
    return 1;
  }
}

/*.......................................................................
 * 
 * Function fit_poly
 *
 * Fits a function with nparam parameters to a series of points (x,y), using
 *  the Numerical Recipes Singular Value Decomposition method.
 *
 * Inputs: Fluxrec *flux       light curve in Fluxrec format
 *         int npoints         number of points in lightcurve
 *         int nparam          number of parameters to fit
 *         int verbose         flag set to 1 for verbose output
 *
 * Output: *float a            polynomial coefficients
 *
 */

float *fit_poly(Fluxrec *flux, int npoints, int nparam, int verbose)
{
  int i;             /* Looping variable */
  int no_error=1;    /* Flag set to 0 on error */
  float chisq;       /* Chi-sq goodness of fit for polynomial */
  float *x=NULL;     /* Day coordinate of light curve */
  float *y=NULL;     /* Flux coordinate of light curve */
  float *sig=NULL;   /* Error coordinate of light curve */
  float *a=NULL;     /* Coefficients of the polynomial fit */
  float **u=NULL;    /* Workspace for svdfit routine */
  float **v=NULL;    /* Workspace for svdfit routine */
  float *w=NULL;     /* Workspace for svdfit routine */
  Fluxrec *fptr;     /* Pointer to navigate flux */

  /*
   * Allocate memory for 1-d arrays using NR routines to make 
   *  Fortran offsets.
   */

  if(!(x = vector(1,npoints))) {
    fprintf(stderr,"ERROR: fit_poly\n");
    return NULL;
  }

  if(!(y = vector(1,npoints)))
    no_error = 0;

  if(no_error)
    if(!(sig = vector(1,npoints)))
      no_error = 0;

  if(no_error)
    if(!(a = vector(1,nparam)))
      no_error = 0;

  if(no_error)
    if(!(w = vector(1,nparam)))
      no_error = 0;

  /*
   * Allocate memory for 2-d arrays
   */

  if(no_error)
    if(!(u = matrix(1,npoints,1,nparam)))
      no_error = 0;

  if(no_error)
    if(!(v = matrix(1,nparam,1,nparam)))
      no_error = 0;

  /*
   * Assign values to x, y, and sig from flux
   */

  if(no_error) {
    for(i=1,fptr=flux; i<=npoints; i++,fptr++) {
      x[i] = fptr->day;
      y[i] = fptr->flux;
      sig[i] = fptr->err;
    }
  }

  /*
   * Call the SVDFIT function
   */

  if(no_error)
    svdfit(x,y,sig,npoints,a,nparam,u,v,w,&chisq,poly);

  /*
   * Print output
   */

  if(verbose && no_error) {
    printf("\n chisq = %f\n\n",chisq);
    for(i=1; i<=nparam; i++)
      printf("a_%d = %14.10f\n",i-1,a[i]);
  }

  /*
   * Clean up and return
   */

  free_matrix(u,1,npoints,1,nparam);
  free_matrix(v,1,nparam,1,nparam);
  free_vector(x,1,npoints);
  free_vector(y,1,npoints);
  free_vector(sig,1,npoints);
  free_vector(w,1,nparam);

  if(no_error)
    return a+1;
  else {
    fprintf(stderr,"ERROR: fit_poly.\n");
    return NULL;
  }
}

/*.......................................................................
 *
 * Function poly
 *
 * Generates the proper "x"s for polynomial fit of degree ndeg-1.
 * NB:  This routine was taken from Numerical Recipes, and thus uses
 *       the NR Fortran format for arrays.
 *
 * Inputs: float x             value at which polynomial is to be evaluated
 *         float p[]           polynomial "x"s
 *         int ndeg            1 + degree of polynomial
 *
 * Output: none
 *
 */

void poly(float x, float p[], int ndeg)
{
  int j;

  p[1] = 1.0;
  for(j=2; j<=ndeg; j++)
    p[j] = p[j-1]*x;
}

/*.......................................................................
 *
 * Function sinpoly
 *
 * Generates the proper "x"s for a function of sin(nx), where
 *  the maximum value of n is (ndeg-1).  The first "x"s is a constant term.
 * The recursive formula used to calculate sin(nx) comes from Section 5.5
 *  of Numerical Recipes.
 * NB:  This routine is a modification of one taken from Numerical Recipes, 
 *      and thus uses the NR Fortran format for arrays.
 *
 * Inputs: float x             value at which polynomial is to be evaluated
 *         float p[]           function "x"s
 *         int ndeg            n_max + 1
 *
 * Output: none
 *
 */

void sinpoly(float x, float p[], int ndeg)
{
  int j;        /* Looping variable */
  float cosx;   /* Cosine of x */
		   

  p[1] = 1.0;
  p[2] = sin(x/229.0);
  if(ndeg > 2)
    p[3] = sin(2 * x / 229.0);
  if(ndeg > 3) {
    cosx = cos(x/229.0);
    for(j=3; j<=ndeg; j++) {
      p[j] = 2 * cosx * p[j-1] - p[j-2];
    }
  }
}

/*.......................................................................
 *
 * Function chisq_fit
 *
 * Calculates the reduced chisq between a light curve and a polynomial
 *  function that has been fit to the curve.  This functions calculates
 *  the value of the polynomial curve at each measured point in the
 *  light curve.
 *
 * Inputs: Fluxrec *flux       light curve
 *         int npoints         number of points in light curve
 *         float *polyx        temporary holder for "x" values in polynomial
 *         float *a            polynomial coefficients
 *         int ncoeff          number of coefficients
 *         int nhidden         number of parameters used in the fit that
 *                              are "hidden" from fit_poly and thus aren't
 *                              included in ncoeff
 *         int verbose         flag set to 1 for verbose output
 *
 * Output: float chisq         reduced chisq.  Set to a value < 0 on error.
 *
 */

float chisq_fit(Fluxrec *flux, int npoints, float *a, int ncoeff, int nhidden,
		int verbose)
{
  int i,j;            /* Looping variables */
  float chisq=0.0;    /* Reduced chisq */
  float y;            /* Value of polynomial fit at light curve points */
  float *polyx=NULL;  /* Temporary container for  x values */
  Fluxrec *fptr;      /* Pointer to navigate flux */

  /*
   * Allocate memory for container for polynomial values
   */

  if(!(polyx = new_array(ncoeff,1))) {
    fprintf(stderr,"ERROR: chisq_fit.\n");
    return -999.0;
  }

  /*
   * Loop through light curve, calculating value of polynomial at each step
   */

  for(i=0,fptr=flux; i<npoints; i++,fptr++) {

    /*
     * Calculate "x" values
     */

    sinpoly(fptr->day,polyx-1,ncoeff);

    /*
     * Multiply x values by coefficients to get y
     */

    y = 0;
    for(j=0; j<ncoeff; j++)
      y += a[j] * polyx[j];

    /*
     * Add chisq for this point to running sum
     */

    chisq += (fptr->flux - y) * (fptr->flux - y) / (fptr->err * fptr->err);
  }

  /*
   * Convert chisq into reduced chisq
   */

  chisq /= (npoints - ncoeff - nhidden);
  if(verbose)
    printf("\nReduced chisq = %f\n\n",chisq);

  return chisq;
}

/*.......................................................................
 *
 * Function fit_parab
 *
 * Fits a parabola to 3 points (x1,y1), (x2,y2) and (x3,y3).  This can
 *  be done analytically because the equation for a parabola is parameterized
 *  by 3 values (a, b and c).  Thus you have 3 equations in 3 unknowns:
 *
 *     a x1^2 + b x1 + c = 0
 *     a x2^2 + b x2 + c = 0
 *     a x3^2 + b x3 + c = 0.
 *
 *  which gives:
 *
 *     a = ((y1-y3)(x2-x3) - (y2-y3)(x1-x3))/(x1-x2)(x2-x3)(x1-x3)
 *     b = (y2-y3)/(x2-x3) - a(x2+x3)
 *     c = y3 - a*x3^2 - b*x3
 *
 * Inputs: float x1
 *         float y1
 *         float x2
 *         float y2
 *         float x3
 *         float y3
 *         float *a            fitted value of a (set by this function)
 *         float *b            fitted value of b (set by this function)
 *         float *c            fitted value of c (set by this function)
 *         int doprint         flag set to 1 for verbose output
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int fit_parab(float x1, float y1, float x2, float y2, float x3, float y3,
	      float *a, float *b, float *c, int doprint)
{
  /*
   * Error checking
   */

  if(x1 == x2 || x1 == x3 || x2 == x3) {
    fprintf(stderr,"ERROR: fit_parab.  Repeated x values.\n");
    return 1;
  }

  *a = ((y1-y3)*(x2-x3) - (y2-y3)*(x1-x3))/((x1-x2)*(x2-x3)*(x1-x3));
  *b = (y2-y3)/(x2-x3) - *a * (x2+x3);
  *c = y3 - *a * x3 * x3 - *b * x3;

  if(doprint) {
    printf("For the 3 points (%f,%f) (%f,%f) (%f,%f)\n",x1,y1,x2,y2,x3,y3);
    printf(" the parabola has parameters a=%g, b=%f, c=%f\n\n",*a,*b,*c);
    printf("The minimum value occurs at x = %f\n\n",-(*b)/(2 * *a));
  }

  return 0;
}

/*.......................................................................
 *
 * Function daycmp
 *
 * Compares the day fields of two Fluxrec structures and returns 1 if
 *  the first day is greater than the second, 0 if they're equal, and
 *  -1 if the first is less than the second.  This function is called
 *  by qsort.
 *
 * Inputs: void *v1            first structure (cast to void)
 *         void *v2            second structure (cast to void)
 *
 * Output: int (-1,0,1)        as described above
 *
 */

int daycmp(const void *v1, const void *v2)
{
  Fluxrec *f1 = (Fluxrec *) v1;  /* Fluxrec casting of v1 */
  Fluxrec *f2 = (Fluxrec *) v2;  /* Fluxrec casting of v2 */

  /*
   * Do the comparison
   */

  if(f1->day > f2->day)
    return 1;
  else if(f1->day == f2->day)
    return 0;
  else
    return -1;
}

/*.......................................................................
 *
 * Function fluxcmp
 *
 * Compares the flux fields of two Fluxrec structures and returns 1 if
 *  the first day is greater than the second, 0 if they're equal, and
 *  -1 if the first is less than the second.  This function is called
 *  by qsort.
 *
 * Inputs: void *v1            first structure (cast to void)
 *         void *v2            second structure (cast to void)
 *
 * Output: int (-1,0,1)        as described above
 *
 */

int fluxcmp(const void *v1, const void *v2)
{
  Fluxrec *f1 = (Fluxrec *) v1;  /* Fluxrec casting of v1 */
  Fluxrec *f2 = (Fluxrec *) v2;  /* Fluxrec casting of v2 */

  /*
   * Do the comparison
   */

  if(f1->flux > f2->flux)
    return 1;
  else if(f1->flux == f2->flux)
    return 0;
  else
    return -1;
}

/*.......................................................................
 *
 * Function make_compos
 *
 * Creates a composite lightcurve from two to four input 
 *  light curves.
 *
 * Inputs: Fluxrec *flux[]     input light curves
 *         int ncurves         number of curves to be used
 *         int *npoints        number of points in each light curve
 *         int *index          array defining which curves are to be used
 *         float *lag          array of time delays
 *         float *mu           array of scaling factors (scale[i] = 1/mu[i])
 *         int *ncompos        number of points in composite curve (set by
 *                              this function)
 *         int verbose         flag set to 1 for verbose output
 *
 * Output: Fluxrec *compos     composite curve
 *
 * v03Sep2002 CDF, Changed passed variable npoints from a scalar (i.e.,
 *                  all curves have the same number of points) to an
 *                  array (each curve may have a different number of points).
 */

Fluxrec *make_compos(Fluxrec *flux[], int ncurves, int *npoints, 
		     int *index, float *lag, float *mu, int *ncompos,
		     int verbose)
{
  int i,j;                /* Looping variables */
  Fluxrec *compos=NULL;   /* Composite light curve */
  Fluxrec *compptr;       /* Pointer to navigate compos */
  Fluxrec *fptr;          /* Pointer to navigate flux arrays */

  /*
   * Calculate the number of points in the composite curve
   */

  *ncompos = 0;
  for(i=0; i<ncurves; i++)
    *ncompos += npoints[i];
  if(verbose)
    printf("\nmake_compos: %d points in composite curve.\n\n",*ncompos);

  /*
   * Allocate memory for the composite curve
   */

  if(!(compos = new_fluxrec(*ncompos))) {
    fprintf(stderr,"ERROR: make_compos\n");
    return NULL;
  }

  /*
   * Loop through the individual curves and shift by the appropriate
   *  lags to create the composite curve.
   */

  compptr = compos;
  for(i=0; i<ncurves; i++) {
    for(j=0,fptr=flux[index[i]]; j<npoints[i]; j++,fptr++) {
      if(fptr->match > -1) {
	*compptr = *fptr;
	compptr->day -= lag[index[i]];
	compptr->flux /= mu[index[i]];
	compptr->err /= mu[index[i]];
	compptr->match = index[i];
	compptr++;
      }
    }
  }

  /*
   * Sort compos by calling qsort with the comparison function daycmp
   */

  qsort(compos,*ncompos,sizeof(compos[0]),daycmp);

  /*
   * Clean up and return
   */

  return compos;
}
