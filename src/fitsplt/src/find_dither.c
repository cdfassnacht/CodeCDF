/* find_dither.c
 *
 * Description:  Plots two fits files as greyscale plots and then
 *                overlays the catalogs generated from the images.  The
 *                user clicks on an object in each image.  The catalogs
 *                are then used to determine the offsets between the
 *                images.
 *               NB: This program assumes that there is no rotation between
 *                the dithered images.
 *
 * Usage: find_dither fitsfile1 catfile1 fitsfile2 catfile2 ([setup_filename])
 *
 * Revision history
 *  27Aug2006 Chris Fassnacht (CDF),  First working version
 *  30Aug2006 CDF, Added printing of the median shifts to an output file
 *                  and plotting of the median values on the distribution
 *                  plot
 *  31Aug2006 CDF, Fixed output shift file so that it matches the format used
 *                  by the imcombine task in iraf.
 *  2008Jul16 CDF, Small modifications to incorporate new image reading
 *                  and display-opening functions. 
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "libfits.h"
#include "cpgplot.h"
#include "structdef.h"
#include "coords.h"
#include "dataio.h"
#include "fitsim.h"
#include "plotfuncs.h"

/*.......................................................................
 *
 * Function declarations
 *
 */

void get_init_shift(Secat *cat1, int ncat1, Secat *cat2, int ncat2, 
		    Pos *initshift);
Secat find_closest(Pos cpos, Secat *secat, int ncat);
void markcatobj(Secat object);
Secat *find_shifts(Secat *cat1, int ncat1, Secat *cat2, int ncat2,
		   Pos initshift, Setup *setup2, int *nshift, int verbose);
int calc_shift_stats(Secat *shiftcat, int nshift, Pos *shiftmean,
		     Pos *shiftrms, Pos *shiftmedian);
void doubstats(double *data, int ndata, double *mean, double *sig, 
	       double *median);
int doubcmp(const void *v1, const void *v2);

/*.......................................................................
 *
 * Main program
 *
 */

int main(int argc, char *argv[])
{
  int no_error=1;       /* Flag set to 0 on error */
  int usefile=0;        /* Flag set to 1 to use input setup file */
  int ncat1,ncat2;      /* Number of catalog objects */
  int nshift=0;         /* Number of matching objects */
  int another_loop=1;   /* Flag set to 0 when ending a loop */
  char fitsfile1[MAXC]; /* The name of the 1st FITS file to load */
  char catfile1[MAXC];  /* The name of the 1st SExtractor catalog file */
  char fitsfile2[MAXC]; /* The name of the 2nd FITS file to load */
  char catfile2[MAXC];  /* The name of the 2nd SExtractor catalog file */
  char setupfile[MAXC]; /* Name of the optional setup file */
  char line[MAXC];      /* General string to read input */
  Pos initshift;        /* Initial shift between two images */
  Pos shiftmean;        /* Mean shift */
  Pos shiftrms;         /* RMS on shift */
  Pos shiftmed;         /* Median shift */
  Image *image1;        /* The container of the 1st FITS image array */
  Image *image2;        /* The container of the 2nd FITS image array */
  Setup *setup1;        /* Container for the 1st image display info */
  Setup *setup2;        /* Container for the 2nd image display info */
  Secat *catdat1=NULL;  /* SExtractor data for 1st catalog */
  Secat *catdat2=NULL;  /* SExtractor data for 2nd catalog */
  Secat *shiftcat=NULL; /* Catalog of shifts between 1st and 2nd catalogs */

  /*
   * Check the command line
   */

  if(argc < 5) {
    fprintf(stderr,"\nfind_dither fitsfile1 catfile1 fitsfile2 catfile2 ");
    fprintf(stderr,"(setupfile)\n\n");
    return 1;
  }

  /*
   * Get the names of the files.
   */

  strcpy(fitsfile1,argv[1]);
  strcpy(catfile1,argv[2]);
  strcpy(fitsfile2,argv[3]);
  strcpy(catfile2,argv[4]);
  if(argc == 6) {
    usefile = 1;
    strcpy(setupfile,argv[5]);
  }
  else
    sprintf(setupfile,"");

  /*
   * Read the 1st data array.
   */

  if(!(image1 = new_Image(fitsfile1))) {
    fprintf(stderr,"ERROR.  Exiting program. Could not open %s\n\n",
	    fitsfile1);
    return 1;
  }

  if(load_image_data(image1) == ERROR)
    no_error = 0;

  /*
   * Create setup structure and fill with information used to set plotting
   *  parameters.
   */

  if(no_error) 
    if(!(setup1 = fill_setup(image1,setupfile,usefile,PIXEL)))
      no_error = 0;

  /*
   * Loop over displaying image until happy with results
   */

  if(no_error) {
    open_plot_window();
    if(display_image(image1,setup1))
      no_error = 0;
    else
      while(no_error && another_loop) {
	switch(another_loop = setup_menu(image1,setup1)) {
	case 1:
	  if(display_image(image1,setup1))
	    no_error = 0;
	  break;
	case 0:
	  break;
	case -1:
	  no_error = 0;
	  break;
	default:
	  break;
	}
      }
  }

  another_loop = 1;

  /*
   * Read the 2nd data array.
   */

  if(!(image2 = new_Image(fitsfile2))) {
    fprintf(stderr,"ERROR.  Exiting program. Could not open %s\n\n",
	    fitsfile1);
    return 1;
  }

  if(load_image_data(image2) == ERROR)
    no_error = 0;

  /*
   * Create setup structure and fill with information used to set plotting
   *  parameters.
   */

  if(no_error)
    if(!(setup2 = fill_setup(image2,setupfile,usefile,PIXEL)))
      no_error = 0;

  /*
   * Loop over displaying image until happy with results
   */

  if(no_error) {
    open_plot_window();
    if(display_image(image2,setup2))
      no_error = 0;
    else
      while(no_error && another_loop) {
	switch(another_loop = setup_menu(image2,setup2)) {
	case 1:
	  if(display_image(image2,setup2))
	    no_error = 0;
	  break;
	case 0:
	  break;
	case -1:
	  no_error = 0;
	  break;
	default:
	  break;
	}
      }
  }

  another_loop = 1;

  /*
   * Read in the first SExtractor catalog file
   */

  printf("---------------------------------------------------------------\n");
  if(no_error)
    if(!(catdat1 = read_secat(catfile1,'#',&ncat1,6)))
      no_error = 0;

  /*
   * Plot the catalog
   */

  if(no_error) {
    cpgslct(1);
    printf("\nPlotting SExtractor positions for first catalog\n");
    if(plot_secat(setup1,catdat1,ncat1,2,3,2.0))
      no_error = 0;
  }

  /*
   * Read in second SExtractor catalog file
   */

  if(no_error)
    if(!(catdat2 = read_secat(catfile2,'#',&ncat2,6)))
      no_error = 0;

  /*
   * Plot the catalog
   */

  if(no_error) {
    cpgslct(2);
    printf("\nPlotting SExtractor positions for second catalog\n");
    if(plot_secat(setup2,catdat2,ncat2,2,3,2.0))
      no_error = 0;
  }

  /*
   * Get initial shift
   */

  if(no_error) {
    get_init_shift(catdat1,ncat1,catdat2,ncat2,&initshift);
  }

  /*
   * Calculate shifts between first and second catalogs
   */

  if(no_error)
    if(!(shiftcat = find_shifts(catdat1,ncat1,catdat2,ncat2,initshift,
				setup2,&nshift,0)))
      no_error = 0;
    else
      printf("\nCalculated shifts between %d pairs.\n",nshift);

  /*
   * Get statistics on the shifts and plot results
   */

  if(no_error)
    plot_shifts(shiftcat,nshift);

  /*
   * Clean up and exit
   */

  cpgend();
  cpgend();
  image1 = del_Image(image1);
  setup1 = del_setup(setup1);
  catdat1 = del_secat(catdat1);
  image2 = del_Image(image2);
  setup2 = del_setup(setup2);
  catdat2 = del_secat(catdat2);
  shiftcat = del_secat(shiftcat);
  
  if(no_error) {
    printf("\nProgram find_dither.c completed.\n\n");
    return 0;
  }
  else {
    fprintf(stderr,"\nERROR. Exiting find_dither.c\n\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function get_init_shift
 *
 * Gets the initial (x,y) shift between the two catalogs by having the
 *  user click on an object seen in both images.
 *
 * Inputs:
 * 
 */

void get_init_shift(Secat *cat1, int ncat1, Secat *cat2, int ncat2, 
		    Pos *initshift)
{
  int no_error = 1;    /* Flag set to 0 on error */
  float cx,cy;         /* Cursor position */
  char cchar='0';      /* Character returned by cursor command */
  Pos cpos;            /* Cursor position */
  Secat bestmatch1;    /* Closest match to cursor position in image 1 */
  Secat bestmatch2;    /* Closest match to cursor position in image 2 */


  /*
   * Get the desired position in image 1
   */

  printf("---------------------------------------------------------------\n");
  if(no_error) {
    cx = 0.0;
    cy = 0.0;
    printf("\n");
    printf("Find a compact object that is clearly detected in both images\n");
    printf("Click the mouse on the object in the FIRST image\n");
    cpgslct(1);
    if(cpgcurs(&cx,&cy,&cchar)>0) {
      printf("\nMouse clicked at: %7.2f %7.2f\n",cx,cy);
      cpos.x = cx;
      cpos.y = cy;
      bestmatch1 = find_closest(cpos,cat1,ncat1);
      printf("Position of closest match is: %7.2f %7.2f\n",
	     bestmatch1.x,bestmatch1.y);
      markcatobj(bestmatch1);
    }
    else {
      fprintf(stderr,"*** Invalid cursor input ***\n");
      no_error = 0;
    }
  }

  /*
   * Now get position in image 2
   */

  if(no_error) {
    printf("\n Now click the mouse on the object in the SECOND image\n");
    cpgslct(2);
    if(cpgcurs(&cx,&cy,&cchar)>0) {
      printf("\nMouse clicked at: %7.2f %7.2f\n",cx,cy);
      cpos.x = cx;
      cpos.y = cy;
      bestmatch2 = find_closest(cpos,cat2,ncat2);
      printf("Position of closest match is: %7.2f %7.2f\n\n",
	     bestmatch2.x,bestmatch2.y);
      markcatobj(bestmatch2);
    }
    else {
      fprintf(stderr,"*** Invalid cursor input ***\n");
      no_error = 0;
    }
  }

  /*
   * Calculate x,y offsets from selected position in image 1
   */

  printf("---------------------------------------------------------------\n");
  if(no_error) {
    initshift->x = bestmatch2.x - bestmatch1.x;
    initshift->y = bestmatch2.y - bestmatch1.y;
    printf("\nInital estimate of shift between images: %8.2f %8.2f\n",
	   initshift->x,initshift->y);
  }

}

/*.......................................................................
 *
 * Function find_closest
 *
 * Finds the object giving the closest match between the input postion and 
 *  the members of a SExtractor catalog.
 *  The closest matching positions are returned in a Secat structure.
 *
 * Inputs: Pos cpos            input position
 *         Secat *secat        catalog based on fits image (x,y) positions
 *         int ncat            number of members in secat
 *
 * Output: Secat bestmatch     best match
 *
 */

Secat find_closest(Pos cpos, Secat *secat, int ncat)
{
  int i,j;                /* Looping variables */
  double dx,dy,dpos;      /* Offsets between (x,y) pairs */
  double mindpos;         /* Running minimum position offset */
  Secat bestmatch;        /* Container for best-matching xy object */
  Secat *sptr;            /* Pointers to navigate astcat and xycat */

  /*
   * Loop through catalog
   */

  dx = cpos.x - secat->x;
  dy = cpos.y - secat->y;
  bestmatch.dpos = sqrt(dx*dx + dy*dy);
  for(i=0,sptr=secat; i < ncat; i++,sptr++) {
    dx = cpos.x - sptr->x;
    dy = cpos.y - sptr->y;
    if((dpos = sqrt(dx*dx + dy*dy)) < bestmatch.dpos) {
      bestmatch = *sptr;
      bestmatch.dpos = dpos;
    }
  }
  printf("For this position, the best match is ID=%d at",
	 bestmatch.id);
  printf(" mindpos = %6.3f pix\n",bestmatch.dpos);

  /*
   * Exit
   */

  return bestmatch;
}

/*.......................................................................
 *
 * Function find_shifts
 *
 * Given an initial shift, this function first finds matches between a 
 *  shifted version of catalog 1 and catalog 2.  For each match, the
 *  function calculates the shifts between the original catalog 1 position
 *  and the catalog 2 position.  
 *  between the matched pairs.  The
 * The shifts are returned in a Secat array.
 *
 * Inputs: 
 *   Secat *cat1               first catalog
 *   int ncat1                 number of objects in first catalog
 *   Secat *cat2               second catalog
 *   int ncat2                 number of objects in second catalog
 *   Pos initshift             initial guess for shift
 *   Setup *setup2             information about image 2
 *   int *nshift               number of output shifts (set by this function)
 *   int verbose               controls amount of output
 *
 * Output:
 *    Secat *shiftcat          catalog containing shifts between catalog 1
 *                              and catalog 2
 *
 */

Secat *find_shifts(Secat *cat1, int ncat1, Secat *cat2, int ncat2,
		   Pos initshift, Setup *setup2, int *nshift, int verbose)
{
  int i,j;                /* Looping variables */
  double shiftx,shifty;   /* Shifted version of cat1 positions */
  double dx,dy,dpos;      /* Offsets between (x,y) pairs */
  double mindpos;         /* Running minimum position offset */
  double matchlim=3.0;    /* Maximum pixel offset for a "good" match */
  Secat bestmatch;        /* Container for best-matching xy object */
  Secat *shiftcat=NULL;   /* Container for the shifts */
  Secat *ptr1,*ptr2;      /* Pointers to navigate cat1 and cat2 */
  Secat *sptr;            /* Pointer to navigate shiftcat */

  /*
   * Allocate memory for output structure
   */

  if(!(shiftcat = new_secat(ncat2))) {
    fprintf(stderr,"ERROR: find_shifts.\n");
    return NULL;
  }

  /*
   * Loop through first catalog, shifting each of its positions by
   *  initshift.
   */

  sptr = shiftcat;
  for(i=0,ptr1=cat1; i < ncat1; i++,ptr1++) {
    shiftx = ptr1->x + initshift.x;
    shifty = ptr1->y + initshift.y;

    /*
     * Only search for a match if the shifted object actually falls
     *  within the observed image.
     */

    if(shiftx >= 1 && shiftx <= (setup2->xsize - 1) && 
       shifty >= 1 && shifty <= (setup2->ysize - 1)) {
    
      bestmatch.dpos = 1.0 * setup2->xsize;
      for(j=0,ptr2=cat2; j<ncat2; j++,ptr2++) {
	dx = shiftx - ptr2->x;
	dy = shifty - ptr2->y;
	if((dpos = sqrt(dx*dx + dy*dy)) < bestmatch.dpos) {
	  bestmatch = *ptr2;
	  bestmatch.dx = ptr2->x - ptr1->x;
	  bestmatch.dy = ptr2->y - ptr1->y;
	  bestmatch.dpos = dpos;
	}
	if(bestmatch.dpos < matchlim)
	  break;
      }

      /*
       * Only accept best match if it is within the tolerance
       */

      if(bestmatch.dpos < matchlim) {
	*sptr = bestmatch;
	if(verbose) {
	  printf("find_shifts: For %d, mindpos = %6.3f pix, shift = %f %f\n",
		 ptr1->id,sptr->dpos,sptr->dx,sptr->dy);
	}
	sptr++;
	(*nshift)++;
      }
    }
  }

  /*
   * Exit
   */

  return shiftcat;
}

/*.......................................................................
 *
 * Function plot_shifts
 *
 * Calculates statistics of shifts and then plots the shifts.  Also
 *  outputs the median shifts to an output file called tmp.offsets.
 *
 * Inputs:
 *  Secat shiftcat             catalog containing shifts
 *  int nshift                 number of shifts
 *
 * Output:
 *  int (0 or 1)               0 on success, 1 on error
 *
 * v30Aug2006, Added output of median shifts to file and plotting of
 *              median values on plot of distribution.
 */

int plot_shifts(Secat *shiftcat, int nshift)
{
  int i;                     /* Looping variable */
  int no_error=1;            /* Flag set to 0 on error */
  float x1,x2,y1,y2;         /* Limits on plot */
  float *fdx=NULL;           /* float version of x offsets */
  float *fdy=NULL;           /* float version of y offsets */
  float *fxptr,*fyptr;       /* Navigation pointers */
  double *dx=NULL;           /* x offsets */
  double *dy=NULL;           /* y offsets */
  double *xptr,*yptr;        /* Navigation pointers */
  double xmean, xsig, xmed;  /* Statistics on dx */
  double ymean, ysig, ymed;  /* Statistics on dx */
  Secat *sptr;               /* Pointer to navigate shiftcat */
  FILE *ofp=NULL;            /* Output file pointer */

  /*
   * Allocate memory for dx and dy arrays
   */

  if(!(dx = new_doubarray(nshift))) {
    fprintf(stderr,"ERROR: calc_shift_stats\n");
    return 1;
  }
  if(!(dy = new_doubarray(nshift)))
    no_error = 0;
  if(!(fdx = new_array(nshift,1)))
    no_error = 0;
  if(!(fdy = new_array(nshift,1)))
    no_error = 0;

  if(no_error) {

    /*
     * Transfer info to new arrays
     */

    for(i=0,sptr=shiftcat,xptr=dx,yptr=dy,fxptr=fdx,fyptr=fdy; 
	i<nshift; i++,sptr++,xptr++,yptr++,fxptr++,fyptr++) {
      *xptr = sptr->dx;
      *yptr = sptr->dy;
      *fxptr = (float) sptr->dx;
      *fyptr = (float) sptr->dy;
    }

    /*
     * Calculate statistics on dx and dy
     */

    doubstats(dx,nshift,&xmean,&xsig,&xmed);
    doubstats(dy,nshift,&ymean,&ysig,&ymed);

    /*
     * Give output values
     */

    printf("\nStatistics on x shift:\n");
    printf("  mean = %f\n",xmean);
    printf("  rms = %f\n",xsig);
    printf("  median = %f\n",xmed);
    printf("Statistics on y shift:\n");
    printf("  mean = %f\n",ymean);
    printf("  rms = %f\n",ysig);
    printf("  median = %f\n",ymed);

    /*
     * Set the limits and median
     */


    x1 = xmed - 5.0 * xsig;
    x2 = xmed + 5.0 * xsig;
    y1 = ymed - 5.0 * ysig;
    y2 = ymed + 5.0 * ysig;

    /*
     * Plot distribution
     */

    cpgslct(2);
    cpgenv(x1,x2,y1,y2,0,1);
    cpglab("x shift","y shift","Calculated Shifts");
    cpgpt(nshift,fdx,fdy,9);

    /*
     * Plot median
     */

    cpgsci(2);
    cpgslw(5);
    fdx[0] = fdx[1] = xmed;
    fdy[0] = y1;
    fdy[1] = y2;
    cpgline(2,fdx,fdy);
    fdy[0] = fdy[1] = ymed;
    fdx[0] = x1;
    fdx[1] = x2;
    cpgline(2,fdx,fdy);
    cpgsci(1);
    cpgslw(1);
  }

  /*
   * Write median shifts to output file -- NB: for these to be the
   *  proper shifts for an iraf imcombine offsets file, the value
   *  need to be the negative of what the above calculation gives.
   */

  if(!(ofp = open_writefile("tmp.offsets")))
    no_error = 0;
  else
    fprintf(ofp,"%8.2f %8.2f\n",-xmed,-ymed);

  /*
   * Clean up and exit
   */

  dx = del_doubarray(dx);
  dy = del_doubarray(dy);
  fdx = del_array(fdx);
  fdy = del_array(fdy);
  if(ofp)
    fclose(ofp);
  if(no_error)
    return 0;
  else {
    fprintf(stderr,"ERROR: calc_shift_stats\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function doubstats
 *
 * Calculates simple statistics on an array of doubles
 *
 * Inputs:
 *  double *data           input array
 *  int ndata              size of input array
 *  double *mean           mean of data (set by this function)
 *  double *sig            RMS of data (set by this function)
 *  double *median         median of data (set by this function)
 *
 * Output:
 *  (none)
 */

void doubstats(double *data, int ndata, double *mean, double *sig, 
	       double *median)
{
  int i;            /* Looping variable */
  int medindex;     /* Array index used for calculating median */
  double *dptr;     /* Pointer to navigate data */
  double sum=0.0;   /* Running sum */
  double diff;      /* Temporary variable for RMS calculation */

  /*
   * Calculate mean
   */

  for(i=0,dptr=data; i<ndata; i++,dptr++)
    sum += *dptr;
  *mean = sum / ndata;

  /*
   * Calculate RMS
   */

  sum = 0.0;
  for(i=0,dptr=data; i<ndata; i++,dptr++) {
    diff = *dptr - *mean;
    sum += diff * diff;
  }
  *sig = sqrt(sum / (ndata - 1));

  /*
   * Start calculation of median by sorting data
   */

  qsort(data,ndata,sizeof(data[0]),doubcmp);

  /*
   * Find median by taking middle value of array.  The definition of
   *  "middle value" depends on whether the number of elements in
   *  the array (image->ntotal) is even or odd.
   */

  if((ndata/2.0) - (int) (ndata/2.0) == 0) {
    medindex = ndata / 2;
    *median = (data[medindex-1] + data[medindex])/2.0;
  }
  else {
    medindex = (ndata - 1) / 2;
    *median = data[medindex];
  }

}

/*.......................................................................
 *
 * Function doubcmp
 *
 * Compares the two double values and returns 1 if the first is greater 
 *  than the second, 0 if they're equal, and -1 if the first is less 
 *  than the second.  This function is called by qsort.
 *
 * Inputs: void *v1            first structure (cast to void)
 *         void *v2            second structure (cast to void)
 *
 * Output: int (-1,0,1)        as described above
 *
 */

int doubcmp(const void *v1, const void *v2)
{
  double *d1 = (double *) v1;  /* double casting of v1 */
  double *d2 = (double *) v2;  /* double casting of v2 */

  /*
   * Do the comparison
   */

  if(*d1 > *d2)
    return 1;
  else if(*d1 == *d2)
    return 0;
  else
    return -1;
}

/*.......................................................................
 *
 * Function markcatobj
 *
 * Marks position of object contained in the passed secat structure.
 *
 * Input:  Secat object        object to be marked
 *
 * Output: (none)
 *
 */

void markcatobj(Secat object)
{
  cpgsci(2);
  cpgsfs(2);
  cpgslw(3);
  cpgcirc(object.x,object.y,3.0*object.fwhm);
  cpgsci(1);
  cpgslw(1);
}



int MAIN_(void)
{
  return 0;
}
