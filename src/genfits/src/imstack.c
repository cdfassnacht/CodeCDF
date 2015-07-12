/* imstack.c
 *
 * Usage: imstack fitsfile incat format_incat astcat format_astcat outfile
 *
 * Description:  Takes a list of input fits files and median combines them.
 *
 * Revision history:
 *  2010_01_09: First attempt.  Chris Fassnacht (CDF)
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "libfits.h"
#include "structdef.h"
#include "coords.h"
#include "dataio.h"
#include "fitsim.h"

/*.......................................................................
 *
 * Function declarations
 *
 */

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
  int i;
  int no_error=1;          /* Flag set to 0 on error */
  int usefile=0;           /* Flag set to 1 to use input setup file */
  int another_loop=1;      /* Flag set to 0 to stop looping */
  int nfile;               /* Number input images */
  int nflag;               /* Number input flag files */
  int nmatch=0;            /* Number of matches between incat and astcat */
  float **imdat={NULL};    /* The data cube container */
  char fitslist[MAXC];     /* File that contains the list of image files */
  char flaglist[MAXC];     /* File that contains the list of flag files */
  char comb_method[MAXC];  /* Method used to combine the data */
  char outfile[MAXC];      /* Output file */
  char line[MAXC];         /* General string to read input */
  Image **imcube={NULL};   /* The container for the input images */
  Image **flagcube={NULL}; /* The container for the associated flag images */

  /*
   * Check the command line
   */

  printf("%d\n",argc);
  if(argc < 3) {
    fprintf(stderr,"\nUsage: imstack fitslist comb_method [flaglist]\n\n");
    fprintf(stderr,"Inputs:\n");
    fprintf(stderr,
	    "  1. file containing a list of fits files to combine\n");
    fprintf(stderr,"  2. Select 'med', 'ave', or 'sum'\n");
    fprintf(stderr,"  3. [OPTIONAL] file containing a list of fits files to");
    fprintf(stderr," be used as flags\n\n");
    return 1;
  }

  /*
   * Get the names of the files.
   */

  strcpy(fitslist,argv[1]);
  strcpy(comb_method,argv[2]);
  if(argc>2)
    strcpy(fitslist,argv[3]);
#if 0
  /*
   * Load the image and print out basic image information
   */

  if(!(image = new_Image(infits))) {
    fprintf(stderr,"\n\nERROR: Could not read header info from %s\n\n",
	    infits);
    return 1;
  }

  /*
   * Actually load the image data
   */

  if(load_image_data(image) == ERROR)
    no_error = 0;


  /*
   * Clean up and exit
   */

  cpgend();
  image = del_Image(image);
  setup = del_setup(setup);
  incat = del_secat(incat);
  astcat = del_secat(astcat);
  matchcat = del_secat(matchcat);
#endif
  if(no_error) {
    printf("\nProgram imstack.c completed.\n\n");
    return 0;
  }
  else {
    fprintf(stderr,"\nERROR. Exiting imstack.c\n\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function full_display
 *
 * Redisplays image with overlaid astrometric catalog, once the catalog
 *  (x,y) positions have been recalculated.  Returns the catalog of
 *  matches between the astrometric catalog and the input catalog
 *
 * Inputs:
 *  Image *image               fits file and associated info
 *  Setup *setup               plotting parameters
 *  Secat *astcat              astrometric catalog
 *  int nast                   numbers of members in the astrometric catalog
 *  Secat *incat               catalog based on fits image (x,y) positions
 *  int ncat                   number of members in incat
 *  double *matchpix           maximum offset, in pixels, for a "good" match
 *  int *nmatch                number of incat members that match the
 *                              predicted positions (set by this function)
 *
 * Output: Secat *matchcat     catalog containing matched (x,y) positions,
 *                              plus corresponding (RA,Dec) coordinates.
 *
 * Output:
 *  int (0 or 1)               0 ==> success, 1 ==> error
 */

Secat *full_display(Image *image, Setup *setup, Secat *astcat, int nast,
		    Secat *incat, int ncat, double *matchpix, int *nmatch)
{
  char line[MAXC];       /* Generic string for reading input */
  Secat *matchcat=NULL;  /* Catalog of matches */

  /*
   * Using WCS info, convert astrometric data into (x,y) positions
   */

  astcat2xy(astcat,nast,image->wcs);

  /*
   * Display the image
   */

  cpgslct(1);
  if(display_image(image,setup))
    return NULL;

  /*
   * Plot the converted astrometric catalog
   */

  printf("\nPlotting objects in astrometric catalog in green\n");
  if(plot_secat(setup,astcat,nast,2,GREEN,2.0))
    return NULL;

  /*
   * Find the matches between the predicted x,y positions from the
   *  astrometric catalog and the actual measured x,y positions in incat
   */

  if(!(matchcat = find_match(image,astcat,nast,incat,ncat,nmatch,*matchpix,1)))
    return NULL;

  /*
   * Refine maximum offset
   */

  printf("\n");
  printf("Change maximum offset for a good match?\n");
  printf("Enter new value, or just hit return to keep current value ");
  printf("(%5.1f pix = %7.4f arcsec)\n",*matchpix,
	 *matchpix * 3600.0 * image->wcs.pixscale[0]);
  printf("Enter the maximum offset: [%5.1f] ",*matchpix);
  fgets(line,MAXC,stdin);
  if(line[0] != '\n') {
    while(sscanf(line,"%lf",matchpix) != 1) {
      fprintf(stderr,"Invalid number.  Enter value again: ");
      fgets(line,MAXC,stdin);
    }
  }

  return matchcat;
}

/*.......................................................................
 *
 * Function astcat2xy
 *
 * Converts the astrometric catalog RA and Dec positions into (x,y) positions
 *  on the CCD using the wcs information in the CD matrix, etc.
 *
 * Inputs:
 *  Secat *astcat              astrometric catalog
 *  int nast                   number of members in the catalog
 *  WCSinfo wcs                WCS information
 *
 * Output:
 *  int (0 or 1)               0 ==> success, 1 ==> error
 *
 */

int astcat2xy(Secat *astcat, int nast, WCSinfo wcs)
{
  int i;        /* Looping variable */
  Pos offset;   /* Offset in arcsec */
  Secat *aptr;  /* Pointer to navigate astcat */
  WCSinfo tmp;  /* Temporary wcs structure with inverted CD matrix */

  /*
   * Invert the CD matrix
   */

  if(invert_cdmatrix(wcs,&tmp)) {
    fprintf(stderr,"ERROR: astcat2xy");
    return 1;
  }

  /*
   * Loop through astrometric catalog, converting each RA and Dec into
   *  an (x,y) pair, which is stored in the x and y members of the
   *  astcat structure.
   */

  for(i=0,aptr=astcat; i<nast; i++,aptr++) {

    /*
     * First convert the coordinates into RA and Dec offsets, in arcsec
     */

    offset = ddeg2xy(wcs.crval[0],wcs.crval[1],aptr->alpha,aptr->delta);

    /*
     * Convert offsets into degrees
     */

    offset.x /= 3600.0;
    offset.y /= 3600.0;

    /*
     * Now use the inverted CD matrix to convert the offsets into (x,y)
     *  positions on the CCD
     */

    aptr->x = wcs.crpix[0] + offset.x*tmp.cdcol1[0] + offset.y*tmp.cdcol2[0];
    aptr->y = wcs.crpix[1] + offset.x*tmp.cdcol1[1] + offset.y*tmp.cdcol2[1];
#if 0
    printf("%6.2f %6.2f\n",aptr->x,aptr->y);
#endif
  }

  return 1;
}

/*.......................................................................
 *
 * Function find_match
 *
 * Finds matches between the possibly-rotated input astrometric catalog and 
 *  the xy file constructed by running SExtractor on the input image.
 *  The closest matching positions are returned in a Secat array.
 *
 * Inputs: 
 *  Image *image               image plotting information
 *  Secat *astcat              astrometric catalog containing (RA,Dec) and
 *                              _predicted_ (x,y) positions of stars
 *  int nast                   number of members of astcat
 *  Secat *xycat               catalog based on fits image (x,y) positions
 *  int nxy                    number of members in xycat
 *  int *nmatch                number of xycat members that match the
 *                              predicted positions (set by this function)
 *  double matchpix            maximum offset for a "good" match, in pixels
 *  int verbose                flag set to 1 for verbose output
 *
 * Output: Secat *matchcat     catalog containing matched (x,y) positions,
 *                              plus corresponding (RA,Dec) coordinates.
 *
 */

Secat *find_match(Image *image, Secat *astcat, int nast, Secat *xycat, int nxy,
		  int *nmatch, double matchpix, int verbose)
{
  int i,j;                /* Looping variables */
  double dx,dy,dpos;      /* Offsets between (x,y) pairs */
  double mindpos;         /* Running minimum position offset */
  double matchasec;       /* Maximum offset, in arcsec, for a "good" match */
  Secat bestmatch;        /* Container for best-matching xy object */
  Secat *matchcat=NULL;   /* Container for the matches */
  Secat *aptr,*xptr;      /* Pointers to navigate astcat and xycat */
  Secat *mptr;            /* Pointer to navigate matchcat */

  /*
   * Allocate memory for output structure
   */

  if(!(matchcat = new_secat(nast))) {
    fprintf(stderr,"ERROR: find_match.\n");
    return NULL;
  }

  /*
   * Set the maximum pixel offset for a "good" match
   */

  matchasec = matchpix * (3600.0 * image->wcs.pixscale[0]);
  if(verbose) {
    printf("\nfind_match: ");
    printf("Using a %5.1f pixel (%5.2f arcsec) limit for a good match\n",
	   matchpix,matchasec);
  }

  /*
   * Loop through astrometric stars
   */

  *nmatch = 0;
  mptr = matchcat;
  for(i=0,aptr=astcat; i < nast; i++,aptr++) {

    /*
     * Only search for a match if the astrometric object actually falls
     *  within the observed image.
     */

    if(aptr->x >= 1 && aptr->x <= (image->nx - 1) && 
       aptr->y >= 1 && aptr->y <= (image->ny - 1)) {
    
      bestmatch.dpos = 1.0 * image->nx;
      for(j=0,xptr=xycat; j<nxy; j++,xptr++) {
	dx = aptr->x - xptr->x;
	dy = aptr->y - xptr->y;
	if((dpos = sqrt(dx*dx + dy*dy)) < bestmatch.dpos) {
	  bestmatch = *xptr;
	  bestmatch.dx = dx;
	  bestmatch.dy = dy;
	  bestmatch.dpos = dpos;
	}
      }

      /*
       * Only take match if it is closer than limit
       */

      if(bestmatch.dpos < matchpix) {
	if(verbose)
	  printf("find_match: For %s, mindpos = %6.3f pix\n",aptr->name,
		 bestmatch.dpos);
	*mptr = bestmatch;
	mptr->skypos = aptr->skypos;
	strcpy(mptr->name,aptr->name);
	mptr++;
	(*nmatch)++;
      }
    }
  }
  if(verbose) {
    printf("\nfind_match: ");
    printf("Found %d good matches within %5.2f arcsec = %5.2f pix\n",
	   *nmatch,matchasec,matchpix);
  }

  /*
   * Exit
   */

  return matchcat;
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

  bestmatch = *secat;
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

/*.......................................................................
 *
 * Function catxy2dpos
 *
 * Calculates offsets from some central position in (x,y) space.  Converts
 *  the offsets into (r,theta) format and puts results into the dpos and
 *  dpostheta parameters in the input catalog.
 *
 * Inputs: Pos cpos            central position
 *         Secat *secat        catalog based on fits image (x,y) positions
 *         int ncat            number of members in secat
 *
 * Output: (none)
 *
 */

void catxy2dpos(Pos cpos, Secat *secat, int ncat)
{
  int i,j;                /* Looping variables */
  double dx,dy,dpos;      /* Offsets between (x,y) pairs */
  double mindpos;         /* Running minimum position offset */
  Pos tmppos;             /* Temporary container for (x,y) offset */
  Secat *sptr;            /* Pointers to navigate astcat and xycat */

  /*
   * Loop through catalog
   */

  for(i=0,sptr=secat; i < ncat; i++,sptr++) {
    sptr->dx = sptr->x - cpos.x;
    sptr->dy = sptr->y - cpos.y;
    tmppos.x = sptr->dx;
    tmppos.y = sptr->dy;
    xy2rth(tmppos,&sptr->dpos,&sptr->dpostheta,1);
  }
}

/*.......................................................................
 *
 * Function get_click_object
 *
 * Returns the object from the input catalog that is closest to the
 *  (x,y) position in a pgplot display window where the user has clicked the 
 *  mouse.
 *
 * Inputs:
 *  Secat *incat               input catalog
 *  int ncat                   number of catalog members
 *  Secat *bestmatch           best matching object (set by this function)
 *
 * Output:
 *  int (0 or 1)               0 ==> success, 1==> error
 */

int get_click_object(Secat *incat, int ncat, Secat *bestmatch)
{
  float cx,cy;            /* Cursor position */
  char cchar='0';         /* Character returned by cursor command */
  Pos cpos;               /* Cursor position */
  
  if(cpgcurs(&cx,&cy,&cchar)>0) {
    printf("\nMouse clicked at: %7.2f %7.2f\n",cx,cy);
    cpos.x = cx;
    cpos.y = cy;
    *bestmatch = find_closest(cpos,incat,ncat);
    printf("Position of closest match is: %7.2f %7.2f\n",
	   bestmatch->x,bestmatch->y);
    markcatobj(*bestmatch);
  }
  else {
    fprintf(stderr,"*** Invalid cursor input ***\n");
    return 1;
  }
  return 0;
}

/*.......................................................................
 *
 * Function update_shift
 *
 * Gets the shift to apply to the WCS reference pixel in an interactive
 *  manner.  The user will click on an astrometric object and then on the
 *  corresponding object seen in the greyscale image.  Given the initial
 *  guess of the shift, found by calculating delta_x and delta_y between
 *  those two objects, then use the full catalog to refine the shift.
 *  Finally apply the shift to the crpix values in wcs.
 *  user click on an object seen in both images.
 *
 * Inputs: 
 *  Image *image               image data
 *  Secat *astcat              astrometric catalog
 *  int nast                   number of astrometric objects
 *  Secat *xycat               input SExtractor catalog
 *  incat nxy                  number of objects in input catalog
 *  WCSinfo *wcs               WCS info to be updated
 *  int interactive            flag set to 1 for interactive determination
 *                              of the initial shift
 *
 * Output:
 *  int (0 or 1)               0 ==> success, 1 ==> error
 *  
 */

int update_shift(Image *image, Secat *astcat, int nast, Secat *xycat, int nxy,
		 int interactive)
{
  int no_error=1;       /* Flag set to 0 on error */
  int nshift=0;         /* Number of matches found between catalogs */
  char line[MAXC];
  Pos initshift;        /* Initial guess of shift */
  Pos finalshift;       /* Final shift in x,y */
  Secat bestmatchxy;    /* Selected object in xy catalog */
  Secat bestmatchast;   /* Selected object in astrometric catalog */
  Secat *shiftcat=NULL; /* Catalog of offsets */

  /*
   * Using WCS info, convert astrometric data into (x,y) positions
   */

  astcat2xy(astcat,nast,image->wcs);

  /*
   * If interactive flag is set, get the initial shift interactively
   */

  if(interactive) {

    /*
     * Get object in astrometric catalog
     */

    printf("\nClick on an object in the astrometric catalog (GREEN circles)\n");
    if((get_click_object(astcat,nast,&bestmatchast))) {
      fprintf(stderr,"ERROR: update_shift.\n");
      return 1;
    }
    markcatobj(bestmatchast);

    /*
     * Get object in (x,y)  catalog
     */

    printf("\nClick on an object in the SExtractor catalog (YELLOW circles)\n");
    if((get_click_object(xycat,nxy,&bestmatchxy))) {
      fprintf(stderr,"ERROR: update_shift.\n");
      return 1;
    }
    markcatobj(bestmatchxy);

    /*
     * Convert to initial guess for shift
     */

    initshift.x = bestmatchxy.x - bestmatchast.x;
    initshift.y = bestmatchxy.y - bestmatchast.y;
    printf("\nInital estimate of shift between images: %8.2f %8.2f\n",
	   initshift.x,initshift.y);
  }

  /*
   * If not interactive, then set initial guess to (0,0)
   */

  else {
    initshift.x = 0.0;
    initshift.y = 0.0;
  }

  /*
   * Refine shifts by using full catalogs
   */

  if(!(shiftcat = refine_shifts(image,astcat,nast,xycat,nxy,initshift,
				&nshift,15.0,1)))
    no_error = 0;
  else
    printf("\nCalculated shifts between %d pairs.\n",nshift);

  /*
   * Get statistics on the shifts, plot results, and apply them to
   *  crpix
   */

  if(no_error) {
    plot_shifts(shiftcat,nshift,&finalshift);
    printf("\nFinal shifts:\n");
    printf("  Pixels: %8.2f %8.2f\n",finalshift.x,finalshift.y);
    printf("  Arcsec: %8.2f %8.2f\n",finalshift.x*3600.0*image->wcs.pixscale[0],
	   finalshift.y*3600.0*image->wcs.pixscale[1]);
    printf("Applying shifts to CRPIX values.\n");
    image->wcs.crpix[0] += finalshift.x;
    image->wcs.crpix[1] += finalshift.y;
  }

  /*
   * Clean up and exit
   */

  printf("Enter return to continue: ");
  fgets(line,MAXC,stdin);
  shiftcat = del_secat(shiftcat);
  return 0;
}

/*.......................................................................
 *
 * Function refine_shifts
 *
 * Given an initial shift, this function first finds matches between a 
 *  shifted version of catalog 1 and catalog 2.  For each match, the
 *  function calculates the shifts between the original catalog 1 position
 *  and the catalog 2 position.  
 * The shifts are returned in a Secat array.
 *
 * Inputs: 
 *   Image *image              image information    
 *   Secat *cat1               first catalog
 *   int ncat1                 number of objects in first catalog
 *   Secat *cat2               second catalog
 *   int ncat2                 number of objects in second catalog
 *   Pos initshift             initial guess for shift
 *   int *nshift               number of output shifts (set by this function)
 *   double matchlim           maximum offset, in pixels, for a "good" match
 *   int verbose               controls amount of output
 *
 * Output:
 *    Secat *shiftcat          catalog containing shifts between catalog 1
 *                              and catalog 2
 *
 */

Secat *refine_shifts(Image *image, Secat *cat1, int ncat, Secat *cat2, 
		     int ncat2, Pos initshift, int *nshift, double matchlim, 
		     int verbose)
{
  int i,j;                /* Looping variables */
  double shiftx,shifty;   /* Shifted version of cat1 positions */
  double dx,dy,dpos;      /* Offsets between (x,y) pairs */
  Secat bestmatch;        /* Container for best-matching xy object */
  Secat *shiftcat=NULL;   /* Container for the shifts */
  Secat *ptr1,*ptr2;      /* Pointers to navigate cat1 and cat2 */
  Secat *sptr;            /* Pointer to navigate shiftcat */

  /*
   * Allocate memory for output structure
   */

  if(!(shiftcat = new_secat(ncat2))) {
    fprintf(stderr,"ERROR: refine_shifts.\n");
    return NULL;
  }

  /*
   * Loop through first catalog, shifting each of its positions by
   *  initshift.
   */

  sptr = shiftcat;
  for(i=0,ptr1=cat1; i < ncat; i++,ptr1++) {
    shiftx = ptr1->x + initshift.x;
    shifty = ptr1->y + initshift.y;

    /*
     * Only search for a match if the shifted object actually falls
     *  within the observed image.
     */

    if(shiftx >= 1 && shiftx <= (image->nx - 1) && 
       shifty >= 1 && shifty <= (image->ny - 1)) {
    
      bestmatch.dpos = 1.0 * image->nx;
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
	  printf("refine_shifts: For %d, mindpos = %6.3f pix, shift = %f %f\n",
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
 *  Pos *medshift              median shifts (set by this function)
 *
 * Output:
 *  int (0 or 1)               0 on success, 1 on error
 *
 * v30Aug2006, Added output of median shifts to file and plotting of
 *  median values on plot of distribution.
 * v24Jul2007, Now returns the median shifts through the new passed
 *  variable, medshift.
 */

int plot_shifts(Secat *shiftcat, int nshift, Pos *medshift)
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

  medshift->x = xmed;
  medshift->y = ymed;
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


int MAIN_(void)
{
  return 0;
}
