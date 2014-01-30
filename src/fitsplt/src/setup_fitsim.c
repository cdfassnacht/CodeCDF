/* Library: setup_fitsim.c
 *
 * Contains setup section of fits plotting package.  The functions
 *  in this library include those which allocate and free memory for
 *  the setup structure, those which read setup information from the
 *  setup file, and the interactive setup information function.  The
 *  actual functions that query the user for the setup information are
 *  in the get_params.c library.
 *
 * Revision history:
 *  05Dec1999 CDF,  Split out of fitsim.c
 *  v15Feb2001 CDF, Changed default initial display values in setup_interact
 *                   from min --> max to (mean - 10 sigma) --> (mean + 20 sigma)
 *                   where sigma is the rms that comes out of the sigma
 *                   clipping routine.
 *  v22Apr2002 CDF, Modified set_setup_defaults to set more parameters based
 *                   on the image properties in addition to the instrument
 *                   chosen.  This change will take over a lot of the 
 *                   default-setting functionality of the various interactive
 *                   setup tasks.
 *  v04Mar2003 CDF, Changed default initial display values in setup_interact
 *                   to (mean - 2 sigma) --> (mean + 20 sigma)
 *                   where sigma is the rms that comes out of the sigma
 *                   clipping routine.
 *  v12Jul2003 CDF, Added the init_setup function to do initial display before
 *                   displaying to the screen.
 *  v28May2004 CDF, Added a new fill_setup function to do all of the 
 *                   setup-structure creation and filling steps that were
 *                   done in the main calling programs (e.g., fitsplt.c or
 *                   astrom_rot.c)
 *  v29Aug2006 CDF, Changes affecting setting of height of plot labels.  Now
 *                   the user is only queried if the cheight has not yet
 *                   been set, as opposed to being queried every time.  The
 *                   affected functions are new_setup, setup_defaults, and
 *                   setup_interact. 
 *  v09Jul2007 CDF, Small change to setup_interact.
 *  v14Jul2008 CDF, Added a read_slitmask function
 *                  Introduced the setup_menu function, which will eventually
 *                   replace the setup_interact function.
 *  v2008Jul28 CDF, FINALLY added a title keyword to the optional setup file
 *  v2008Aug01 CDF, Made setup_menu more comprehensive.  Modifications to
 *                   setup_defaults to aid this.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "structdef.h"
#include "dataio.h"
#include "libfits.h"
#include "fitswcs.h"
#include "fitsfuncs.h"
#include "fitsim.h"



/*.......................................................................
 *
 * Function new_setup
 *
 * Allocates memory for a Setup structure array.
 *
 * Input:  int size            size of array
 *
 * Output: Setup *newsetup     new setup array
 *
 */

Setup *new_setup(int size)
{
  Setup *newsetup;

  newsetup = (Setup *) malloc((sizeof(Setup) * size));
  if(!newsetup) {
    fprintf(stderr,"Insufficient memory for Setup array.\n");
    return NULL;
  }

  /*
   * Initialize the setup parameters.
   */

  newsetup->border = UNSET;
  newsetup->axislab = ARCSEC;
  newsetup->drawcompass = UNSET;
  newsetup->drawmarker = UNSET;
  newsetup->markx = newsetup->marky = -1.0;
  newsetup->instrument = UNSET;
  newsetup->labheight = -1.0;
  newsetup->labwidth = -1.0;
  newsetup->lo = newsetup->hi = 0.0;
  sprintf(newsetup->title,"");
  newsetup->ncircle = UNSET;
  newsetup->docontour = UNSET;
  newsetup->ncont = 0;
  newsetup->cmul = 1.0;
  newsetup->ndraw = UNSET;
  newsetup->lweight = UNSET;
  newsetup->dolabel = UNSET;
  newsetup->nintlab = UNSET;
  newsetup->intlabs = NULL;
  newsetup->cheight = -1.0;
  newsetup->pixset = UNSET;
  newsetup->pixscale = 0.0;
  newsetup->trans = UNSET;
  newsetup->centset = UNSET;
  newsetup->sizeset = UNSET;
  newsetup->axcentset = UNSET;
  sprintf(newsetup->slitfile,"#");
  newsetup->slitmask = NULL;
  newsetup->plotwfpc = 0;
  newsetup->wfpccent.x = newsetup->wfpccent.y = 0.0;
  newsetup->wfpcpa = 0.0;

  return newsetup;
}

/*.......................................................................
 *
 * Function del_setup
 *
 * Frees memory associated with Setup array
 *
 * Input:  Setup *setup        array to be freed
 *
 * Output: NULL
 *
 */

Setup *del_setup(Setup *setup)
{
  if(setup)
    free(setup);

  return NULL;
}

/*.......................................................................
 *
 * Function setup_help
 *
 * Prints out helpful suggestions for setting up setup file
 *
 * Inputs: none
 *
 * Output: none
 *
 */

void setup_help()
{
  char spaces[MAXC]="                     ";

  printf("\n*******************************************************\n\n");
  printf("USAGE: fitsplt fits_filename [setup_filename]");
  printf(" where the setup file is optional\n\n");
  printf("The setup file consists of lines containing a keyword and then\n");
  printf(" one or more associated values.\n");
  printf("Any line beginning with a # is considered a comment line and");
  printf(" will be ignored.\n\n");

  printf("Currently supported keywords are:\n");

  printf(" border %%d         --");
  printf(" Flag for putting border and external labels on plot.\n");
  printf(" %s1 ==> put border on plot, 0 ==> no border.\n",spaces);

  printf(" circle %%f %%f %%f [%%s] --");
  printf(" Locations, radii, and _optional_ labels\n");
  printf(" %sof circles to be placed on plot.  Parameters are x, y,\n",spaces);
  printf(" %sradius (in pixels) and _optional_ label text.\n",spaces);
  printf(" %sThis requires one line per circle.\n",spaces);

  printf(" cmul %%f           -- Contour multiplication factor.\n");
  printf(" %sUsed, for example, in the case when all contours\n",spaces);
  printf(" %sare multiples of the rms noise in the image\n",spaces);

  printf(" compass %%d        --");
  printf(" Flag for putting compass on plot.\n");
  printf(" %s1 ==> put compass on plot, 0 ==> no compass.\n",spaces);

  printf(" cont %%f           -- Contour level.\n");
  printf(" %sCurrently there must be one line of the form\n",spaces);
  printf(" %s\"cont %%f\" for each desired contour level\n",spaces);

  printf(" draw %%f %%f %%f %%f  --");
  printf(" Start and end positions for internal line segments.\n");
  printf(" %sParameters are x_start, y_start, x_end, y_end.\n",spaces);
  printf(" %sOne line per line segment is required.\n",spaces);

  printf(" frange %%f %%f      -- Range of displayed values.\n");
  printf(" %sParameters are min and max levels\n",spaces);

  printf(" imcent %%f %%f      -- Central pixel location (x and y).\n");

  printf(" imsize %%d %%d      -- Image size (xsize and ysize).\n");

  printf(" lab %%f %%f %%s      --");
  printf(" Internal labels.  Parameters are x, y and label text.\n");
  printf(" %sThis requires one line per label.\n",spaces);

  printf(" labh %%f           -- Character height for internal labels\n");

  printf(" labw %%d           -- Character line width for internal labels.\n");
  printf(" %sMust be an integer in the range 1 - 210.\n",spaces);

  printf(" lweight %%d        -- Line width used in the draw command.\n");
  printf(" %sMust be an integer in the range 1 - 210.\n",spaces);

  printf(" mark %%f %%f        -- Position for crosshair marker.\n");
  printf(" %sParameters are x and y\n",spaces);

  printf(" nocontour         -- No contours are desired (turns off query)\n");

  printf(" nolabel           -- No labels are desired (turns off query)\n");

  printf(" pixscale %%f       -- Pixel scale in arcsec\n");

  printf(" title %%s          -- Title for plot\n");

  printf(" trans %%d          -- Image transfer function type.\n");
  printf(" %sParameter is in the range 0 - 2, for linear (0), log(1)\n",spaces);
  printf(" %sor sqrt(2) transfer function\n",spaces);

  printf("\nFor the above keywords, %%f = floating point, %%d = integer");
  printf(" and %%s = string.\n\n");
}

/*.......................................................................
 *
 * Function fill_setup
 *
 * Creates a new setup structure and fills its parameters with values set
 *  by the input setup file or interactively by the user.
 *
 * Inputs: Image *image         image container
 *         char *setupfile      name of optional setup file
 *         int usefile          set to 1 to use setup file
 *         int autopars         set to 1 to set many parameters automatically
 *
 * Output: Setup *setup         filled setup container
 *
 */

Setup *fill_setup(Image *image, char *setupfile, int usefile, int autopars)
{
  int no_error=1;     /* Flag set to 0 on error */
  Setup *setup=NULL;  /* Container for plotting parameters */

  /*
   * Allocate memory for Setup container
   */

  if(!(setup = new_setup(1))) {
    fprintf(stderr,"ERROR: fill_setup\n");
    return NULL;
  }

  /*
   * Set up image display parameters
   */

  if(no_error)
    if(set_setup_defaults(image,setup))
      no_error = 0;

  if(autopars) 
    if(init_setup(image,setup,autopars))
      no_error = 0;

  if(no_error && usefile)
    if(setup_file(image,setup,setupfile))
      no_error = 0;

  if(no_error)
    setup_interact(image,setup);

  /*
   * Return filled setup container if no errors
   */

  if(no_error)
    return setup;
  else {
    fprintf(stderr,"ERROR: fill_setup\n");
    return del_setup(setup);
  }
}

/*.......................................................................
 *
 * Function init_setup
 *
 * Sets the setup values for an initial image display with no interaction
 *  from the user.
 *
 * Inputs: Image *image        image container
 *         Setup *setup        setup structure to be filled
 *         int autopars        flag controlling whether the axis is requested
 *                              to be in pixels
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 * v29Aug2006, Added a new passed parameter, autopars, that controls whether
 *              the axes are requested to be in pixels.  This will affect
 *              the plotting of the axis labels later.
 */

int init_setup(Image *image, Setup *setup, int autopars)
{
  printf("\ninit_setup: Setting default values.\n");
  if(autopars == PIXEL && image->wcs.validwcs == 0) {
    setup->pixset = TRUE;
    setup->pixscale = 1.0;
    setup->axislab = PIXEL;
  }
  sprintf(setup->title,"%s",image->filename);
  if(image->wcs.validwcs == 0)
    setup->drawcompass = FALSE;
  setup->drawmarker = FALSE;
  setup->border = TRUE;
  setup->cheight = 1.2;
  setup->sizeset = TRUE;
  setup->trans = 0;
  setup->docontour = FALSE;
  setup->dolabel = FALSE;
  setup->ndraw = 0;

  return 0;
}

/*.......................................................................
 *
 * Function set_setup_defaults
 *
 * Sets default values for setup structure depending on the image
 *  properties and on the instrument chosen.
 *
 * Inputs: Image *image        image container
 *         Setup *setup        setup structure to be filled
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 * v28May2004 CDF, Changed default cheight to 1.2
 * v29Aug2006 CDF, Moved cheight default to setup_interact
 * v17Jul2007 CDF, Moved setting of default image size from setup_file to
 *  this function.
 */

int set_setup_defaults(Image *image, Setup *setup)
{
  /*
   * Do the instrument-dependent defaults
   */

  read_instrument(image,setup);
  switch(setup->instrument) {
  case DSS:
    setup->pixset = TRUE;
    if(image->wcs.validwcs == 0)
      setup->pixscale = PIXDSS;
    setup->drawcompass = TRUE;
    setup->drawmarker = AUTO;
    setup->border = TRUE;
    setup->axislab = ARCMIN;
    break;
  case NIRC:
    if(image->wcs.validwcs == 0)
      setup->pixscale = PIXNIRC;
    setup->drawcompass = TRUE;
    break;
  case CCD13:
    if(image->wcs.validwcs == 0)
      setup->pixscale = PIXCCD13;
    setup->drawcompass = TRUE;
    break;
  case LRIS:
    if(image->wcs.validwcs == 0)
      setup->pixscale = PIXLRIS;
    break;
  case RADIO:
    if(image->wcs.validwcs == 0)
      setup->pixscale = PIXDEF;
    setup->drawcompass = TRUE;
    break;
  case WF:
    if(image->wcs.validwcs == 0)
      setup->pixscale = PIXWF;
    setup->drawcompass = UNSET;
    break;
  case PC:
    if(image->wcs.validwcs == 0)
      setup->pixscale = PIXPC;
    setup->drawcompass = UNSET;
    break;
  case NIC1:
    if(image->wcs.validwcs == 0)
      setup->pixscale = PIXNIC1;
    setup->drawcompass = UNSET;
    break;
  case OTHER:
    if(image->wcs.validwcs == 0)
      setup->pixscale = PIXDEF;
    break;
  default:
    if(image->wcs.validwcs == 0)
      setup->pixscale = PIXDEF;
    break;
  }

  /*
   * Set the pixel scale and other defaults if valid wcs has been read from
   *  image
   */

  if(image->wcs.validwcs) {
    setup->pixscale = 3600.0 * image->wcs.pixscale[0];
    setup->drawcompass = TRUE;
  }

  /*
   * Set image size defaults to full image size
   */
#if 0
  if(image->nx > 2048)
    setup->xsize = 2048;
  else
    setup->xsize = image->nx;
  if(image->ny > 2048)
    setup->ysize = 2048;
  else
    setup->ysize = image->ny;
#endif
  setup->xsize = image->nx;
  setup->ysize = image->ny;
  printf("Image default size: %d %d --> %d %d\n",image->nx,image->ny,
	 setup->xsize,setup->ysize);
  setup->centx = (image->nx - 1)/2.0;
  setup->centy = (image->ny - 1)/2.0;
  printf("Image default center: %6.1f %6.1f\n\n",setup->centx,setup->centy);
  setup->markx = setup->centx;
  setup->marky = setup->centy;

  /*
   * Convert center+size into corner positions
   */

  setup->rxstart = setup->centx - (setup->xsize / 2.0);
  setup->rxend   = setup->centx + (setup->xsize / 2.0);
  setup->rystart = setup->centy - (setup->ysize / 2.0);
  setup->ryend   = setup->centy + (setup->ysize / 2.0);

  /*
   * Set frange defaults through a series of steps.  First get image
   *  statistics:
   */

  if(imstat(image)) {
    fprintf(stderr,"ERROR: set_setup_defaults. *****\n");
    fprintf(stderr,"      Could not get image statistics.\n");
  }

  /*
   * Then sigma-clip at 3-sigma level
   */

  if(sigclip(image->data,image->ntotal,&image->pixmean,&image->pixrms,
	     3.0,0.1) < 0) {
    fprintf(stderr,"ERROR: set_setup_defaults. Could not get sigma-clip.\n");
  }

  /*
   * Finally set the defaults in terms of the mean and rms of the clipped
   *  data set.
   */

  setup->lo = image->pixmean - 1.5 *image->pixrms;
  setup->hi = image->pixmean + 10.0 * image->pixrms;

  printf("\nImage parameters:\n");
  printf("--------------------------------\n");
  printf("Minimum:       %9.3f\n",image->datamin);
  printf("Maximum:       %9.3f\n",image->datamax);
  printf("Clipped mean:  %9.3f\n",image->pixmean);
  printf("Clipped rms:   %9.3f\n\n",image->pixrms);
  printf("Default display range: %9.3f %9.3f\n\n",setup->lo,setup->hi);

  return 0;
}

/*.......................................................................
 *
 * Function setup_file
 *
 * Fills in setup container with information from a file.
 *
 * Inputs: Image *image        image data
 *         char *inname        name of file containing info
 *         Setup *setup        setup structure to be filled
 *
 * Output:int (0 or 1)        0 ==> success, 1 ==> error
 *
 * v17Jul2007 CDF, Moved setting of default image size from this function to
 *  setup_defaults
 */

int setup_file(Image *image, Setup *setup, char *inname)
{
  int no_error=1;     /* Flag set to 0 on error */
  float *contptr;     /* Pointer to navigate setup->cont */
  char keyword[MAXC]; /* Keyword at the beginning of each input line */
  char line[MAXC];    /* General string for reading input */
  char *cptr;         /* Pointer to navigate line (for reading title) */
  char *tptr;         /* Pointer to navigate setup->title (for writing title) */
  Pos *startptr;      /* Pointer to starting positions of line segments */
  Pos *endptr;        /* Pointer to starting positions of line segments */
  Labinfo *labptr;    /* Pointer to navigate setup->intlabs */
  Labinfo *circptr;   /* Pointer to navigate setup->circle */
  FILE *ifp=NULL;     /* Input file pointer */

  /*
   * Open input file
   */

  if(!(ifp = open_readfile(inname))) {
    fprintf(stderr,"ERROR: setup_file.\n");
    return 1;
  }

  /*
   * Read input from setup file.  For quantities like the contour levels
   *  and image labels that need to be put in arrays, count the number
   *  of array members first.
   */

  printf("\nReading setup info from file %s\n",inname);
  while(fgets(line,MAXC,ifp) != NULL && no_error) {
    if(line[0] >= 32 && line[0] != '#') {
      switch(read_setup_line(line,keyword)) {
      case SETUPERR:
	no_error = 0;
	break;
      case AXISCENT:
	if(sscanf(line,"%s %f %f",keyword,&setup->axcentx,&setup->axcenty) 
	   != 3) {
	  fprintf(stderr,"ERROR: setup_file.  Bad input for axiscent.\n");
	  fprintf(stderr," Axes will be centered in center of image.\n");
	  setup->axcentset = FALSE;
	}
	else
	  setup->axcentset = TRUE;
	break;
      case BORDER:
	if(sscanf(line,"%s %d",keyword,&setup->border) != 2 &&
	   (setup->border < 0 || setup->border > 1)) {
	  fprintf(stderr,"ERROR: setup_file.  Bad input for border\n");
	  fprintf(stderr,"Taking border = 1\n");
	  setup->border = TRUE;
	}
	break;
      case CIRCLE:
	if(setup->ncircle == UNSET)
	  setup->ncircle = 1;
	else
	  setup->ncircle++;
	break;
      case CMUL:
	if(sscanf(line,"%s %f",keyword,&setup->cmul) != 2) {
	  fprintf(stderr,"ERROR: setup_file.  Bad input for cmul\n");
	  fprintf(stderr,"Taking cmul = 1.0\n");
	  setup->cmul = 1.0;
	}
	break;
      case COMPASS:
	if(sscanf(line,"%s %d",keyword,&setup->drawcompass) != 2 &&
	   (setup->drawcompass < 0 || setup->drawcompass > 1)) {
	  fprintf(stderr,"ERROR: setup_file.  Bad input for compass\n");
	  fprintf(stderr,"Taking compass = 0\n");
	  setup->drawcompass = FALSE;
	}
	break;
      case CONT:
	setup->docontour = TRUE;
	if(setup->ncont == UNSET)
	  setup->ncont = 1;
	else
	  setup->ncont++;
	break;
      case DRAW:
	if(setup->ndraw == UNSET)
	  setup->ndraw = 1;
	else
	  setup->ndraw++;
	break;
      case FRANGE:
	if(sscanf(line,"%s %f %f",keyword,&setup->lo,&setup->hi) != 3) {
	  fprintf(stderr,"ERROR: setup_file.  Bad input for frange\n");
	  setup->lo = 0.0;
	  setup->hi = 0.0;
	}
	break;
      case IMCENT:
	if(sscanf(line,"%s %f %f",keyword,&setup->centx,&setup->centy) != 3) {
	  fprintf(stderr,"ERROR: setup_file.  Bad input for imcent.\n");
	  fprintf(stderr," Setting imcent to half of image size.\n");
	  setup->centx = setup->xsize / 2.0;
	  setup->centy = setup->ysize / 2.0;
	  setup->centset = TRUE;
	}
	else if(setup->centx < 1.0 || setup->centx > 1.0*image->nx) {
	  fprintf(stderr,"ERROR: setup_file. Bad x value of imcent.\n");
	  fprintf(stderr," Setting x value to half of xsize.\n");
	  setup->centx = setup->xsize / 2.0;
	  setup->centset = TRUE;
	}
	else if(setup->centy < 1.0 || setup->centy > 1.0*image->ny) {
	  fprintf(stderr,"ERROR: setup_file. Bad y value of imcent.\n");
	  fprintf(stderr," Setting y value to half of ysize.\n");
	  setup->centy = setup->ysize / 2.0;
	  setup->centset = TRUE;
	}
	else
	  setup->centset = TRUE;
	break;
      case IMSIZE:
	if(sscanf(line,"%s %d %d",keyword,&setup->xsize,&setup->ysize) != 3) {
	  fprintf(stderr,"ERROR: setup_file.  Bad input for imsize.\n");
	  fprintf(stderr," Setting imsize to full image size.\n");
	  setup->xsize = image->nx;
	  setup->ysize = image->ny;
	  setup->sizeset = TRUE;
	}
	else {
	  sprintf(line,"%d %d",setup->xsize,setup->ysize);
	  if(check_imsize(image,setup,line,0)) {
	    fprintf(stderr,"ERROR: setup_file.\n");
	    no_error = 0;
	  }
	  else
	    setup->sizeset = TRUE;
	}
	break;
      case LAB:
	setup->dolabel = TRUE;
	if(setup->nintlab == UNSET)
	  setup->nintlab = 1;
	else
	  setup->nintlab++;
	break;
      case LABH:
	if(sscanf(line,"%s %f",keyword,&setup->labheight) != 2) {
	  fprintf(stderr,"ERROR: setup_file.  Bad input for labh\n");
	  fprintf(stderr,"Taking labh = 2.0\n");
	  setup->labheight = 2.0;
	}
	break;
      case LABW:
	if(sscanf(line,"%s %d",keyword,&setup->labwidth) != 2) {
	  fprintf(stderr,"ERROR: setup_file.  Bad input for labw\n");
	  fprintf(stderr,"Taking labw = 10\n");
	  setup->labwidth = 10;
	}
	break;
      case LWEIGHT:
	if(sscanf(line,"%s %d",keyword,&setup->lweight) != 2) {
	  fprintf(stderr,"ERROR: setup_file.  Bad input for lweight\n");
	  fprintf(stderr,"Taking lweight = 1\n");
	  setup->lweight = 1;
	}
	break;
      case MARK:
	if(sscanf(line,"%s %f %f",keyword,&setup->markx,&setup->marky) != 3) {
	  fprintf(stderr,"ERROR: setup_file.  Bad input -- marker position\n");
	  fprintf(stderr,"No marker will be drawn.\n");
	  setup->drawmarker = FALSE;
	}
	else {
	  setup->drawmarker = TRUE;
	}
	break;
      case NOCONTOUR:
	setup->docontour = FALSE;
	break;
      case NOLABEL:
	setup->dolabel = FALSE;
	break;
      case PIXSCALE:
	if(sscanf(line,"%s %f",keyword,&setup->pixscale) != 2) {
	  fprintf(stderr,"ERROR: setup_file.  Bad input for pixel scale.\n");
	  setup->pixscale = 0.0;
	  setup->pixset = UNSET;
	}
	else {
	  setup->pixset = TRUE;
	}
	break;
      case PLOTWFPC:
	if(sscanf(line,"%s %d",keyword,&setup->plotwfpc) != 2 &&
	   setup->plotwfpc < 0) {
	  fprintf(stderr,"ERROR: setup_file.  Bad input for plotwfpc\n");
	  fprintf(stderr,"Taking plotwfpc = 0 (no overlay)\n");
	  setup->plotwfpc = 0;
	}
	break;
      case SLITMASK:
	if(sscanf(line,"%s %s",keyword,setup->slitfile) != 2) {
	  fprintf(stderr,"ERROR: setup_file. Bad input for slitmask name\n");
	  sprintf(setup->slitfile,"#");
	}
	break;
      case TITLE:
	if((cptr = strstr(line,"title")) == NULL && 
	   (cptr = strstr(line,"TITLE")) == NULL)
	  fprintf(stderr,"ERROR: Title keyword not found in title line\n");
	else {
	  cptr += 6;
	  while(*cptr == ' ')
	    cptr++;
	  tptr = setup->title;
	  while(*cptr != '\n' && *cptr != '\0') {
	    *tptr = *cptr;
	    cptr++;
	    tptr++;
	  }
	  *tptr = '\0';
	}
	break;
      case TRANS:
	if(sscanf(line,"%s %d",keyword,&setup->trans) != 2 &&
	   (setup->trans < 0 || setup->trans > 2)) {
	  fprintf(stderr,"ERROR: setup_file.  Bad input for trans\n");
	  fprintf(stderr,"Taking trans = 0 (linear scaling)\n");
	  setup->trans = 0;
	}
	break;
      default:
	printf("***WARNING: Not yet taking file info for keyword %s.\n",
	       keyword);
      }
    }
  }

  /*
   * Allocate arrays
   */

  if(no_error && setup->ncont != UNSET) {
    if(!(setup->cont = new_array(setup->ncont,1)))
      no_error = 0;
    else {
      printf("  setup_file: Found %d contours in setup file.\n",
	     setup->ncont);
      contptr = setup->cont;
    }
  }
  if(no_error && setup->nintlab != UNSET) {
    if(!(setup->intlabs = new_labinfo(setup->nintlab)))
      no_error = 0;
    else {
      printf("  setup_file: Found %d internal labels in setup file.\n",
	     setup->nintlab);
      labptr = setup->intlabs;
    }
  }
  if(no_error && setup->ncircle != UNSET) {
    if(!(setup->circle = new_labinfo(setup->ncircle)))
      no_error = 0;
    else {
      printf("  setup_file: Found %d circle positions in setup file.\n",
	     setup->ncircle);
      circptr = setup->circle;
    }
  }
  if(no_error && setup->ndraw != UNSET) {
    if(!(setup->drawstart = new_pos(setup->ndraw,1)))
      no_error = 0;
    else 
      startptr = setup->drawstart;
    if(!(setup->drawend = new_pos(setup->ndraw,1)))
      no_error = 0;
    else 
      endptr = setup->drawend;
  }

  /*
   * Fill the arrays and get other info
   */

  rewind(ifp);
  
  while(fgets(line,MAXC,ifp) != NULL && no_error) {
    if(line[0] > 32 && line[0] != '#') {
      switch(read_setup_line(line,keyword)) {
      case SETUPERR:
	no_error = 0;
	break;
      case CIRCLE:
	switch(sscanf(line,"%s %f %f %f %s",keyword,&circptr->x,&circptr->y,
		      &circptr->size,circptr->text)) {
	case 4:
	  strcpy(circptr->text,"");
	  circptr++;
	  break;
	case 5:
	  circptr++;
	  break;
	default:
	  no_error = 0;
	}
	break;
      case CONT:
	if(sscanf(line,"%s %f",keyword,contptr) != 2)
	  no_error = 0;
	else {
	  *contptr *= setup->cmul;
	  contptr++;
	}
	break;
      case DRAW:
	if(sscanf(line,"%s %lf %lf %lf %lf",keyword,&startptr->x,&startptr->y,
		  &endptr->x,&endptr->y) != 5)
	  no_error = 0;
	else {
	  startptr++;
	  endptr++;
	}
	break;
      case LAB:
	if(sscanf(line,"%s %f %f %s",keyword,&labptr->x,&labptr->y,
		  labptr->text) != 4)
	  no_error = 0;
	else {
	  labptr++;
	}
	break;
      default:
	break;
      }
    }
  }

  /*
   * Clean up and exit
   */

  if(ifp)
    fclose(ifp);

  if(no_error)
    return 0;
  else {
    fprintf(stderr,"ERROR: setup_file\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function read_setup_line
 *
 * Reads one line of the setup file and checks the keyword on the line.
 * Returns an integer value identifying the keyword value.
 *
 * Input:  char *line          input line
 *         char *keyword       keyword read from input file (modified by this
 *                              function)
 *
 * Output: int keyval          keyword value (see enumeration in fitsim.h)
 *                             SETUPERR for error
 *
 */

int read_setup_line(char *line, char *keyword)
{
  /*
   * Check to make sure that there is a keyword on the line.
   */

  if(sscanf(line,"%s",keyword) != 1) {
    fprintf(stderr,"ERROR: read_setup_line\n");
    keyword = NULL;
    return SETUPERR;
  }

  /*
   * If there is a keyword, read it and return the appropriate value
   */

  if(strcmp(keyword,"axiscent") == 0 || 
     strcmp(keyword,"AXISCENT") == 0)
    return AXISCENT;
  if(strcmp(keyword,"border") == 0 || 
     strcmp(keyword,"BORDER") == 0)
    return BORDER;
  if(strcmp(keyword,"circle") == 0 || 
     strcmp(keyword,"CIRCLE") == 0)
    return CIRCLE;
  if(strcmp(keyword,"cmul") == 0 || 
     strcmp(keyword,"CMUL") == 0)
    return CMUL;
  if(strcmp(keyword,"compass") == 0 || 
     strcmp(keyword,"COMPASS") == 0)
    return COMPASS;
  if(strcmp(keyword,"cont") == 0 || 
     strcmp(keyword,"CONT") == 0)
    return CONT;
  if(strcmp(keyword,"draw") == 0 || 
     strcmp(keyword,"DRAW") == 0)
    return DRAW;
  if(strcmp(keyword,"frange") == 0 ||
     strcmp(keyword,"FRANGE") == 0)
    return FRANGE;
  if(strcmp(keyword,"imcent") == 0 ||
     strcmp(keyword,"IMCENT") == 0)
    return IMCENT;
  if(strcmp(keyword,"imsize") == 0 ||
     strcmp(keyword,"IMSIZE") == 0)
    return IMSIZE;
  if(strcmp(keyword,"lab") == 0 || 
     strcmp(keyword,"LAB") == 0)
    return LAB;
  if(strcmp(keyword,"labh") == 0 || 
     strcmp(keyword,"LABH") == 0)
    return LABH;
  if(strcmp(keyword,"labw") == 0 || 
     strcmp(keyword,"LABW") == 0)
    return LABW;
  if(strcmp(keyword,"lweight") == 0 || 
     strcmp(keyword,"LWEIGHT") == 0)
    return LWEIGHT;
  if(strcmp(keyword,"mark") == 0 || 
     strcmp(keyword,"MARK") == 0)
    return MARK;
  if(strcmp(keyword,"nocontour") == 0 || 
     strcmp(keyword,"NOCONTOUR") == 0)
    return NOCONTOUR;
  if(strcmp(keyword,"nolabel") == 0 || 
     strcmp(keyword,"NOLABEL") == 0)
    return NOLABEL;
  if(strcmp(keyword,"pixscale") == 0 || 
     strcmp(keyword,"PIXSCALE") == 0)
    return PIXSCALE;
  if(strcmp(keyword,"plotwfpc") == 0 || 
     strcmp(keyword,"PLOTWFPC") == 0)
    return PLOTWFPC;
  if(strcmp(keyword,"slitmask") == 0 || 
     strcmp(keyword,"SLITMASK") == 0)
    return SLITMASK;
  if(strcmp(keyword,"title") == 0 || 
     strcmp(keyword,"TITLE") == 0)
    return TITLE;
  if(strcmp(keyword,"trans") == 0 || 
     strcmp(keyword,"TRANS") == 0)
    return TRANS;

  /*
   * If none of the above checks have been satisfied, return the
   *  default value
   */

  return DEFAULT;
}

/*.......................................................................
 *
 * Function setup_menu
 *
 * An interactive menu to change setup parameters.  This will eventually
 *  take the place of setup_interact.
 *
 * Inputs: 
 *  Image *image               fits file
 *  Setup *setup               plotting info
 *
 * Output: 
 *  int result                 -1 ==> error
 *                              0 ==> stop display loop
 *                             +1 ==> continue display loop
 */

int setup_menu(Image *image, Setup *setup)
{
  char line[MAXC];    /* Generic string for reading input */

  printf("\nAdjust display menu\n");
  printf("-------------------\n");
  printf(" Enter one of the following codes:\n");
  printf("   f - Change display range\n");
  printf("   p - Change pixel scale of image\n");
  printf("   t - Change title\n");
  printf("   z - Change display center and size (crude zoom)\n");
  printf("       (%d %d centered at %f %f)\n",setup->xsize,setup->ysize,
	 setup->centx,setup->centy);
  printf("   + - Toggle compass display on/off\n");
  printf("   x - Add/move the cross-hair marker\n");
  printf("   q - Quit out of display adjustment and continue with program\n");
  printf("\nEnter option: [q] ");
  fgets(line,MAXC,stdin);
  if(line[0] == '\n')
    return 0;
  else {
    switch(line[0]) {
    case 'f': case 'F':
      get_frange(image,setup);
      break;
    case 'p': case 'P':
      get_pixscale(image,setup);
      break;
    case 't': case 'T':
      get_title(image,setup);
      break;
    case 'z': case 'Z':
      get_imsize(image,setup);
      break;
    case '+':
      toggle_compass(setup);
      break;
    case 'x':
      setup->drawmarker = TRUE;
      get_marker_info(setup);
      break;
    case 'q': case 'Q':
      return 0;
      break;
    default:
      fprintf(stderr,"\nsetup_menu: Not a valid option.\n");
      break;
    }
  }

  return 1;
}

/*.......................................................................
 *
 * Function setup_interact
 *
 * Sets up image display parameters in an interactive manner.  The setup
 *  parameters are put into a Setup structure.
 *
 * Inputs: Image *image        image container
 *         Setup *setup        setup structure to be filled
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 * v25Mar1999 CDF, Added an estimation of pixel size for radio images via
 *                  a modification of read_file_pixscale.
 * v09Jul2007 CDF, Changed call to get_frange to reflect change in get_params.c
 */

int setup_interact(Image *image, Setup *setup)
{
  char line[MAXC];     /* General string variable for getting input */


  /*
   * Get pixel scale if not already set
   */

  if(image->wcs.validwcs == 0) {
    if(get_pixscale(image, setup)) {
      fprintf(stderr,"ERROR: setup_interact\n");
      return 1;
    }
  }
#if 0
  /*
   * Get image size
   */

  if(get_imsize(image,setup)) {
    fprintf(stderr,"ERROR: setup_interact\n");
    return 1;
  }

  /*
   * Mark an object?  Note that the mark position is set to be the
   *  requested center pixel in the get_imsize function
   */

  if(setup->drawmarker != TRUE)
    get_marker_info(setup);

  /*
   * Get display range.  Note that display range for DSS plots is
   *  set automatically.
   */

  printf("\nSetting limits for image display.\n");
  /*
   * Get overall image statistics.
   */

  if(setup->lo != setup->hi) {
    printf("  *** Using display values set in setup file ***\n");
    printf("      Min = %f\n",setup->lo);
    printf("      Max = %f\n",setup->hi);
  }
  else {
    setup->lo = image->pixmean - 2.0 *image->pixrms;
    setup->hi = image->pixmean + 20.0 * image->pixrms;
    printf("   Initial guess for range: %9.3f %9.3f\n",setup->lo,
	   setup->hi);
    if(get_frange(image,setup)) {
      fprintf(stderr,"ERROR: setup_interact\n");
      return 1;
    }
  }

  /*
   * Get image title
   */

  get_title(image,setup);
#endif
  /*
   * Get image scaling information
   */

  if(setup->instrument != DSS && setup->trans == UNSET)
    if(get_transfunc(setup))
      fprintf(stderr,"ERROR: setup_interact\n");
  
  /*
   * Get contour info
   */

  if(setup->instrument != DSS && setup->docontour != FALSE)
    if(get_contours(setup))
      fprintf(stderr,"ERROR: setup_interact\n");

  /*
   * Get internal label info
   */

  if(setup->instrument != DSS && setup->dolabel != FALSE)
    if(get_intlab(image,setup))
      fprintf(stderr,"ERROR: setup_interact\n");

  /*
   * Get image scaling information
   */

  if(setup->nintlab > 0 && setup->instrument != DSS)
    if(get_labscale(setup))
      fprintf(stderr,"ERROR: setup_interact\n");
  
  /*
   * Draw a compass?
   */

  if(setup->drawcompass == UNSET) {
    printf("Draw a compass? [y] ");
    fgets(line,MAXC,stdin);
    if(line[0] == 'n' || line[0] == 'N')
      setup->drawcompass = FALSE;
    else
      setup->drawcompass = TRUE;
  }

  /*
   * Put border and external labels on the plot?
   */

  if(setup->instrument != DSS)
    if(get_border_info(image,setup))
      fprintf(stderr,"ERROR: setup_interact\n");

  /*
   * Get the character height for plot labels if not set.  The default
   *  value is 1.2, but the user can change this.
   */

  if(setup->instrument != DSS && setup->cheight<0.0) {
    setup->cheight = 1.2;
    printf("\nEnter character height for plot labels: [%4.1f]  ",
	   setup->cheight);
    fgets(line,MAXC,stdin);
    if(line[0] != '\n') {
      while(sscanf(line,"%f",&setup->cheight) != 1 || setup->cheight < 0.0) {
	fprintf(stderr,"ERROR: Bad input format.  Enter character height: ");
	fgets(line,MAXC,stdin);
      }
    }
  }

  /*
   * If a WFPC2 overlay is desired, get the overlay parameters.
   */

  if(setup->plotwfpc)
    if(get_wfpc_overlay(setup))
      fprintf(stderr,"ERROR: setup_interact\n");
    
  return 0;
}

/*.......................................................................
 *
 * Function read_slitmask
 *
 * Reads slitmask information from the input file into the Secat structure
 *  that will store the information.
 *
 * The input file has probably been produced by the find_prevslits.sh
 *  script.  It should have the following format:
 *
 *     name ra_hr ra_min ra_sec dec_deg dec_amin dec_asec mag priority
 *      %s   %d    %d     %lf    %d       %d      %lf     %f    %d
 *
 *  However, this can easily be approximated by the format=9 for read_secat,
 *   which has the format: name,hr,min,sec,deg,amin,asec,mag,magerr
 *   and then just assign the magerr values to priority.
 *
 * Inputs: 
 *  Image *image               fits file
 *  Setup *setup               plotting info
 *
 * Output: 
 *  int (SUCCESS or ERROR)
 *
 */

int read_slitmask(Image *image, Setup *setup)
{
  /*
   * Read slit positions into secat
   */

  if(!(setup->slitmask = read_secat(setup->slitfile,'#',&setup->nmask,9))) {
    fprintf(stderr,"ERROR: read_slitmask\n");
    return ERROR;
  }

  /*
   * Convert the (RA,Dec) coordinates in the secat container into (x,y)
   *  positions on the chip.
   */

  if(radec_to_ccdxy(image->fits,setup->slitmask,setup->nmask) == ERROR) {
    fprintf(stderr,"ERROR: read_slitmask");
    setup->slitmask = del_secat(setup->slitmask);
    return ERROR;
  }
  else
    return SUCCESS;
}

