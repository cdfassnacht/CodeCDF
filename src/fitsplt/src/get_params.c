/*
 * get_params.c
 *
 * This is a collection of functions that get the parameters that are
 *  put into the Setup container for fitsplt.c
 *
 * Function list:
 *  get_instrument
 *  read_instrument
 *  get_object
 *  read_object
 *  get_pixscale
 *  read_dss_pixscale
 *  read_file_pixscale
 *  get_imsize
 *  check_imsize
 *  get_title
 *  toggle_compass
 *  get_marker_info
 *  get_frange
 *  get_epoch
 *  get_zooms
 *  get_transfunc
 *  get_intlab
 *  get_labscale
 *  get_border_info
 *  get_wfpc_overlay
 *
 * 09Feb2000 CDF,  Split off from fitsim.c.  For early history, see fitsim.c
 *                  documentation.
 * v22Apr2002 CDF, Added the new read_instrument, read_object, and
 *                  read_dss_pixscale functions.
 *                 Small modifications to read_file_pixscale and get_pixscale.
 * v26Aug2002 CDF, Modifications to get_pixscale
 * v20Apr2004 CDF, Modifications to read_file_pixscale and get_border_info.
 * v28May2004 CDF, Modifications to get_imsize, get_marker_info.
 * v11Apr2005 CDF, Modification to read_instrument.
 * V06Jul2007 CDF, Moved much of the functionality of read_file_pixscale
 *                  into the new read_wcs_info in fitswcs.c
 * v09Jul2007 CDF, Changed get_frange passed parameters to eliminate
 *                  duplication
 * v17Jul2007 CDF, Updated get_imsize to prepare for more flexible usage
 * v2008Jul28 CDF, Added a get_title function.
 * v2008Jul31 CDF, Added a toggle_compass function and updated get_marker_info.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "structdef.h"
#include "libfits.h"
#include "fitswcs.h"
#include "fitsim.h"

/*.......................................................................
 *
 * Function get_instrument
 *
 * Sets observing instrument interactively
 *
 * Inputs: Setup *setup        setup info
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 * v06Nov99 CDF, Added "quit" option.
 *               Improved user interface.
 * v22Apr02 CDF, Split off initialization into new read_instrument function.
 */

int get_instrument(Image *image, Setup *setup)
{
  char line[MAXC];   /* General string variable for getting input */

  /*
   * Initialize to default value
   */

  read_instrument(image,setup);

  /*
   * Get instrument value
   */

  printf("Observation made with:\n");
  printf("  %d. DSS/POSS2\n",DSS);
  printf("  %d. NIRC\n",NIRC);
  printf("  %d. CCD13\n",CCD13);
  printf("  %d. LRIS\n",LRIS);
  printf("  %d. VLA/VLBA\n",RADIO);
  printf("  %d. WFPC2-WF\n",WF);
  printf("  %d. WFPC2-PC\n",PC);
  printf("  %d. NICMOS-NIC1\n",NIC1);
  printf("  %d. Other\n",OTHER);
  printf(" %d. Quit\n",QUIT);
  printf("Enter choice [%d]: ",setup->instrument);
  fgets(line,MAXC,stdin);
  if(line[0] != '\n') {
    while(sscanf(line,"%d",&setup->instrument) != 1 || 
	  setup->instrument < QUIT || setup->instrument > NIC1) {
      fprintf(stderr,"ERROR: invalid input.  Enter choice again:  ");
      fgets(line,MAXC,stdin);
    }
  }

  return 0;
}

/*.......................................................................
 *
 * Function read_instrument
 *
 * Reads instrument values from the image header and transfers that
 *  information into the setup container.
 *
 * Inputs: Image *image        image container
 *         Setup *setup        setup container
 *
 * Output: (none)
 *
 * v11Apr2005 CDF, Added a NIRC option
 */

void read_instrument(Image *image, Setup *setup)
{
  /*
   * Use instrument or telescope value to set setup->instrument
   */

  if(strcmp(image->telescope,"P60") == 0 && 
     strcmp(image->instrument,"CCD 13") == 0)
    setup->instrument = CCD13;
  else if(strcmp(image->telescope,"VLA") == 0 ||
	  strcmp(image->telescope,"VLBA") == 0)
    setup->instrument = RADIO;
  else if(strncmp(image->telescope,"Keck",4) == 0) {
    if(strncmp(image->instrument,"NIRC",4) == 0)
      setup->instrument = NIRC;
    else
      setup->instrument = LRIS;
  }
  else if(strncmp(image->instrument,"WFPC2",5) == 0)
    setup->instrument = WF;
  else
    setup->instrument = OTHER;
}

/*.......................................................................
 *
 * Function get_object
 *
 * Fills in the source variable in the setup structure.  If the observations
 *  come from the DSS, the default name is the name of the fits file
 *  minus the ".fits" extension.  Otherwise, the function attempts to read
 *  the "object" header card and will take that name if successful.
 *
 * Inputs: Image *image        image container
 *         Setup *setup        setup structure to be filled
 *         char *filename      FITS file name
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 * v22Apr02 CDF, Moved reading of object name into new read_object function
 */

int get_object(Image *image, Setup *setup, char *filename)
{
  int nchars;      /* Number of characters in filename */
  int nsufchar=0;  /* Number of characters in file suffix */
  char *cptr;      /* Pointer used to navigate filename */
  char line[MAXC]; /* General string for reading input */

  /*
   * For the DSS the source name will NOT be in the object header card AND
   *  for speed the user interface is not used.  So take the fits filename
   *  (assumed to be [sourcename].fits) and extract the sourcename from it.
   * If the filename does not end with ".fits", take the whole filename as
   *  the source name.
   */

  if(setup->instrument == DSS) {
    nchars = strlen(filename);
    if(nchars) {
      cptr = filename + nchars - 1;
      while(*cptr != '.' && cptr != filename) {
	nsufchar++;
	cptr--;
      }
      if(*cptr == '.' && *(cptr+1) == 'f' && *(cptr+2) == 'i' &&
	 *(cptr+3) == 't') {
	strncpy(setup->source,filename,nchars-nsufchar-1);
      }
      else {
	strcpy(setup->source,filename);
      }
      printf("\nUsing source name %s\n",setup->source);
    }
    else
      strcpy(setup->source,"");
  }

  /*
   * If not a DSS file, read the source name from the object header
   *  card in the image.
   */

  else {
    read_object(image,setup);

    /*
     * Now, get the source name, taking the FITS object header as the
     *  default
     */

    printf("Enter source name [%s]:  ",setup->source);
    gets(line);
    if(strcmp(line,"") != 0) {
      strcpy(setup->source,line);
    }
  }

  return 0;
}

/*.......................................................................
 *
 * Function read_object
 *
 * Reads in the object name from the "object" header card.
 *
 * Inputs: Image *image        image container
 *         Setup *setup        setup structure to be filled
 *
 * Output: (none)
 *
 */

void read_object(Image *image, Setup *setup)
{
  Phdu *phdu;      /* Primary data-header unit in the fits file */

  /*
   * Set phdu to point at the appropriate place
   */

  phdu = (Phdu *) image->fits->hdu;

  /*
   * Read the name.
   */

  if(phdu->object) {
    printf("\nSource name in file is %s\n",phdu->object);
    strcpy(setup->source,phdu->object);
  }
  else {
    printf("\nNo source name found\n");
    strcpy(setup->source,"");
  }
}

/*.......................................................................
 *
 * Function get_pixscale
 *
 * Function to interactively set the pixel scale for the plot.
 *
 *
 * Inputs: Image *image         image to be displayed
 *         Setup *setup         image display information
 *                               to the desired units/pixel.
 *
 * Output: int (0 or 1)         0 ==> success, 1 ==> error
 *
 * v22Apr02 CDF, Split off DSS pixel scale calculation into new
 *                read_dss_pixscale
 *               Made the user choice of pixel scale easier.
 * v26Aug03 CDF, Took out call to read_file_pixscale since that is already
 *                done in image_info function
 */

int get_pixscale(Image *image, Setup *setup)
{
  int scalechoice=1;   /* Counting variable */
  int foundflag=0;     /* Flag set to 1 if able to read pixscale from header */
  char line[MAXC];     /* General string variable for getting input */

  /*
   * First check to see if the pixel scale has been set in the setup
   *  file.  If it has, then exit the program.  However, if the
   *  image is from a DSS or POSS2 plate, check the header cards
   *  for information.
   */

  if(setup->pixset == TRUE) {
    if(setup->instrument == DSS)
      read_dss_pixscale(image,setup);
    else
      printf("\nUsing pixel scale from setup file: %5.3f arcsec/pix.\n",
	     setup->pixscale);
    return 0;
  }

  /*
   * If the pixel scale has not been set, first try to read it from the
   *  image header.
   * This step is now done in the image_info function, and the results
   *  are stored in the image->pixscale parameter.  If the 
   *  read_file_pixscale call in image_info was successful, the
   *  image->pixscale parameter will have a non-zero value.
   */

  printf("Choose pixel scale:\n");
  printf(" ------------------------------------------------------\n");
  printf("  1. Default pixel scale             -- ");
  printf("%6.4f arcsec/pix\n",setup->pixscale);
  if(image->pixscale > 0) {
    printf("  2. Pixel scale from image header   -- ");
    printf("%6.4f arcsec/pix\n",image->pixscale);
    scalechoice = 2;
  }
  else {
    printf("  2. (not available)");
    printf(" -- no pixel scale in image headers!\n");
  }
  printf("  3. Another value.\n");
  printf("  ---------------------------------------------\n");
  printf("  Enter choice: [%d] ",scalechoice);
  fgets(line,MAXC,stdin);
  if(line[0] != '\n') {
    while(sscanf(line,"%d",&scalechoice) != 1) {
      fprintf(stderr,"  ERROR: bad input format.  Enter choice:  ");
      fgets(line,MAXC,stdin);
    }
  }

  /*
   * Set the pixel scale depending on choice above.  NB: By choosing
   *  option 1 above (i.e., scalechoice == 1), the user wants the default
   *  pixel scale, in which case setup->pixscale already has the correct
   *  value.
   */

  if(scalechoice == 2)
    setup->pixscale = image->pixscale;
  else if(scalechoice != 1) {
    printf("\n User input -- Enter pixel scale in arcsec/pix: ");
    fgets(line,MAXC,stdin);
    if(line[0] != '\n') {
      while(sscanf(line,"%f",&setup->pixscale) != 1) {
	fprintf(stderr,"  ERROR: bad input format.  Enter pixel scale:  ");
	fgets(line,MAXC,stdin);
      }
    }
  }

  return 0;
}

/*.......................................................................
 *
 * Function read_dss_pixscale
 *
 * Reads, if possible, pixel scale from image headers for DSS or
 *  POSS2 images.
 *
 * Inputs: Image *image         image to be displayed
 *         Setup *setup         image display information
 *
 * Output: (none)
 */

void read_dss_pixscale(Image *image, Setup *setup)
{
  /*
   * Read pixel scale from CDELT1 and CDELT2 headers
   */

  printf("\nSetting pixel scale...\n");
  if(read_file_pixscale(image,60.0)) {
    setup->pixscale = PIXDSS;
    fprintf(stderr,"***** WARNING: Could not read pixel scale from ");
    fprintf(stderr," file.  Using default scale of %5.3f arcmin/pix\n",
	    setup->pixscale);
  }
  else {
    if(image->pixscale == 0.0) {
      fprintf(stderr,"\nERROR.  File pixel scale is 0.0!");
      fprintf(stderr,
	      " Using default pixel scale of %5.3f arcmin/pix.\n",
	      setup->pixscale);
    }
    else {
      printf(" Using pixel scale %5.3f arcmin/pix\n",image->pixscale);
      setup->pixscale = image->pixscale;
    }
  }
}

/*.......................................................................
 *
 * Function read_file_pixscale
 *
 * Reads the pixel scale of the fits file contained in image.  This
 *  information should be stored in the CDELT1 and CDELT2 header cards.  
 * The value in the header cards, which should be in degrees, is converted
 *  to the desired scale (arcsec or arcmin) via the passed parameter "scale".
 * If the CDELTn header cards do not exist, then check for CDn_m headers.
 * If those header cards do not exist, then the function will return the 
 *  error value of 1.
 *
 * Inputs: Image *image         image to be displayed
 *         float scale          factor to change scale from degrees/pixel
 *                               to the desired units/pixel.
 *
 * Output: int (0 or 1)         0 ==> success, 1 ==> error
 *
 * v22Apr2002 CDF, Moved output pixel scale into the image container
 *                  rather than having it in its separate variable.
 * v10May2003 CDF, Added calculation of rotation from CDm_n values, 
 *                  a la IDL processing.
 * v26Aug2003 CDF, Added checks for invalid CTYPEs.
 * v20Apr2004 CDF, Added rotation and cdmatx parameters to Image structure.
 * v02Sep2006 CDF, Fixed bug where wrong pixel scale was returned in case
 *                  of unequal pixel scales on the two axes.  What was being
 *                  returned was the average of the two pixel scales (an OK
 *                  approach, but needs to get fixed in a future version),
 *                  but in degrees/pix rather than in arcsec/pix.
 *                 Took out some of the warning messages if CDELT* headers
 *                  weren't there, since CD matrices are much more common
 *                  these days.
 * v06Jul2007 CDF, Pushed most of functionality into new read_wcs_info
 *                  function in fitswcs.c
 */

int read_file_pixscale(Image *image, float scale)
{
  int i;              /* Looping variable */
  int no_error=1;     /* Flag set to 0 on error */

  /*
   * Read in WCS information
   */

  if(read_wcs_info(image->fits,&image->wcs))
    no_error = 0;
  else {
    image->rotation = image->wcs.rotation;
  }

  /*
   * Print out values in arcsec
   */

  if(no_error) {
    printf("  read_file_pixscale: Using CDELT1 = %7.4f arcsec/pix\n",
	   image->wcs.cdelt[0]*3600.0);
    printf("  read_file_pixscale: Using CDELT2 = %7.4f arcsec/pix\n",
	   image->wcs.cdelt[1]*3600.0);
    printf("  read_file_pixscale: Using rotation = %7.2f\n",
	   (180.0/PI) * image->rotation);
  }

  /*
   * Check the values
   */

  if(no_error) {
    if(image->wcs.cdelt[0] + image->wcs.cdelt[1] > 
       fabs(0.001 * image->wcs.cdelt[1])) {
      image->pixscale = 
	scale * (fabs(image->wcs.cdelt[0]) + fabs(image->wcs.cdelt[1]))/2.0;
      fprintf(stderr,"   WARNING: read_file_pixscale. Unequal pixel scales.\n");
      fprintf(stderr,
	      "     Setting pixscale to the average of the x and y values.\n");
    }
    else {
      image->pixscale = scale * fabs(image->wcs.cdelt[0]);
    }
  }

  /*
   * Clean up and return
   */

  if(no_error)
    return 0;
  else
    return 1;
}

/*.......................................................................
 *
 * Function get_imsize
 *
 * Queries the user for the image dimensions to display and does some
 *  checking of the validity of the choices.
 *
 * Inputs: Image *image         image to be displayed
 *         Setup *setup         image display information
 *
 * Output: int (0 or 1)         0 ==> success, 1 ==> error
 *
 * v26May1999 CDF, Fixed bug in which marker position was not being correctly
 *  set in some situations.
 * v28May2004 CDF, Moved setting of rxstart, etc. from display_image to
 *  this function.
 * v09Jul2007 CDF, Changed assumptions - now assumes that default values 
 *  are whatever setup->xsize, setup->ysize, setup->centx, and setup->centy
 *  are set to before calling this function.  This gives the function more
 *  versatility.
 */

int get_imsize(Image *image, Setup *setup)
{
  int no_error=1;   /* Flag set to 0 on error */
  char line[MAXC];  /* General string for reading input */

  /*
   * Get and check the image size (x and y), unless this is a DSS finding
   *  chart, in which case we just take the whole image.
   */

  printf("\nSetting image size:\n");
  printf("-------------------\n");
    printf(" The full image has dimensions %d x %d.\n",image->nx,image->ny);

  if(setup->instrument == DSS) {
    printf(" POSS2/DSS: The plotted image will have dimensions %d x %d.\n",
	   setup->xsize,setup->ysize);
  }
  else {
    printf(" Enter \"0 0\" to get the full image size.\n");
    printf(" Enter the size of the image in the form \"xsize ysize\"  ");
    printf("[%d %d]:  ",setup->xsize,setup->ysize);
    fgets(line,MAXC,stdin);
    if(line[0] != '\n') {
      if(check_imsize(image,setup,line,1)) {
	fprintf(stderr," ERROR: get_imsize\n");
	return 1;
      }
    }
  }

  /*
   * If the full image is requested, put the appropriate info into
   *  the setup container.
   */

  if(setup->xsize == image->nx && setup->ysize == image->ny) {
    setup->centx = (image->nx - 1)/2.0;
    setup->centy = (image->ny - 1)/2.0;
  }

  /*
   * If less than the full image is requested, find out where the
   *  image subsection is to be centered.
   */

  else {
    printf("\n Enter the center pixel position:  [%6.2f %6.2f] ",
	   setup->centx,setup->centy);
    fgets(line,MAXC,stdin);
    if(line[0] != '\n') {
      while(sscanf(line,"%f %f",&setup->centx,&setup->centy) != 2 ||
	    setup->centx < 1.0 || setup->centx > 1.0*image->nx || 
	    setup->centy < 1.0 || setup->centy > 1.0*image->ny) {
	fprintf(stderr,"  ERROR. Not a valid position. ");
	printf("Enter the position again [x y]:  ");
	fgets(line,MAXC,stdin);
      }
    }
  }

  /*
   * Convert center+size into corner positions
   */

  setup->rxstart = setup->centx - (setup->xsize / 2.0);
  setup->rxend = setup->centx + (setup->xsize / 2.0);
  setup->rystart = setup->centy - (setup->ysize / 2.0);
  setup->ryend= setup->centy + (setup->ysize / 2.0);

  return 0;
}

/*.......................................................................
 *
 * Function check_imsize
 *
 * Checks to see if requested image size is reasonable.
 *
 * Inputs: Image *image         image to be displayed
 *         Setup *setup         image display information
 *         char *line           line containing xsize and ysize
 *         int interactive      flag set to 1 if interactive
 *
 * Output: int (0 or 1)         0 ==> success, 1 ==> error
 *
 */

int check_imsize(Image *image, Setup *setup, char *line, int interactive)
{
  if(interactive) {
    while(sscanf(line,"%d %d",&setup->xsize,&setup->ysize) != 2) { 
      fprintf(stderr,"ERROR: check_imsize.  Bad input values\n");
      fprintf(stderr,"Enter xsize and ysize again:  ");
      fgets(line,MAXC,stdin);
    }
    if(setup->xsize == 0 || setup->ysize == 0) {
      printf(" Image size = 0 ==> full image size\n");
      setup->xsize = image->nx;
      setup->ysize = image->ny;
    }
  }
  else {
    if(setup->xsize < 1 || setup->xsize > image->nx) {
      fprintf(stderr,"ERROR: check_imsize. Bad x value of imsize.\n");
      fprintf(stderr," Setting x value to full image value.\n");
      setup->centx = image->nx;
    }
    else if(setup->ysize < 1 || setup->ysize > image->ny) {
      fprintf(stderr,"ERROR: check_imsize. Bad y value of imsize.\n");
      fprintf(stderr," Setting y value to full image value.\n");
      setup->centy = image->ny;
    }
  }

  return 0;
}

/*.......................................................................
 *
 * Function get_title
 *
 * Gets the title for the displayed plot
 *
 * Inputs:
 *  Image *image            container for image data and header info
 *  Setup *setup            setup container (includes title)
 *
 * Output:
 *  (none)
 */

void get_title(Image *image, Setup *setup)
{
  int not_at_end=1;  /* Flag set to 0 at end of string */
  float epoch=0.0;   /* Epoch of observations */
  char line[MAXC];   /* General string for reading input */
  char *cptr;        /* Pointer to navigate line */

  /*
   * Automatically set the title for DSS images
   */

  if(setup->instrument == DSS) {
    if(get_epoch(image,&epoch) == 0 && epoch > 1960.0)
      sprintf(setup->title,"%s Finding Chart (POSS2)",setup->source);
    else {
      epoch = 0.0;
      sprintf(setup->title,"%s Finding Chart (DSS)",setup->source);
    }
  }
  else {
    if(strcmp(setup->title,"") == 0)
      printf("\nEnter the title for the plot [(none)]:  ");
    else
      printf("\nEnter the title for the plot [%s]:  ",setup->title);
    fgets(line,MAXC,stdin);
    if(line[0] != '\n') {
      cptr = line;
      while(not_at_end) {
	if(*cptr == '\0')
	  not_at_end = 0;
	else if(*cptr == '\n') {
	  *cptr = '\0';
	  not_at_end = 0;
	}
	else
	  cptr++;
      }
      strcpy(setup->title,line);
    }
  }
}

/*.......................................................................
 *
 * Function toggle_compass
 *
 * Switches the state of the setup->drawcompass parameter
 *
 * Inputs:
 *  Setup *setup            setup container (includes drawcompass parameter)
 *
 * Output:
 *  (none)
 */

void toggle_compass(Setup *setup)
{
  switch(setup->drawcompass) {
  case TRUE:
    printf("toggle_compass: Turning compass display OFF\n");
    setup->drawcompass = FALSE;
    break;
  case FALSE: case UNSET:
    printf("toggle_compass: Turning compass display ON\n");
    setup->drawcompass = TRUE;
    break;
  default:
    printf("toggle_compass: drawcompass is in an unrecognized state.\n");
    printf("Turning compass display OFF\n");
    setup->drawcompass = FALSE;
    break;
  }
}

/*.......................................................................
 *
 * Function get_marker_info
 *
 * Gets information about whether to plot a cross-hair marker and, if so,
 *  where to plot it.
 *
 * Inputs: Setup *setup         image display information
 *
 * Output: (none)
 *
 * v28May2004 CDF, Fixed bug where setup->drawmarker==FALSE was not 
 *                  checked.
 * v31Jul2008 CDF, Made more compatible with new setup_menu function in
 *                  setup_fitsim.c
 */

void get_marker_info(Setup *setup)
{
  char line[MAXC];   /* General string variable for reading input */

  if(setup->drawmarker == FALSE)
    return;

  /*
   * For the POSS and POSS2 plates, a marker is automatically drawn
   *  in the center of the images, so the default values are fine.
   * Otherwise, query if a marker is desired
   */

  if(setup->drawmarker == UNSET) {
    setup->markx = setup->centx;
    setup->marky = setup->centy;
    printf("\nDraw a cross-hair marker? [n] ");
    fgets(line,MAXC,stdin);
    if(line[0] != 'y' && line[0] != 'Y') {
      setup->drawmarker = FALSE;
      return;
    }
  }

  /*
   * Get marker position
   */

  printf("  Enter position for marker in format \"x y\": [%f %f] ",
	 setup->markx,setup->marky);
  fgets(line,MAXC,stdin);
  if(line[0] != '\n') {
    while(sscanf(line,"%f %f",&setup->markx,&setup->marky) != 2) {
      fprintf(stderr,"ERROR: Bad input format.  Enter position again:  ");
      fgets(line,MAXC,stdin);
    }
  }

  setup->drawmarker = TRUE;
  printf("\nCross-hair marker will be drawn at %f %f.\n",setup->markx,
	 setup->marky);
}

/*.......................................................................
 *
 * Function get_frange
 *
 * Gets the range of flux values to be plotted.
 *
 * Input: Image *image         container for image data and header info
 *        Setup *setup         image display information (modified by this
 *                              function).
 *        float f_low          current lower limit for flux display
 *        float f_high         current upper limit for flux display
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 * v09Jul2007 CDF, Eliminated f_low and f_high passed parameters since
 *  they easily could be included in the passed setup structure.
 */

int get_frange(Image *image, Setup *setup)
{
  char line[MAXC];    /* General string for reading input */


  /*
   * If this is a POSS or POSS2 plate, set the display to the
   *  pre-determined values.
   */

  if(setup->instrument == DSS) {
    setup->lo = image->datamin+((image->datamax-image->datamin)/20.0);
    setup->hi = image->datamax-((image->datamax-image->datamin)/20.0);
  }

  /*
   * Otherwise give the user a chance to change from the current
   *  display limits set by f_low and f_high.
   */

  else {
    printf("\n--------------------------------------------------\n\n");
    printf("Set Display Limits:\n");
    printf("  Enter low and high limits to display.");
    printf("  Just hit return to use current limits.\n\n");
    printf("  Current display limits:  %f  %f\n",
	   setup->lo,setup->hi);
    printf("  Enter limits in format \"low high\": [use current limits] ");
    fgets(line,MAXC,stdin);

    /*
     * Default -- just use f_low and f_high as display limits.
     */

    if(line[0] != '\n') {
      while(sscanf(line,"%f %f",&setup->lo,&setup->hi) != 2) {
	fprintf(stderr,"  ERROR: bad input format.  Enter low and high ");
	fprintf(stderr,"limits again:  ");
	fgets(line,MAXC,stdin);
      }
    }
  }

  return 0;
}

/*.......................................................................
 *
 * Function get_epoch
 *
 * Reads the EPOCH header card, if it exists, and returns the value.
 *
 * Input: Image *image         container for image data and header info
 *        float *epoch         value of EPOCH card (set by this function)
 *
 * Output: int (0 or 1)        0 ==> EPOCH card found, 1 ==> not found
 *
 */

int get_epoch(Image *image, float *epoch)
{
  Fitkey key;         /* Used for reading FITS header */

  /*
   * Check for the EPOCH header card
   */

  if(get_key(image->fits, image->fits->hdu, "EPOCH", DAT_DBL, LOOP_SEEK, &key))
    return 1;
  else {
    *epoch = KEYDBL(key);
    printf("\nget_epoch: Epoch = %6.1f\n",*epoch);
    return 0;
  }
}

/*.......................................................................
 *
 * Function get_zooms
 *
 * Gets the zoom factors for each axis (defaults are 1.0 for each).
 *
 * Input:  Setup *setup         image display information
 *
 * Output: int (0 or 1)         0 ==> success, 1 ==> error
 *
 */

int get_zooms(Setup *setup)
{
  float zoom1=1.0;   /* Zoom factors for x axis */
  float zoom2=1.0;   /* Zoom factors for y axis */
  char line[MAXC];   /* General string for reading input */

  printf("\nEnter zoom factors for x and y axes: [%3.1f %3.1f] ",zoom1,zoom2);
  gets(line);
  if(strcmp(line,"") != 0) {
    while(sscanf(line,"%f %f",&zoom1,&zoom2) != 1 && zoom1 <= 0.0 && 
	  zoom2 <= 0.0) {
      fprintf(stderr,"ERROR: bad input.  Enter zoom factors again:  ");
      gets(line);
    }
  }

  setup->zoom[0] = zoom1;
  setup->zoom[1] = zoom2;

  return 0;
}


/*.......................................................................
 *
 * Function get_transfunc
 *
 * Gets image transfer function for the greyscale plotting.  The
 *  three possible choices are linear, log or sqare-root.
 *
 * Input:  Setup *setup         image display information
 *
 * Output: int (0 or 1)         0 ==> success, 1 ==> error
 *
 */

int get_transfunc(Setup *setup)
{
  int tchoice=0;     /* Transfer function value (set to default) */
  char line[MAXC];   /* General string for reading input */

  printf("\nImage transfer function choices\n");
  printf("  0. Linear\n");
  printf("  1. Logarithmic\n");
  printf("  2. Square root\n");
  printf("---------------------------------\n");
  printf("Enter choice [%d]:  ",tchoice);
  fgets(line,MAXC,stdin);
  if(line[0] != '\n') {
    while(sscanf(line,"%d",&tchoice) != 1 && tchoice < 0 && tchoice > 3) {
      fprintf(stderr,"ERROR: bad input.  Enter choice:  ");
      fgets(line,MAXC,stdin);
    }
  }

  setup->trans = tchoice;

  return 0;
}

/*.......................................................................
 *
 * Function get_intlab
 *
 * Gets interal labels, if needed
 *
 * Inputs: Image *image        image to be displayed
 *         Setup *setup        image display information
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 * v12Oct99 CDF, Added some error checking and cleaned up code a bit.
 * v06Nov99 CDF, Added "reset" variable to allow this function to
 *                be used for resetting the label positions as well
 *                as setting them in the first place.
 */

int get_intlab(Image *image, Setup *setup)
{
  int i;              /* Looping variable */
  int reset=0;        /* Flag set to 1 if re-setting the label positions */
  char line[MAXC];    /* General string for reading input */
  Labinfo *labptr;    /* Pointer for navigating setup->intlabs */

  /*
   * See if this function is being used to reset the label position
   *  or if it is being called for the first time.  An easy way to do
   *  this is to see if setup->intlabs is NULL or not.  If it is
   *  _not_ NULL, then it means that this function has been called
   *  before (or that the label values have been set in the input
   *  setup file, and thus that this function is being used to reset
   *  the values.
   */
   
  if(setup->intlabs) {
     reset = 1;
  }

  /*
   * Get number of desired internal labels.  If this number is
   *  set to zero (the default choice) then return.
   */

  if(reset == 0) {
    if(setup->nintlab == UNSET)
      setup->nintlab = 0;
    printf("\nNumber of internal labels: [%d] ",setup->nintlab);
    fgets(line,MAXC,stdin);
    if(line[0] != '\n') {
      while(sscanf(line,"%d",&setup->nintlab) != 1 || setup->nintlab < 0) {
        fprintf(stderr,"ERROR: Bad input format.  Enter number again:  ");
        fgets(line,MAXC,stdin);
      }
    }
    if(setup->nintlab == 0)
      return 0;

    /*
     * Allocate memory for requested internal label positions
     */

    if(!(setup->intlabs = new_labinfo(setup->nintlab))) {
      fprintf(stderr,"ERROR: get_intlab\n");
      return 1;
    }
  }
  
  /*
   * Fill the arrays
   */

  for(i=0,labptr=setup->intlabs; i<setup->nintlab; i++,labptr++) {
    printf("  Enter x y and text for internal label %d, separated ",i+1);
    if(reset)
      printf("by spaces: [%f %f %s] ",labptr->x,labptr->y,labptr->text);
    else
      printf("by spaces:  ");
    fgets(line,MAXC,stdin);
    if(reset == 0 || (reset == 1 && line[0] != '\n')) {
      while(sscanf(line,"%f %f %s",&labptr->x,&labptr->y,
  		 labptr->text) != 3) {
        fprintf(stderr,"ERROR: incorrect input format.\n");
        fprintf(stderr,"Enter positions and text as 3 values, separated by ");
        fprintf(stderr,"spaces, as \"25 134 A\":  ");
        fgets(line,MAXC,stdin);
      }
    }
  }

  return 0;
}

/*.......................................................................
 *
 * Function get_labscale
 *
 * Gets internal label scaling in the form of character height and
 *  character line width.
 *
 * Input:  Setup *setup         image display information
 *
 * Output: int (0 or 1)         0 ==> success, 1 ==> error
 *
 * v12Oct99 CDF, Cleaned up code and corrected problems with default choices.
 */

int get_labscale(Setup *setup)
{
  char line[MAXC];    /* General string for reading input */

  /*
   * Get character height, first setting the value to default
   */

  if(setup->labheight < 0.0) {
    setup->labheight = 2.0;
    printf("\nEnter character height for internal labels.\n");
    printf(" Input is a floating point number:  [%4.1f] ",setup->labheight);
    fgets(line,MAXC,stdin);
    if(line[0] != '\n') {
      while(sscanf(line,"%f",&setup->labheight) != 1 || setup->labheight < 0) {
	fprintf(stderr,"ERROR: Bad input format.   Enter value again:  ");
	fgets(line,MAXC,stdin);
      }
    }
  }
  else {
    printf("\nUsing character height = %6.2f from setup file\n",
	   setup->labheight);
  }

  /*
   * Get character width, first setting the value to default
   */

  if(setup->labwidth < 0.0) {
    setup->labwidth = 10;
    printf("\nEnter line width for internal label characters.\n");
    printf(" Input is an integer in range 1 - 210:  [%d] ",
	   setup->labwidth);
    fgets(line,MAXC,stdin);
    if(line[0] != '\n') {
      while(sscanf(line,"%d",&setup->labwidth) != 1 || 
	    setup->labwidth < 1 || setup->labwidth > 210) {
	fprintf(stderr,"ERROR: Bad input format.   Enter value again:  ");
	fgets(line,MAXC,stdin);
      }
    }
  }
  else {
    printf("\nUsing character width = %d from setup file\n",
	   setup->labwidth);
  }

  return 0;
}

/*.......................................................................
 *
 * Function get_border_info
 *
 * Gets the parameters governing the plotting of the borders.
 *
 * Inputs: Image *image        image to be displayed
 *         Setup *setup        setup info
 *         
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 * v20Apr2004 CDF, Changed labeling such that image gets the "Relative R.A."
 *                  and "Relative Dec" label if it is within 1 degree of
 *                  north-up/east-left.  Otherwise, label the axes as
 *                  delta_x and delta_y
 */

int get_border_info(Image *image, Setup *setup)
{
  int labchoice;        /* Flag used to set labeling style */
  char line[MAXC];      /* General string for reading input */
  char scalestr[MAXC];  /* String used in axis labeling for scales */

  if(setup->border == UNSET) {
    printf("\nInclude border and external labels? [y] ");
    fgets(line,MAXC,stdin);
    if(line[0] == 'n' || line[0] == 'N') {
      setup->border = FALSE;
      return 0;
    }
    else
      setup->border = TRUE;
  }

  if(setup->axislab == PIXEL) {
    sprintf(setup->xlab,"Pixels");
    sprintf(setup->ylab,"Pixels");
  }
  else {
    if(setup->axislab == ARCMIN)
      sprintf(scalestr,"(arcmin)");
    else
      sprintf(scalestr,"(arcsec)");
    if(setup->instrument==DSS || setup->instrument==RADIO || 
       image->rotation*180.0/PI<1.0) {
      sprintf(setup->xlab,"Relative R.A. %s",scalestr);
      sprintf(setup->ylab,"Relative Dec. %s",scalestr);
    }
    else {
      sprintf(setup->xlab,"\\gD x %s",scalestr);
      sprintf(setup->ylab,"\\gD y %s",scalestr);
    }
  }

  return 0;
}

/*.......................................................................
 *
 * Function get_wfpc_overlay
 *
 * Interactively gets the setup parameters that determine how the WFPC2
 *  overlay will be plotted.
 *
 * Inputs: Setup *setup        setup info
 *         
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int get_wfpc_overlay(Setup *setup)
{
  char line[MAXC];  /* General string for getting input */

  /*
   * Get central position in pixels
   */

  printf("\n--------------------------------------------------\n\n");
  printf("Set WFPC overlay parameters:\n");
  printf("  To keep current values, just hit return.\n");
  printf("  Central position in pixels (x y): ");
  printf("[%5.0f %5.0f] ",setup->wfpccent.x,setup->wfpccent.y);
  fgets(line,MAXC,stdin);
  if(line[0] != '\n') {
    while(sscanf(line,"%lf %lf",&setup->wfpccent.x,&setup->wfpccent.y) != 2) {
      fprintf(stderr,"\nERROR.  Bad input format. Enter position again:  ");
      fgets(line,MAXC,stdin);
    }
  }

  /*
   * Get position angle
   */

  printf("  Position angle (PA_V3) degrees: [%4.1f] ",setup->wfpcpa);
  fgets(line,MAXC,stdin);
  if(line[0] != '\n') {
    while(sscanf(line,"%f",&setup->wfpcpa) != 1) {
      fprintf(stderr,"ERROR.  Bad input format. Enter PA again:  ");
      fgets(line,MAXC,stdin);
    }
  }

  return 0;
}
