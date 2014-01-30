#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "structdef.h"
#include "dataio.h"
#include "coords.h"
#include "libfits.h"
#include "cpgplot.h"
#include "plotfuncs.h"
#include "fitsfuncs.h"
#include "fitswcs.h"
#include "fitsim.h"

/*
 * fitsim.c
 *
 * This is a collection of functions that deal with the reading and
 *  display of images that are stored as FITS files.
 *
 * Function listing:
 * -----------------
 *  fitsplt_copyright_notice
 *  new_Image
 *  del_Image
 *  imstat
 *  im_rms
 *  im_median
 *  sigclip
 *  display_image
 *  adjust_display
 *  calc_ibound
 *  doborder
 *  draw_compass
 *  compass_angle_keck
 *  compass_angle_wfpc2
 *  compass_angle_nicmos
 *  draw_locator
 *  plot_circles
 *  plot_secat
 *  read_slitmask
 *  plot_slitmask
 *  draw_lines_setup
 *  draw_lines_interactive
 *  set_line_attrib
 *  draw_line
 *  plot_wfpc
 *  label_ra_dec
 *  print_image_data
 *  get_plotname
 *  imslice
 *  find_vslice_max
 *  make_vslice
 *  plot_vslice
 *  fit_continuum
 *  pos_lsf
 *  floatcmp
 *
 *-----------------------------------------------------------------------
 * Revision history:
 * -----------------
 *  20Sep96 Chris Fassnacht (CDF) - First working version
 *  v03May97 CDF, Created a structure for the setup information and added
 *                 functions to label the plots and to draw a compass
 *                 to put on NIRC images.
 *  v24Jul97 CDF, Moved the section for getting and checking image size
 *                 from setup_interact to its own function, get_imsize.
 *  v25Jul97 CDF, Imported section from findplt.c that marks an object in
 *                 the field (previously used only for finding charts).
 *  v26Jul97 CDF, Incorporated all of the findplt.c functions in a general
 *                 way so that findplt.c can be deleted and fitsplt.c can
 *                 be run to produce DSS finding charts.
 *                Incorporated all of the nircplt.c functions in a general
 *                 way so that nircplt.c can be deleted and fitsplt.c can
 *                 be run to produce plots for NIRC.
 *                Added capability to produce plots for LRIS and CCD13 with
 *                 default pixel scales, etc. for those instruments.
 *                Added choice of putting in a compass for non-DSS plots.
 *  v01Aug97 CDF, Now takes source name from object header in FITS file, or
 *                 gets the source name from the input file name in the case
 *                 of DSS files.
 *                Added choice of setting contours as multiples of some
 *                 base level.
 *                Better error checking.
 *  v02Aug97 CDF, First implementation of setup file input.  Still very crude.
 *  v12Aug97 CDF, Added option of marking target object.  Program assumes
 *                 that object is in the middle of the field.
 *  v13Aug97 CDF, Change marked object position to be the requested central
 *                 pixel.
 *                Changed marker size in draw_locator from fixed pixel values
 *                 to fraction of image size.
 *  v18Aug97 CDF, Added function to mark positions of objects selected by
 *                 the "autoslit" program.
 *                Set up framework for more general setup file input
 *  v13Oct97 CDF, Fixed bug in get_intlab function.
 *                Created a new Labinfo structure which contains the
 *                 positions and text for internal labels.  Structure is
 *                 defined in structdef.c
 *  v14Oct97 CDF, Added capacity to set image transfer function (to linear,
 *                 log or sqrt), to turn off image border and external
 *                 labels and to change the size of the internal labels.
 *                Added a much more detailed help file.
 *  v22Oct97 CDF, Added WFPC2 compass function.
 *  v12Jan98 CDF, Added capability to draw lines on plot which can be used to
 *                 mark, for example, the boundary of the WFPC2 chips on a
 *                 plot.
 *  v13Jan98 CDF, Fixed bug in draw_lines_interactive.
 *  v19Feb98 CDF, Added marker position option to setup file parameters.
 *                Added AUTO option for logical parameter values.  This will
 *                 be used to differentiate between, the DSS values and other
 *                 values.
 *  v24Feb98 CDF, Changed setup structure from a xstart, xend, ystart, yend
 *                 format to a centx, centy, xsize, ysize format
 *  v25Feb98 CDF, Fixed bug in read_setup_line.
 *                Added capability for pixel scale and compass flag to be read 
 *                 from setup file.
 *  v01Apr98 CDF, Fixed bug with image centering for full-image request.
 *  v07May98 CDF, Added necessary dynamic memory allocation to get_contour.
 *  v12May98 CDF, Linked to plotfuncs library.
 *                Made line widths for border thicker in doborder.
 *                Changed default for setup->title to null string.
 *                Modified plot_slitmask to accept input with no name
 *                 parameter.
 *  v01Jun98 CDF, Fixed bug in get_imsize in which the default central y 
 *                 position was being set to half of the X size rather than 
 *                 half the Y size.
 *  v03Aug98 CDF, Added print_image_data as a diagnostic function that allows
 *                 the printing of the values in image->data to an ASCII
 *                 output file.
 *  v19Aug98 CDF, Made lines wider in draw_locator, plot_slitmask, 
 *                 draw_lines_setup, and draw_lines_interactive.
 *  v20Aug98 CDF, Added WFPC2-PC pixel scale.
 *  v30Sep98 CDF, Fixed bugs with DSS processing.
 *  v09Dec98 CDF, Added read_file_pixscale because SkyView returns a
 *                 constant-sized (in PIXELS) image with varying pixel
 *                 scales set by the desired ANGULAR size of the image. (!)
 *  v25Mar99 CDF, Allowed for FITS images with more than 2 dimensions by 
 *                 truncating these images after the first two dimensions.
 *                Modified read_file_pixscale to deal with non-DSS images.
 *  v08Apr99 CDF, Moved fitsplt copyright notice from fitsplt.c into a 
 *                 new function fitsplt_copyright_notice.
 *                Fixed error in pixscale calculation for DSS files.
 *                Modified read_file_pixscale to deal with situations in which
 *                 a file has CDn_m header cards instead of CDELTn cards (in
 *                 a very crude way).
 *  v02May99 CDF, Added get_epoch function to get the value of the EPOCH
 *                 header card, if it exists.  This value can be used to 
 *                 distinguish between DSS and POSS2 images.
 *  v26May99 CDF, Fixed a bug with the marker position determination in
 *                 get_imsize
 *  v09Sep99 CDF, Added keywords and a function (plot_circles) to plot circles
 *                 of various sizes and with optional labels on the plot.
 *  v12Oct99 CDF, Moved pixel-scale determination into the get_pixscale function
 *                 and added a check for CCD13 in the algorithm.  This 
 *                 additional check is needed because the CCD13 pixel scale
 *                 is set to 1 pixel in the header instead of any useful
 *                 units.
 *                Added a parameter to display_image which sets the PGPLOT
 *                 display device.
 *  v16Oct99 CDF, Moved flux-range determination into get_frange function, 
 *                 which allows the possibility of adjusting the flux range
 *                 multiple times from a new function (adjust_display) in 
 *                 fitsplt.c
 *                Changed draw_lines_interactive to put the interactive
 *                 drawing info into the setup container at the end.  This
 *                 way the user doesn't have to enter the information multiple
 *                 time for the new multiple-display capability.
 *                Put cross-hair marker queries into get_marker_info function
 *                 to fix yet another bug in the marker positioning.
 *  v06Nov99 CDF, Added "quit" option to get_instrument.
 *  v12Nov99 CDF, Added "nolabel" and "nocontour" as possible setup-file inputs.
 *                Moved contour color index into get_contour function and added
 *                 the color index to the Setup structure definition.
 *  v29Nov99 CDF, Because the automatic determination of the pixel scale wasn't
 *                 working well for some instruments (e.g. WFPC2), made the
 *                 get_pixscale function print out both the instrument default
 *                 AND the scale read from the header cards.  The default option
 *                 is now the instrument default rather than the header card
 *                 value.
 *  v04Dec99 CDF, Split off contouring functions into fitsim_contours.c
 *                Fixed bug in reading "nolabel" keyword.
 *  v05Dec99 CDF, Fixed bug in reading internal labels.
 *                Split setup functions into fitsim_setup.c
 *                Added an Imaxis array to Image structure definition for use
 *                 in future features.  This requires some modification of
 *                 read_file_pixscale.
 *  v19Dec99 CDF, Better handling of CDi_j header cards inf read_file_pixscale.
 *                More information to the user in new_Image, including values
 *                 of some standard keywords.  Convert CRVALi pixels to sky
 *                 positions if they are of type RA---??? or DEC--??.
 *                Changed default value in get_instrument to be an actual guess
 *                 of the instrument based on the fits header cards.
 *  v09Feb2000 CDF, Split off all the functions that interactively get 
 *                   information for the Setup container into get_params.c
 *  v15Feb2001 CDF, Added calculation of image mean and rms as well as 
 *                   sigma-clipping of the data through the new functions
 *                   imstat, im_rms, and sigclip.
 *  v11Oct2001 CDF, Added option to change line width (via the new setup
 *                   lweight parameter).
 *  v22Apr2002 CDF, Added new print_image_info function.
 *  v04Mar2003 CDF, Added new get_plotname function.
 *  v25Aug2003 CDF, Added fitsra, fitsdec members to Image structure, which
 *                   allows modification of image_info and print_image_info
 *                   functions.
 *                  Added axis_info function, which modifies image_info.
 *  v20Apr2004 CDF, Changes to draw_compass
 *  v28May2004 CDF, Modifications to display_image, get_imsize, and
 *                   plot_circles.
 *                  Moved adjust_display from fitsplt.c to this library.
 *                  Added the plot_secat function.
 *  v09Jul2007 CDF, Small change in adjust_display
 *  v08Jul2008 CDF, Major changes to plot_slitmask to take as input a RA,Dec
 *                   catalog rather than a x,y catalog.
 *  v10Jul2008 CDF, Moved new_Imaxis and del_Imaxis into fitswcs.c in
 *                   CDFfits directory
 *                  Moved new_Image and del_Image into fitsfuncs.c in
 *                   CDFfits directory
 *                  Combined image_info and print_image_info into the new
 *                   print_fits_info function in fitsfunc.c in CDFfits
 *                  Got rid of axis_info and moved its functionality into
 *                   read_wcs_info in fitswcs.c
 */

/*.......................................................................
 *
 * Function fitsplt_copyright_notice
 *
 * Prints out copyright notice and helpful info.  Moved over from
 *  fitsplt.c on 08Apr99.
 *
 * Inputs: none
 *
 * Output: none
 *
 */

void fitsplt_copyright_notice()
{
  printf("\n*******************************************************");
  printf("****************\n\n");
  printf("  Program: fitsplt     Author: Chris Fassnacht\n");
  printf("  Version: 2008Jul10\n");
  printf("  Copyright (c) 2002-2008 Regents of the University of California\n");
  printf("  Uses fitslib library developed by Martin Shepherd\n");
  printf("  Copyright (c) 1993 California Institute of Technology\n\n");
  printf("  Usage: fitsplt fits_filename (setup_filename)\n\n");
  printf("** For help with setup file, just type \"fitsplt\" and hit ");
  printf("return **\n\n");
  printf("*******************************************************");
  printf("****************\n\n");
}

/*.......................................................................
 *
 * Function plot_fitsfile_initial
 *
 * Reads a fits file, fills the setup container, displays the image, and
 *  lets the user adjust the image.  In other words, calls a lot of 
 *  functions that were called individually before
 *
 * Inputs: char *fitsfile       name of input fits file
 *         char *setupname      name of setup file (can be "")
 *         Image *image         image container (filled by this function)
 *         Setup *setup         setup container (filled by this function)
 *         int usefile          set to 1 to use setup file
 *         int autopars         set to 1 to set many parameters automatically
 *         int screenflag       set to 1 to display image to screen
 *
 * Output: int (SUCCESS or ERROR)
 */

int plot_fitsfile_init(char *fitsfile, char *setupfile, Image *image,
		       Setup *setup, int usefile, int autopars,
		       int screenflag)
{
  int no_error=1;       /* Flag set to 0 on error */
  int another_loop=1;   /* Flag set to 0 to stop looping */

  /*
   * Read the file into the Image container
   */

  if(!(image = new_Image(fitsfile))) {
    fprintf(stderr,"ERROR.  Exiting program. Could not open %s\n\n",
	    fitsfile);
    return ERROR;
  }

  /*
   * Print out basic file information
   */

  if(load_image_data(image) == ERROR)
    no_error = 0;

  /*
   * Create setup structure and fill with information used to set plotting
   *  parameters.
   */

  if(no_error) 
    if(!(setup = fill_setup(image,setupfile,usefile,autopars)))
      no_error = 0;
  printf("Setup instrument = %d\n",setup->instrument);

  /*
   * Loop over displaying image until happy with results
   */

  if(no_error && screenflag) {
    open_plot_window();
    if(display_image(image,setup))
      no_error = 0;
    else
      while(no_error && another_loop) {
	switch(another_loop = setup_menu(image,setup)) {
	case 1:
	  if(display_image(image,setup))
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

  if(no_error)
    return SUCCESS;
  else
    return ERROR;

}

/*.......................................................................
 *
 * Function display_image
 *
 * Displays an image previously read from a FITS file.
 *
 * Inputs: Image *image         image to be displayed
 *         Setup *setup         image display information
 *
 * Output: int (0 or 1)         0 ==> success, 1 ==> error
 *
 * v28May2004 CDF, Moved setting of rxstart, etc. from this function to
 *                  get_imsize in get_params.c.
 *                 Took opening of PGPLOT device out of this function again
 */

int display_image(Image *image, Setup *setup)
{
  int i;
  int ixstart,ixend,iystart,iyend;    /* Start and end of array (int) */
  float tr[6]={0,1,0,0,0,1};          /* Transformation matrix for pggray */
  float epoch=0.0;                    /* Epoch of observation */
  Labinfo *labptr;                    /* Pointer for setup->intlabs */
  
  /*
   * Set up the PGPLOT device
   */

  cpgask(0);
  cpgpap(0.0,1.0);
  cpgpage();
  cpgsch(setup->cheight);
  cpgscf(2);
  cpgvstd();

  /*
   * Set up the image coordinate system as pixels.
   */

  ixstart = calc_ibound(setup->rxstart,image->nx);
  ixend = calc_ibound(setup->rxend,image->nx);
  iystart = calc_ibound(setup->rystart,image->ny);
  iyend = calc_ibound(setup->ryend,image->ny);

  /*
   * Draw the image as a gray-scale array.
   */

  cpgwnad(setup->rxstart,setup->rxend,setup->rystart,setup->ryend);
  if(setup->trans >= 0 && setup->trans < 3)
    cpgsitf(setup->trans);
  else
    cpgsitf(0);
  cpggray(image->data,image->nx,image->ny,ixstart,ixend,
	  iystart,iyend,
	  setup->hi,setup->lo,tr);

  /*
   * Put contours on the image, if wanted.
   */

  if(setup->ncont > 0) {
    cpgsci(setup->contcolor);
    cpgcont(image->data,image->nx,image->ny,ixstart,ixend,
	    iystart,iyend,
	    setup->cont,setup->ncont,tr);
  }

  /*
   * If requested, draw a compass rose on the image
   */

  if(setup->drawcompass == TRUE)
    if(draw_compass(image,setup)) {
      fprintf(stderr,"ERROR: display_image\n");
      return 1;
    }

  /*
   * If requested, put in a marker to locate the object of interest
   */

  if(setup->drawmarker == TRUE) {
    printf("\nDrawing marker at location %7.1f %7.1f\n",setup->markx,
	   setup->marky);
    draw_locator(setup);
  }

  /*
   * Put in internal labels
   */

  cpgsci(1);
  if(setup->nintlab) {
    if(setup->labheight < 0)
      setup->labheight = 2.0;
    if(setup->labwidth < 1 || setup->labwidth > 201)
      setup->labwidth = 10;
    cpgsch(setup->labheight);
    cpgslw(setup->labwidth);
    for(i=0,labptr=setup->intlabs; i<setup->nintlab; i++,labptr++) {
      cpgptxt(labptr->x,labptr->y,0.,0.,labptr->text);
    }
  }

  /*
   * Plot circles on image
   */

  if(setup->ncircle) {
    if(plot_circles(setup)) {
      fprintf(stderr,"ERROR: display_image\n");
      return 1;
    }
  }

  /*
   * Mark objects in slitmask
   */

  if(strcmp(setup->slitfile,"#") != 0) {
    if(plot_slitmask(image,setup)) {
      fprintf(stderr,"ERROR: display_image\n");
      return 1;
    }
  }

  /*
   * Label the boundaries
   * Note that this has to be done after drawing the gray-scale since the
   * gray-scale is opaque and otherwise would erase the tick marks.
   */

  cpgslw(1);
  cpgsch(1.);
  if(setup->border > 0) {
    doborder(setup);
    cpgwnad(setup->rxstart,setup->rxend,setup->rystart,setup->ryend);
  }

  /*
   * Plot WFPC2 outline on image, if requested
   */

  if(setup->plotwfpc)
    plot_wfpc(setup);

  /*
   * Draw any additional lines on the plot if requested
   */

  if(setup->instrument != DSS) {
    switch(setup->ndraw) {
    case UNSET:
      draw_lines_interactive(setup);
      break;
    case 0:
      break;
    default:
      draw_lines_setup(setup);
    }
  }

  /*
   * Put RA and Dec in the bottom right corner for DSS images
   *  (not yet implemented).
   */

  if(setup->instrument == DSS) {
    label_ra_dec(image);
  }

  return 0;
}

/*.......................................................................
 *
 * Function adjust_display
 *
 * Adjusts the display limits through calls to get_frange.  This allows
 *  the user to fiddle with the display limits and see the effects on
 *  the screen before plotting the image to a file.  The new values for
 *  the display range are stored in the Setup container so that they
 *  can be used in the final plot.
 *
 * Inputs: Image *image        container for image data and header info
 *         Setup *setup        image display information (modified by this
 *                              function).
 * Output: 0 or 1              0 ==> success, 1 ==> error
 *
 * v09Jul2007 CDF, Changed call to get_frange to reflect change in get_params.c
 */

int adjust_display(Image *image, Setup *setup)
{
  int no_error=1;      /* Flag set to 0 on error */
  int do_adjust=1;     /* Flag set to 0 when no more adjustment is needed */
  char line[MAXC];     /* String variable for reading input */

  /*
   * First see if the user wants to change the current values, unless
   *  the image is a DSS or POSS2 image.
   */

  if(setup->instrument == DSS)
    do_adjust = 0;

  while(do_adjust && no_error) {
    printf("\n--------------------------------------------------\n\n");
    printf("\nNow you have the option of changing the display parameters\n");
    printf(" such as the display limits (image stretch)");
    if(setup->plotwfpc)
      printf(" and the WFPC2 overlay parameters.\n");
    else
      printf(".\n");
    printf("\nDo you want to change the image parameters? [y] ");
    fgets(line,MAXC,stdin);
    if(line[0] == 'n' || line[0] == 'N')
      do_adjust = 0;
    else {
      if(get_frange(image,setup))
	no_error = 0;
      if(setup->plotwfpc)
	if(get_wfpc_overlay(setup))
	  no_error = 0;
      if(no_error)
	if(display_image(image,setup))
	  no_error = 0;
    }
  }

  if(no_error)
    return 0;
  else {
    fprintf(stderr,"ERROR: adjust_display\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function radplt
 *
 * Given an (x,y) position, plots pixel values as a function of distance
 *  from the position.
 *
 * Inputs: Image *image        image container
 *         Pos centpos         position at center
 *         float rmax          maximum radius
 *
 * Output: int (SUCCESS or ERROR)
 */

int radplt(Image *image, Pos centpos, float rmax)
{
  int i,j,k;          /* Looping variables */
  int x,y;            /* Current pixel position */
  float tmpdat;       /* Current pixel value */
  int arrsize;        /* Array size */
  int xmin,xmax;      /* Range of good x values */
  int ymin,ymax;      /* Range of good y values */
  int tmp;            /* Temporary container */
  float dx,dy;        /* Distance of current pixel from centpos */
  float *r=NULL;      /* Container for distances from centpos */
  float *flux=NULL;   /* Container for data in the valid pixels */
  float *dptr,*fptr;  /* Pointers to navigate arrays */

  /*
   * Allocate memory for quantities to be plotted
   */

  arrsize = 2*rmax+1;
  if(!(r = new_array(arrsize,arrsize))) {
    fprintf(stderr,"ERROR: radplt\n");
    return ERROR;
  }
  if(!(flux = new_array(arrsize,arrsize))) {
    fprintf(stderr,"ERROR: radplt\n");
    return ERROR;
  }
  arrsize *= arrsize;

  /*
   * Set limits
   */

  printf("%f %f\n",centpos.x,centpos.y);
  xmin = ((tmp = (int) centpos.x - rmax) > 0) ? tmp : 0;
  xmax = ((tmp = (int) centpos.x + rmax) < image->nx) ? tmp : image->nx;
  ymin = ((tmp = (int) centpos.y - rmax) > 0) ? tmp : 0;
  ymax = ((tmp = (int) centpos.y + rmax) < image->ny) ? tmp : image->ny;
  printf("radplt: xmin, xmax, ymin, ymax = %d %d %d %d\n",xmin,xmax,ymin,ymax);
  printf("radplt: %d %d %d\n",image->nx,image->ny,image->ntotal);

  /*
   * Initialize pointers
   */

  dptr = r;
  fptr = flux;

  /*
   * Fill arrays by finding valid pixels in the data
   */
  tmp = 1;
  for(x=xmin; x<=xmax; x++) {
    for(y=ymin; y<=ymax; y++) {
      i = image->nx * (y-1) + x - 1;
      dx = x - centpos.x;
      dy = y - centpos.y;
      *dptr = sqrt(dx*dx + dy*dy);
      dptr++;
      *fptr = image->data[i];
      fptr++;
      tmp++;
    }
  }

  /*
   * Plot flux vs. r
   */

  plot_xy(r,flux,arrsize,"Radius (pix)","Counts","Radial Plot",0,0);

  /*
   * Clean up and exit
   */

  r = del_array(r);
  flux = del_array(flux);

  return SUCCESS;
}


/*.......................................................................
 *
 * Function calc_ibound
 *
 * Calculates one of the integer boundary value to be passed to cpggrey
 *  given the real world value requested.
 *
 * Input:  float rval          real world value for boundary
 *         int imax            maximum possible integer value
 *
 * Output: int ival            integer boundary value
 *
 */

int calc_ibound(float rval, int imax)
{
  int ival;   /* Integer boundary value */
  int tmp;    /* Temoporary value for test */

  tmp = (int) rval;

  if(tmp > imax)
    ival = imax;
  else if(tmp < 1)
    ival = 1;
  else
    ival = tmp;

  return ival;
}

/*.......................................................................
 *
 * Function doborder
 *
 * Labels the boundary of the area being displayed
 *
 * Inputs: Setup *setup        image display info
 *
 * Output: None
 *
 */

void doborder(Setup *setup)
{
  float arrxstart,arrxend,arrystart,arryend;

  /*
   * Set up variables used in the labelling
   */

  if(setup->axislab == PIXEL) {
    arrxstart = 1.0;
    arrxend = setup->xsize;
    arrystart = 1.0;
    arryend = setup->ysize;
  }
  else {
    if(setup->axcentset == TRUE) {
      arrxstart = (setup->rxstart - setup->axcentx) * setup->pixscale;
      arrxend = (setup->rxend - setup->axcentx) * setup->pixscale;
      arrystart = (setup->rystart - setup->axcenty) * setup->pixscale;
      arryend = (setup->ryend - setup->axcenty) * setup->pixscale;
    }
    else {
      arrxstart=(setup->xsize/2.0)*setup->pixscale;
      arrystart=(-setup->ysize/2.0)*setup->pixscale;
      arrxend=(-setup->xsize/2.0)*setup->pixscale;
      arryend=(setup->ysize/2.0)*setup->pixscale;
    }
  }

  /*
   * Label the boundaries
   */

  cpgslw(2);
  cpgsch(setup->cheight);
  cpgwnad(arrxstart,arrxend,arrystart,arryend);
  cpgbox("BCNST", 0, 0, "BCNST", 0, 0);
  cpglab(setup->xlab,setup->ylab,setup->title);
  cpgsch(1.0);
  cpgslw(1);
}

/*.......................................................................
 *
 * Function draw_compass
 *
 * Puts arrows indicating North and East on the plot
 *
 * Inputs: Image *image        image being displayed
 *         Setup *setup        display information
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 * v20Apr2004 CDF, Make the default behavior to use the rotation angle
 *                  determined from the CD matrix, if available.
 *
 */

int draw_compass(Image *image, Setup *setup)
{
  int flip=0;         /* Flag set to 1 if EW axis is flipped (as in WFPC2) */
  double alpha=0.0;   /* Angle of rotation (CCW) of N from top of array */
  float compx,compy;  /* Position for compass */
  float dc;           /* Shift from center of compass box to compass origin */
  float dnx,dny;      /* Change in x and y positions to end of N axis */
  float dex,dey;      /* Change in x and y positions to end of E axis */
  char line[MAXC];    /* General string for reading input */

  /*
   * If the file has a CD matrix, use the image rotation calculated
   *  from the matrix.
   */

  if(image->wcs.validwcs) {
    alpha = image->wcs.rotation;
    printf("\ndraw_compass: Using alpha = %7.2f\n",alpha*180.0/PI);
  }

  /*
   * If there is no CD matrix, get rotation from the instrument-specific
   *  keywords.  It is assumed that DSS and CCD13 images have no rotation.  
   *  NIRC and LRIS images need to have the proper rotation calculated 
   *  from the ROTPOSN header card in the FITS file.
   */

  else {
    switch(setup->instrument) {
    case CCD13: case DSS: case RADIO:
      alpha = 0.0;
      flip = 0;
      break;
    case NIRC: case LRIS:
      if(compass_angle_keck(image,setup->instrument,&alpha)) {
	fprintf(stderr,"ERROR: draw_compass\n");
	return 1;
      }
      flip = 0;
      break;
    case WF: case PC:
      if(compass_angle_wfpc2(image,&alpha)) {
	fprintf(stderr,"ERROR: draw_compass\n");
	return 1;
      }
      flip = 1;
      break;
    case NIC1:
      if(compass_angle_nicmos(image,&alpha)) {
	fprintf(stderr,"ERROR: draw_compass\n");
	return 1;
      }
      flip = 0;
      break;
    case OTHER:
      flip = 0;
      break;
    default:
      flip = 0;
    }

    /*
     * Set the rotation parameter to the alpha derived from the image header
     */

    image->rotation = alpha;
  }

  /*
   * Calculate position for compass
   */

  compx = setup->centx - 7.0*setup->xsize/18.0;
  compy = setup->centy + 7.0*setup->ysize/18.0;
  dc = setup->xsize * 0.8 / 12.0;
  if(sin(alpha) >= 0.0 && cos(alpha) >= 0.0) {
    if(flip) {
      compx += dc + 4*dc*(alpha - PI/2.0)/PI;
      compy -= dc;
    }
    else {
      compx += dc;
      compy += -dc + 4*dc*alpha/PI;
    }
  }
  else if(sin(alpha) >= 0.0 && cos(alpha) < 0.0) { 
    if(flip) {
      compx += dc;
      compy += -dc + 4*dc*(alpha - PI/2.0)/PI;
    }
    else {
      compx += dc - 4*dc*(alpha - PI/2.0)/PI;
      compy += dc;
    }
  }
  else if(sin(alpha) < 0.0 && cos(alpha) < 0.0) {
    if(flip) {
      compx += dc - 4*dc*(alpha - PI)/PI;
      compy += dc;
    }
    else {
      compx -= dc;
      compy += -3*dc - 4*dc*alpha/PI;
    }
  }
  else {
    if(flip) {
      compx -= dc;
      compy += -3*dc - 4*dc*(alpha-PI/2.0)/PI;
    }
    else {
      compx += dc + 4*dc*alpha/PI;
      compy -= dc;
    }
  }

  /*
   * Calculate end positions of axes
   */

  dc *= 1.5;
  dnx = -dc * sin(alpha);
  dny = dc * cos(alpha);
  if(flip) {
    dex = dny;
    dey = -dnx;
  }
  else {
    dex = -dny;
    dey = dnx;
  }

  /*
   * Draw the compass
   */

  cpgslw(5);
  cpgmove(compx,compy);
  cpgdraw(compx+dnx,compy+dny);
  cpgmove(compx,compy);
  cpgdraw(compx+dex,compy+dey);

  /*
   * Label the compass
   */

  dnx *= 1.35;
  dny *= 1.35;
  dex *= 1.35;
  dey *= 1.35;
  cpgsch(1.0);
  cpgptxt(compx+dnx,compy+dny,0.,0.5,"N");
  cpgptxt(compx+dex,compy+dey,0.,0.5,"E");
  cpgslw(1);

  return 0;
}

/*.......................................................................
 *
 * Function compass_angle_keck
 *
 * Calculates the compass angle based on the FITS header card ROTPOSN
 *  and the instrument (NIRC or LRIS), since there are different
 *  formulae used to convert ROTPOSN to P.A. for NIRC and LRIS.
 *
 * Inputs: Image *image        image being displayed
 *         int instrument      instrument used
 *         double *alpha       P.A. (set by this function)
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int compass_angle_keck(Image *image, int instrument, double *alpha)
{
  float rotposn;      /* Value of the ROTPOSN header in FITS file */
  char line[MAXC];    /* General string used for reading input */
  Fitkey key;         /* Used for reading FITS header */

  /*
   * Get the ROTPOSN header value from the FITS file.  If the keyword
   *  isn't there (get_key returns a value > 0), then prompt for it.
   */

  if(get_key(image->fits,image->fits->hdu,"ROTPOSN",DAT_DBL,LOOP_SEEK,
	     &key)) {
    printf("\ncompass_angle_keck: Warning -- ROTPOSN head was not found!\n");
    printf("Enter the value of the ROTPOSN header:  [0.0]");
    rotposn = 0.0;
    fgets(line,MAXC,stdin);
    if(line[0] != '\n') {
      while(sscanf(line,"%f",&rotposn) != 1) {
	fprintf(stderr,"ERROR: Bad input.  Enter ROTPOSN again:  \n");
	fgets(line,MAXC,stdin);
      }
    }
  }
  else
    rotposn = KEYDBL(key);

  /*
   * Convert ROTPOSN to an angle (north through east) in radians
   *  using the formula appropriate to the instrument.
   * NB: The math library function fmod returns the remainder
   *      remainder in a division.  Thus the first calculation of
   *      alpha will give a value between -2*PI and 2*PI.
   */

  switch(instrument) {
  case NIRC:
    *alpha = 2.0 * PI * (fmod((180.0 - rotposn),360.0))/360.0;
    break;
  case LRIS:
    *alpha = 2.0 * PI * (fmod((270.0 - rotposn),360.0))/360.0;
    break;
  default:
    fprintf(stderr,"\nERROR: compass_angle_keck.  Instrument not recognized.\n");
    return 1;
  }

  /*
   * Put alpha into the range -PI to PI
   */

  if(*alpha > PI)
    *alpha -= 2.0 * PI;
  else if(*alpha < -1.0*PI)
    *alpha += 2.0 * PI;

  printf("\ncompass_angle_keck: ");
  printf("ROTPOSN = %6.1f deg ==> alpha = %6.1f deg = %6.3f rad\n",rotposn,
	 *alpha * 180.0/PI,*alpha);

  return 0;
}

/*.......................................................................
 *
 * Function compass_angle_wfpc2
 *
 * Calculates the compass angle based on the FITS header card PA_V3
 *  found in the WFPC2 files.  The header is defined as follows
 *  (according to WFPC2 documentation):
 *
 *    1. The chips are arranged as follows:
 *
 *          2   1    
 *                     where chip 1 is the PC
 *          3   4
 *
 *    1a.  However, it looks like the standard pipeline processing may
 *          flip the image about the y-axis so that you end up with:
 *
 *          1   2    
 *
 *          4   3
 *
 *    2. Draw an axis diagonally through chip 3 to chip 1.  This
 *        axis points south if PA_V3 = 0.
 *
 *    3. PA_V3 increases clockwise from the "north" axis through chip 2.
 *
 *    4. If PA_V3 = 0, the perpendicular axis pointing through chip 4
 *        points east,  i.e., for PA_V3 = 0, the directions are:
 *
 *          S   W
 *
 *          E   N
 *
 * Inputs: Image *image        image being displayed
 *         double *alpha       P.A. (set by this function)
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int compass_angle_wfpc2(Image *image, double *alpha)
{
  float pav3;         /* Value of the PA_V3 header in FITS file */
  char line[MAXC];    /* General string used for reading input */
  Fitkey key;         /* Used for reading FITS header */

  /*
   * Get the PA_V3 header value from the FITS file.  If the keyword
   *  isn't there (get_key returns a value > 0), then prompt for it.
   */

  if(get_key(image->fits,image->fits->hdu,"PA_V3",DAT_DBL,LOOP_SEEK,
	     &key)) {
    printf("Warning -- PA_V3 header was not found!\n");
    printf("Enter the value of the PA_V3 header:  [0.0]");
    pav3 = 0.0;
    fgets(line,MAXC,stdin);
    if(line[0] != '\n') {
      while(sscanf(line,"%f",&pav3) != 1) {
	fprintf(stderr,"ERROR: Bad input.  Enter ROTPOSN again:  \n");
	fgets(line,MAXC,stdin);
      }
    }
  }
  else
    pav3 = KEYDBL(key);

  /*
   * Convert PA_V3 to an angle (north through east) in radians
   *  using the formula appropriate to the instrument.  Remember that
   *  PA_V3 = 0 gives north pointing to 225 degrees.
   * NB: The math library function fmod returns the remainder
   *      remainder in a division.  Thus the first calculation of
   *      alpha will give a value between -2*PI and 2*PI.
   */

  *alpha = 2.0 * PI * (fmod((225.0 - pav3),360.0))/360.0;

  /*
   * Put alpha into the range -PI to PI
   */

  if(*alpha > PI)
    *alpha -= 2.0 * PI;
  else if(*alpha < -1.0*PI)
    *alpha += 2.0 * PI;

  printf("PA_V3 = %6.1f deg ==> alpha = %6.1f deg = %6.3f rad\n",pav3,
	 *alpha * 180.0/PI,*alpha);

  return 0;
}

/*.......................................................................
 *
 * Function compass_angle_nicmos
 *
 * Calculates the compass angle based on the FITS header card ORIENTAT.
 * NB: For the NICMOS cameras, the ORIENTAT keyword gives the position
 *     angle of the image y axis, E of N.  So, that means that the 
 *     PA of the north axis should be -ORIENTAT.
 *
 * Inputs: Image *image        image being displayed
 *         double *alpha       P.A. (set by this function)
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int compass_angle_nicmos(Image *image, double *alpha)
{
  float orientat;     /* Value of the ORIENTAT header in FITS file */
  char line[MAXC];    /* General string used for reading input */
  Fitkey key;         /* Used for reading FITS header */

  /*
   * Get the ROTPOSN header value from the FITS file.  If the keyword
   *  isn't there (get_key returns a value > 0), then prompt for it.
   */

  if(get_key(image->fits,image->fits->hdu,"ORIENTAT",DAT_DBL,LOOP_SEEK,
	     &key)) {
    printf("Warning -- ORIENTAT head was not found!\n");
    printf("Enter the value of the ORIENTAT header:  [0.0]");
    orientat = 0.0;
    fgets(line,MAXC,stdin);
    if(line[0] != '\n') {
      while(sscanf(line,"%f",&orientat) != 1) {
	fprintf(stderr,"ERROR: Bad input.  Enter ORIENTAT again:  \n");
	fgets(line,MAXC,stdin);
      }
    }
  }
  else
    orientat = KEYDBL(key);

  /*
   * Convert ORIENTAT to an angle (north through east) in radians.
   * NB: The math library function fmod returns the remainder
   *      remainder in a division.  Thus the first calculation of
   *      alpha will give a value between -2*PI and 2*PI.
   */

  *alpha = 2.0 * PI * (fmod(-orientat,360.0))/360.0;

  /*
   * Put alpha into the range -PI to PI
   */

  if(*alpha > PI)
    *alpha -= 2.0 * PI;
  else if(*alpha < -1.0*PI)
    *alpha += 2.0 * PI;

  printf("ORIENTAT = %6.1f deg ==> alpha = %6.1f deg = %6.3f rad\n",orientat,
	 *alpha * 180.0/PI,*alpha);

  return 0;
}

/*.......................................................................
 *
 * Function draw_locator
 *
 * Draws a locator on a field, to mark an object of interest.
 *
 * Inputs: Setup *setup        plotting information
 *
 * Outputs: none
 *
 */

void draw_locator(Setup *setup)
{
  int imsize; /* Maximum image dimesion */
  float in;   /* Inner radius of marker cross */
  float out;  /* Outer radius of marker cross */

  /*
   * Increase the line width for visibility
   */

  cpgslw(5);

  /*
   * Set imsize to max image dimension
   */

  if(setup->xsize > setup->ysize)
    imsize = setup->xsize;
  else
    imsize = setup->ysize;

  /*
   * Set up inner and outer positions of cross as fraction of imsize
   */

  in = 0.01*imsize;
  out = 0.03*imsize;

  /*
   * Draw the cross
   */

  cpgmove(setup->markx-in,setup->marky-in);
  cpgdraw(setup->markx-out,setup->marky-out);
  cpgmove(setup->markx-in,setup->marky+in);
  cpgdraw(setup->markx-out,setup->marky+out);
  cpgmove(setup->markx+in,setup->marky+in);
  cpgdraw(setup->markx+out,setup->marky+out);
  cpgmove(setup->markx+in,setup->marky-in);
  cpgdraw(setup->markx+out,setup->marky-out);

  /*
   * Return to normal line width
   */

  cpgslw(1);
}

/*.......................................................................
 *
 * Function plot_circles
 *
 * Marks the positions of objects selected by autoslit and guide stars.
 * The input file should be the "nsetup_1.dat" file that is output from
 *  the autoslit program.
 *
 * Inputs: Setup *setup        plotting info
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 * v28May2004 CDF, removed the passed rxstart, etc. paramters since
 *                  they are now contained within the setup structure
 */

int plot_circles(Setup *setup)
{
  int i;            /* Looping variable */
  int no_error=1;   /* Flag set to 0 on error */
  Labinfo *circptr; /* Pointer to navigate setup->circle */

  /*
   * Set line width, fill style, character font and character height
   */

  cpgslw(2);
  cpgsfs(2);
  cpgsch(0.9);
  cpgsci(3);

  /*
   * Plot circles
   */

  for(i=0,circptr=setup->circle; i<setup->ncircle; i++,circptr++) {
    if(circptr->x > setup->rxstart && circptr->x < setup->rxend &&
       circptr->y > setup->rystart && circptr->y < setup->ryend) {
      cpgcirc(circptr->x,circptr->y,circptr->size);
      cpgtext(circptr->x+0.71*circptr->size,circptr->y+0.71*circptr->size,
	      circptr->text);
    }
  }

  /*
   * Set fill style, line width, and character height back to default values
   */

  cpgsfs(1);
  cpgslw(1);
  cpgsch(1.0);
  cpgsci(1);

  /*
   * Clean up and exit
   */


  if(no_error)
    return 0;
  else {
    fprintf(stderr,"ERROR: plot_circles\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function plot_secat
 *
 * Marks the positions of objects in a SExtractor-produced catalog.
 *  The format paramter describes the format of the input catalog.
 * Mark colors:
 *  0 white (foreground)
 *  1 black (background)
 *  2 red
 *  3 green
 *  4 blue
 *  5 cyan
 *  6 magenta
 *  7 yellow
 *  8 orange
 *  9-15, see http://www.astro.caltech.edu/~tjp/pgplot/fig51.html
 *
 * Inputs: Setup *setup        plotting info
 *         Secat *secat        catalog
 *         int ncat            number of catalog members
 *         int format          format of catalog
 *         int color           color of marks
 *         float scale         amount to scale size of marks
 *
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int plot_secat(Setup *setup, Secat *secat, int ncat, int format, int color,
	       float scale)
{
  int i;            /* Looping variable */
  int no_error=1;   /* Flag set to 0 on error */
  float circsize;   /* Size of circles */
  Secat *sptr;      /* Pointer to navigate secat */

  /*
   * Set line width, fill style, character font,color, and character height
   */

  cpgslw(1);
  cpgsfs(2);
  cpgsci(color);
  cpgsch(0.9);

  /*
   * Plot circles
   */

  for(i=0,sptr=secat; i<ncat; i++,sptr++) {
    if(sptr->x > setup->rxstart && sptr->x < setup->rxend &&
       sptr->y > setup->rystart && sptr->y < setup->ryend) {
      if(sptr->fwhm == 0.0)
	circsize = 0.015 * setup->ysize;
      else
	circsize = scale * sptr->fwhm;
      cpgcirc(sptr->x,sptr->y,circsize);
#if 0
      cpgtext(sptr->x+0.71*sptr->size,sptr->y+0.71*sptr->size,
	      sptr->text);
#endif
    }
  }

  /*
   * Set fill style, line width, color and character height back to 
   *  default values
   */

  cpgsfs(1);
  cpgslw(1);
  cpgsci(1);
  cpgsch(1.0);

  /*
   * Clean up and exit
   */


  if(no_error)
    return 0;
  else {
    fprintf(stderr,"ERROR: plot_secat\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function plot_slitmask
 *
 * Marks the positions of objects selected by autoslit and guide stars.
 *
 * Inputs: 
 *  Image *image               fits file
 *  Setup *setup               plotting info
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int plot_slitmask(Image *image, Setup *setup)
{
  int i;             /* Looping varible */
  int no_error=1;    /* Flag set to 0 on error */
  float r;           /* Radius of circles or half-length of box sides */
  float labfctr;     /* Dictates how far labels are placed from objects */
  Secat *sptr;       /* Pointer to navigate setup->slitmask */

  /*
   * TEMPORARY: Read in slitmask info
   */

  if(read_slitmask(image,setup) == ERROR)
    return ERROR;

  /*
   * Set circle radius to be a fraction of image size
   */

  r = 0.015 * setup->ysize;
  labfctr = 1.0;

  /*
   * Set line width, fill style, character font and character height
   */

  cpgslw(3);
  cpgsfs(2);
  cpgsch(0.9);

  /*
   * Mark slits with circles and stars (priority = -1) with boxes
   */

  if(no_error) {
    for(i=0,sptr=setup->slitmask; i<setup->nmask; i++,sptr++) {

      /*
       * Plot circles if priority > 0, where priority is approximated by
       *  sptr->merr.
       */
	
      if(sptr->merr > 0) {
	cpgcirc(sptr->x,sptr->y,r);
	cpgtext(sptr->x + labfctr*r,sptr->y,sptr->name);
      }

      /*
       * Else plot boxes around alignment stars
       */

      else {
	cpgrect(sptr->x - r,sptr->x + r,sptr->y - r,sptr->y + r);
	cpgtext(sptr->x + labfctr*r,sptr->y,sptr->name);
      }
    }
  }

  /*
   * Set fill style, line width, and character height back to default values
   */

  cpgsfs(1);
  cpgslw(1);
  cpgsch(1.0);

  /*
   * Clean up and exit
   */

  if(no_error)
    return 0;
  else {
    fprintf(stderr,"ERROR: plot_slitmask\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function draw_lines_setup
 *
 * Draws line segments indicated in setup file
 *
 * Inputs: Setup *setup        plotting info
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int draw_lines_setup(Setup *setup)
{
  int i;          /* Looping variable */
  Pos *startptr;  /* Pointer to starting positions */
  Pos *endptr;    /* Pointer to ending positions */

  /*
   * Set up line attributes
   */

  set_line_attrib(setup);

  /*
   * Loop through line segments
   */

  for(i=0,startptr=setup->drawstart,endptr=setup->drawend; i<setup->ndraw;
      i++,startptr++,endptr++) {
    draw_line(*startptr,*endptr);
  }

  /*
   * Return line attributes to default values
   */

  cpgslw(1);

  return 0;
}

/*.......................................................................
 *
 * Function draw_lines_interactive
 *
 * Draws line segments on plot
 *
 * Inputs: Setup *setup        plotting info
 *
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 * v16Oct99 CDF, Modified to put the line segments chosen here into
 *                the Setup container so that they can be plotted
 *                on the final plot, if desired.
 */

int draw_lines_interactive(Setup *setup)
{
  int i=0;           /* Counting variable */
  int no_error=1;    /* Flag set to 0 on error */
  char line[MAXC];   /* General string variable for reading input */
  Pos start[MAXC];   /* Starting positions */
  Pos end[MAXC];     /* Ending positions */
  Pos *sptr,*eptr;   /* Pointers for navigating start and end arrays */
  Pos *ssptr,*septr; /* Pointers for navigating drawstart and drawend */
  

  printf("\nDraw lines on the plot, e.g. to mark chip boundaries? [n] ");
  fgets(line,MAXC,stdin);
  if(line[0] != 'y' && line[0] != 'Y') {
    setup->ndraw = 0;
    return 0;
  }

  /*
   * Set up line attributes
   */

  set_line_attrib(setup);

  printf("\nEntering line drawing mode.\n\n");

  /*
   * Loop through line segments
   */

  printf("Enter start and end positions of lines in the following format: ");
  printf("\nx_start y_start x_end y_end.\n\n");
  printf("At any time, type q to quit out of line drawing mode\n\n");

  sptr = start;
  eptr = end;

  while(i>=0 && i<MAXC) {
    printf("Enter coordinates for line segment %d:  ",i+1);
    fgets(line,MAXC,stdin);
    if(line[0] == 'q')
      break;
    else {
      while(sscanf(line,"%lf %lf %lf %lf",&sptr->x,&sptr->y,
		   &eptr->x,&eptr->y) != 4) {
	fprintf(stderr,"ERROR: Bad input format.  Enter values again:  ");
	fgets(line,MAXC,stdin);
      }
      draw_line(*sptr,*eptr);
      i++;
      sptr++;
      eptr++;
    }
  }

  /*
   * Allocate memory for the setup drawstart and drawend arrays
   */

  setup->ndraw = i;
  if(!(setup->drawstart = new_pos(setup->ndraw,1)))
    no_error = 0;
  if(!(setup->drawend = new_pos(setup->ndraw,1)))
    no_error = 0;

  /*
   * Put the line segment values into the setup container
   */

  if(no_error) {
    for(i=0,ssptr=setup->drawstart,septr=setup->drawend,sptr=start,eptr=end;
	i<setup->ndraw; i++,ssptr++,septr++,sptr++,eptr++) {
      *ssptr = *sptr;
      *septr = *eptr;
    }
  }

  /*
   * Return line attributes to default values
   */

  cpgslw(1);

  if(no_error) {
    printf("\ndraw_lines_interactive: Stored %d lines in setup container.\n",
	   setup->ndraw);
    return 0;
  }
  else {
    fprintf(stderr,"ERROR: draw_lines_interactive.\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function set_line_attrib
 *
 * Sets up line attributes.
 *
 * Inputs: Setup *setup        plotting info
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int set_line_attrib(Setup *setup)
{
  char line[MAXC];   /* General string for reading input */

  if(setup->lweight == UNSET) {
    setup->lweight = 3;
    printf("\n");
    printf("Enter line width as an integer between 1 and 201: [%d] ",
	   setup->lweight);
    fgets(line,MAXC,stdin);
    if(line[0] != '\n') {
      while(sscanf(line,"%d",&setup->lweight) != 1 || setup->lweight < 1
	    || setup->lweight > 201) {
	fprintf(stderr,"ERROR: set_line_attrib.  Bad input.  ");
	fprintf(stderr,"Enter new value:  ");
	fgets(line,MAXC,stdin);
      }
    }
  }
  cpgslw(setup->lweight);

  return 0;
}

/*.......................................................................
 *
 * Function draw_line
 *
 * Given start and end positions for the line segment, will draw a line
 *  with the current line width and style.
 *
 * Inputs: Pos begin           beginning position of line segment
 *         Pos end             end position of line segment
 *
 * Outputs: none
 *
 */

void draw_line(Pos begin, Pos end)
{
  cpgmove(begin.x,begin.y);
  cpgdraw(end.x,end.y);
}

/*.......................................................................
 *
 * Function plot_wfpc
 *
 * Plots the WFPC2 outline onto the image.  The lines are defined in
 *  angular units and the value in setup->pixscale is used to determine
 *  the number of pixels.  The values used to set the outline attributes
 *  are the central position and the position angle (PA_V3, as defined
 *  for WFPC2).
 *
 * Inputs: Setup *setup        setup info
 *         
 * Output: none
 *
 */

void plot_wfpc(Setup *setup)
{
  int no_error=1;   /* Flag set to 0 on error */
  float alpha=0.0;  /* PA in radians */
  float r45;        /* 45 degrees, in radians */
  float wfside;     /* Length of one side of WF chip in pixels */
  float pcside;     /* Length of one side of PC chip in pixels */
  float wfs2,pcs2;  /* Temporary working variables */
  Pos center;       /* Central position (dummy variable) */
  char line[MAXC];  /* General string for getting input */

  /*
   * Use the dummy variable center to save on typing
   */

  center = setup->wfpccent;

  /*
   * Convert PA_V3 to radians
   */

  alpha = setup->wfpcpa * PI / 180.0;
  r45 = 45.0 * PI / 180.0;

  /*
   * Calculate the length of the lines to draw, in pixels.  Do this
   *  by taking the chips to be 800 pixels on a side and multiplying
   *  by the appropriate pixel scale and then dividing by the pixel
   *  scale of the fits image.
   */

  wfside = 800.0 * PIXWF / setup->pixscale;
  pcside = 800.0 * PIXPC / setup->pixscale;
  wfs2 = wfside*sqrt(2.0);
  pcs2 = pcside*sqrt(2.0);

  /*
   * Draw the outlines.  Do PC1 first.
   */

  cpgmove(center.x,center.y);
  cpgdraw(center.x+pcside*sin(alpha),center.y+pcside*cos(alpha));
  cpgdraw(center.x+pcs2*cos(r45-alpha),
	  center.y+pcs2*sin(r45-alpha));
  cpgdraw(center.x+pcside*cos(alpha),center.y-pcside*sin(alpha));

  /*
   * WF2
   */

  cpgmove(center.x,center.y);
  cpgdraw(center.x+wfside*sin(alpha),center.y+wfside*cos(alpha));
  cpgdraw(center.x-wfs2*cos(r45+alpha),
	  center.y+wfs2*sin(r45+alpha));
  cpgdraw(center.x-wfside*cos(alpha),center.y+wfside*sin(alpha));


  /*
   * WF3
   */

  cpgmove(center.x,center.y);
  cpgdraw(center.x-wfside*cos(alpha),center.y+wfside*sin(alpha));
  cpgdraw(center.x-wfs2*cos(r45-alpha),
	  center.y-wfs2*sin(r45-alpha));
  cpgdraw(center.x-wfside*sin(alpha),center.y-wfside*cos(alpha));

  /*
   * WF4
   */

  cpgmove(center.x,center.y);
  cpgdraw(center.x-wfside*sin(alpha),center.y-wfside*cos(alpha));
  cpgdraw(center.x+wfs2*cos(r45+alpha),
	  center.y-wfs2*sin(r45+alpha));
  cpgdraw(center.x+wfside*cos(alpha),center.y-wfside*sin(alpha));
  cpgdraw(center.x,center.y);

}

/*.......................................................................
 *
 * Function label_ra_dec
 *
 * Puts a label in the bottom right-hand corner of DSS images giving the
 *  RA and Dec of the image.  The RA and Dec are calculated from the
 *  CRVALn header cards.
 *
 * Inputs: Image *image        image being displayed
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int label_ra_dec(Image *image)
{
  int i;               /* Looping variable */
  int no_error=1;      /* Flag set to 0 on error */
  double crval[2];     /* Coordinate values of the reference pixel */
  Phdu *phdu;          /* Primary data-header unit in the fits file */
  Imaxis *axis;        /* Image axis */
  Skypos skypos;       /* RA and Dec of the source */

  /*
   * Set phdu to point at the appropriate place
   */

  phdu = (Phdu *) image->fits->hdu;

  /*
   *  Get the appropriate values
   */

  for(i=0; i<2; i++) {
    if(!(axis=get_axis(phdu,i+1))) {
      printf("label_ra_dec:  ***WARNING. CRVAL%d not found.\n",i+1);
      no_error = 0;
    }
    else {
      if(i == 0) {
	printf("\n label_ra_dec: Reading image coordinates from ");
	printf("FITS header.\n");
      }
      printf(" label_ra_dec: CRVAL%d = %f degrees\n",
	     i+1,axis->crval);
      crval[i] = axis->crval;
    }
  }

  /*
   * Convert to RA and Dec
   */

  if(no_error) {
    printf(" label_ra_dec: Converting to hh:mm:ss format\n");
    deg2spos(crval[0],crval[1],&skypos);
    printf(" label_ra_dec: Position is %02d:%02d:%06.3f %+03d:%02d:%05.2f\n",
	   skypos.hr,skypos.min,skypos.sec,
	   skypos.deg,skypos.amin,skypos.asec);
  }

  return 0;
}

/*.......................................................................
 *
 * Function print_image_data
 *
 * Prints out the values in image->data.
 *
 * Inputs: Image *image        image being displayed
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int print_image_data(Image *image)
{
  int i;              /* Looping variable */
  float *dptr;        /* Pointer used to navigate image->data */
  char line[MAXC];    /* General string variable for reading input */
  char outname[MAXC]; /* Name of output file */
  FILE *ofp=NULL;     /* Output file pointer */

  /*
   * Initialize the output name to image.dat
   */

  sprintf(outname,"image.dat");

  /*
   * Get name of output file
   */

  printf("Enter name of output file for print_image_data: [%s]",outname);
  fgets(line,MAXC,stdin);
  if(line[0] != '\n')
    strcpy(outname,line);

  /*
   * Open output file
   */

  if(!(ofp = open_writefile(outname))) {
    fprintf(stderr,"ERROR: print_image_data\n");
    return 1;
  }

  /*
   * Loop through data file, printing out values
   */

  for(i=0,dptr=image->data; i < image->nx * image->ny; i++,dptr++)
    fprintf(ofp,"%f\n",*dptr);

  /*
   * Clean up and exit
   */

  if(ofp)
    fclose(ofp);

  return 0;
}

/*.......................................................................
 *
 * Function get_plotname
 *
 * Gets plot name for output file.  If the image is a DSS image, the
 *  plot name is set automatically from the input source name.
 *
 * Inputs: Image *image        input image
 *         Setup *setup        setup info
 *         char *plottype      what type of plot are we inquiring about
 *                              (postscript or gif)?
 *         char *plotname      output filename -- set by this function
 *
 * Output: int doplot          Set to 1 if a plot is desired, else set
 *                              to 0.
 *
 */

int get_plotname(Image *image, Setup *setup, char *plottype, char *plotname)
{
  int doplot=0;         /* Flag set to 1 if an output plot is desired */
  float epoch;          /* Epoch of observation -- needed for DSS plots */
  char line[MAXC];      /* General string for getting input */
  char ext[MAXC];       /* Default extension for output filename */

  /*
   * Set default extension
   */

  if(strcmp(plottype,"GIF") == 0)
    sprintf(ext,".gif");
  else
    sprintf(ext,".ps");

  /*
   * Handle automatic processing for DSS.  Only plot for postscript and
   *  not for GIF.
   */

  if(setup->instrument == DSS) {
    if(get_epoch(image,&epoch) == 0 && epoch > 1960.0)
      sprintf(plotname,"%s_fc_poss2_6.ps",setup->source);
    else {
      epoch = 0.0;
      sprintf(plotname,"%s_fc_dss_6.ps",setup->source);
    }
    if(strcmp(plottype,"postscript") == 0)
      doplot = 1;
  }

  /*
   * For all other instrument types, query the user to see if plot is
   *  desired
   */

  else {
    printf("\nPlot image to a %s file? [y] ",plottype);
    fgets(line,MAXC,stdin);
    if(line[0] == 'n' || line[0] == 'N') {
      doplot = 0;
    }
    else {
      doplot = 1;
      sprintf(plotname,"%s%s",setup->source,ext);
      printf("Enter output filename: [%s] ",plotname);
      fgets(line,MAXC,stdin);
      if(line[0] != '\n') {
	while(sscanf(line,"%s",plotname) != 1) {
	  fprintf(stderr,"ERROR: Bad input.  Enter filename again:  ");
	  fgets(line,MAXC,stdin);
	}
      }
    }
  }

  return doplot;
}

/*.......................................................................
 *
 * Function imslice
 *
 * Slices a two-d image into one-d cuts.
 *
 * Inputs: Image *image        input image
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int imslice(Image *image)
{
  int x,y;        /* Looping variables */
  int ymax;       /* Value of y associated with slice maximum */
  float slmax;    /* Maximum along slice */
  float sky=0.0;  /* Estimated sky level at ymax */
  float *dptr;    /* Pointer to navigate image->data */
  FILE *ofp=NULL;

  if(!(ofp=open_writefile("foo.dat"))) {
    fprintf(stderr,"ERROR: imslice\n");
    return 1;
  }

  /*
   * Cut vertical slices
   */

  for(x=0; x<image->nx; x++) {
    find_vslice_max(image,x,0,image->ny,&slmax,&ymax);
#if 0
    fit_continuum(image,x,ymax,-30.0,-15.0,15.0,30.0,&sky);
#endif
    slmax -= sky;
    fprintf(ofp,"%d %f %f\n",x,slmax,sky);
  }

  if(ofp)
    fclose(ofp);

  return 0;
}

/*.......................................................................
 *
 * Function find_vslice_max
 *
 * Finds the maximum value on a vertical slice through a 2D image.
 *
 * Inputs: Image *image        input image
 *         int x               column along which slice is taken
 *         int ylo, yhi        range of lines in slice
 *         float *slmax        maximum flux value (set by this function)
 *         int *ymax           y value associated with maximum flux
 *                              value (set by this function)
 *
 * Output: (none)
 */

void find_vslice_max(Image *image, int x, int ylo, int yhi, 
		     float *slmax, int *ymax)
{
  int y;        /* Looping variable */
  float *dptr;  /* Pointer to navigate image->data */

  *slmax = 0.0;
  for(y=ylo; y<yhi; y++) {
    dptr = image->data + x + y*image->nx;
    if(*dptr > *slmax) {
      *slmax = *dptr;
      *ymax = y;
    }
  }

  printf("\nfind_vslice_max: max at %d\n",*ymax);
  
}

/*.......................................................................
 *
 * Function make_vslice
 *
 * Creates a vertical slice through a 2D image given the column number
 *  (x0) and the range of y values (ymin to ymax).  The slice is returned
 *  as a Pos array containing (y,flux) pairs.
 *
 * Inputs: Image *image        input image
 *         int x0              column
 *         int ymin, ymax      range of y values
 *         float *smin,*smax   range of flux values (set by this function)
 *         double *ysmax       y position of smax (set by this function)
 *
 * Output: Pos *vslice         array containing slice.  NULL on error.
 *
 */

Pos *make_vslice(Image *image, int x0, int ymin, int ymax, float *smin,
		 float *smax, double *ysmax)
{
  int y;             /* Looping variable */
  int no_error=1;    /* Flag set to 0 on error */
  float ydiff;       /* ymax - ymin */
  float *dptr;       /* Pointer to navigate image->data */
  Pos *slice=NULL;   /* (y, flux) values */
  Pos *slptr;        /* Pointer to navigate slice */

  /*
   * Initialize
   */

  dptr = image->data + x0 + ymin*image->nx;
  *smin = *smax = *dptr;
  *ysmax = ymin;

  /*
   * Allocate memory for Pos array
   */

  ydiff = ymax - ymin;
  if(!(slice = new_pos(ydiff+1,1))) {
    fprintf(stderr,"ERROR: make_vslice.\n");
    return NULL;
  }
  
  /*
   * Step through slice, assigning (y,flux) pairs to the Pos array
   */

  for(y=ymin,slptr=slice; y<=ymax; y++,slptr++) {
    slptr->x = y;
    dptr = image->data + x0 + y*image->nx;
    slptr->y = *dptr;
    if(slptr->y > *smax) {
      *smax = slptr->y;
      *ysmax = slptr->x;
    }
    if(slptr->y < *smin)
      *smin = slptr->y;
  }

  if(no_error)
    return slice;
  else {
    fprintf(stderr,"ERROR: make_vslice.\n");
    return del_pos(slice);
  }
}


/*.......................................................................
 *
 * Function plot_vslice
 *
 * Plots a Pos array by calling plot_spec.
 *
 * Inputs: Pos *slice          slice to be plotted
 *         double ymin, ymax   range of y values
 *         float smin,smax     range of flux values
 *         int dopoints        flag set to 1 if plotting slice as
 *                              discrete points rather than connected
 *                              histogram style.
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int plot_vslice(Pos *slice, double ymin, double ymax, float smin,
		float smax, int dopoints)
{
  int i;             /* Looping variable */
  int intlo,inthi;   /* Integer values of ymin and ymax */
  float cheight=1.5; /* Character height scale for plot */
  float fmin,fmax;   /* Min and max values in array */
  float sdiff;       /* smax - smin */
  float ydiff;       /* ymax - ymin */
  Pos *sptr;         /* Pointer to navigate slice */

  /*
   * Get integer versions of ymin and ymax.
   */

  intlo = (int) ymin;
  inthi = (int) ymax;

  /*
   * Find min and max values in array
   */

  fmin = fmax = slice->y;

  for(i=intlo,sptr=slice+intlo; i<=inthi; i++,sptr++) {
    if(sptr->flag ==0 && sptr->y > fmax) {
      fmax = sptr->y;
    }
    if(sptr->flag ==0 && sptr->y < fmin)
      fmin = sptr->y;
  }

  /*
   * Initialize
   */

  ydiff = inthi - intlo;
  sdiff = fmax - fmin;

  /*
   * Call plot_spec to plot the spectrum
   */

  cpgscf(2);
  if(dopoints) {
    if(plot_specpts(slice+intlo,ydiff+1,fmin-0.1*sdiff,fmax+0.1*sdiff,
		    intlo,inthi,1,cheight)) {
      fprintf(stderr,"ERROR: plot_vslice.\n");
      return 1;
    }
    else
      return 0;
  }
  else {
    if(plot_spec(slice+intlo,ydiff+1,fmin-0.1*sdiff,fmax+0.1*sdiff,intlo,inthi,
		 1,1,cheight)) {
      fprintf(stderr,"ERROR: plot_vslice.\n");
      return 1;
    }
    else
      return 0;
  }
}

/*.......................................................................
 *
 * Function fit_continuum
 *
 * Computes a continuum level given two background regions 
 *  delineated by (bk_lo,bk_hi)_1 and (bk_lo,bk_hi)_2.
 *  * 
 * Inputs: Pos *slice          slice to have continuum fit
 *         int npoints         number of points in slice
 *         double bklo_1       lower limit of first background region
 *         double bkhi_1       upper limit of first background region
 *         double bklo_2       lower limit of second background region
 *         double bkhi_2       upper limit of second background region
 *         float *m, *b        slope and intercept of background (set
 *                              by this function)
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int fit_continuum(Pos *slice, int npoints, double bklo_1, double bkhi_1, 
		  double bklo_2, double bkhi_2, float *m, float *b)
{
  int y;           /* Looping variable */
  int count=0;     /* Number of values going into continuum calculation */
  Pos *bkgd=NULL;  /* Container for background values */
  Pos *bptr;       /* Pointer to navigate bkgd */
  Pos *sptr;       /* Pointer to navigate slice */

  /*
   * Allocate memory for bkgd.  For simplicity, make it the
   *  same size as slice.
   */

  if(!(bkgd = new_pos(npoints,1))) {
    fprintf(stderr,"ERROR: fit_continuum.\n");
    return 1;
  }

  /*
   * Load background info into container
   */

  bptr = bkgd;
  for(y=0,sptr=slice; y<npoints; y++,sptr++) {
    if((y >= bklo_1 && y<= bkhi_1) || 
       (y >= bklo_2 && y<= bkhi_2)) {
      *bptr = *sptr;
      count++;
      bptr++;
    }
  }

  /*
   * Do a linear least squares fit to background
   */

  if(pos_lsf(bkgd,count,m,b)) {
    fprintf(stderr,"ERROR: fit_continuum\n");
    return 1;
  }

  return 0;
}

/*.......................................................................
 *
 * Function pos_lsf
 *
 * Does a linear least-squares fit to the members of a Pos array.
 *
 * Inputs: Pos *pos            pairs of points in array
 *         int npoints         number of points in the array
 *         float *m            slope of linear fit (set by this 
 *                              function)
 *         float *b            intercept of linear fit (set by this
 *                              function)        set x2 = x * x
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int pos_lsf(Pos *pos, int npoints, float *m, float *b)
{
  int i;     /* Looping variable */
  float weight;
  float wtsum=0.0;
  float x2sum=0.0;
  float xsum=0.0;
  float ysum=0.0;
  float xysum=0.0;
  float del=0.0;
  Pos *pptr;       /* Pointer to navigate pos */

  /*
   * Assume that errors are shot-noise only, and thus that
   *  err^2 is just equal to the flux, i.e. pos->y.  Note that
   *  this isn't quite right since it doesn't take the CCD gain
   *  into account correctly.
   */

  for(i=0,pptr=pos; i<npoints; i++,pptr++) {
    printf("%f\n",pptr->y);
    weight = 1.0 / pptr->y;
    wtsum += weight;
    xsum += pptr->x * weight;
    ysum += pptr->y * weight;
    x2sum += pptr->x * pptr->x * weight;
    xysum += pptr->x * pptr->y * weight;
  }


  del = wtsum * x2sum - xsum * xsum;
  *b = (x2sum * ysum - xsum * xysum) / del;
  *m = (wtsum * xysum - xsum * ysum) / del;

  printf("\n%d %f %f %f %f %f %f\n",npoints,wtsum,xsum,ysum,x2sum,xysum,
	 del);

  return 0;
}
