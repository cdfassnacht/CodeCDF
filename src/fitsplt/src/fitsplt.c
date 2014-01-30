/* fitsplt.c
 *
 * Usage: fitsplt filename [setup_filename]
 *
 * Description:  Takes an image in FITS format and produces a greyscale
 *                plot along with contours and labels, if desired.  The
 *                program can be run either with a setup file (by specifying
 *                the optional setup file name) or interactively.
 *
 * To compile use the appropriate makefile in this directory.
 *
 * Revision History
 * -----------------
 * 1995Oct28 - MCS and CDF - First working version
 * 1996Sep20 - CDF - Moved new_Image and del_Image into fitsim.c so they can be
 *    used by other programs
 * 1997Jul26 - CDF - Incorporated all of the findplt.c functions in a general
 *    way so that findplt.c can be deleted and fitsplt.c can be run to produce
 *    DSS finding charts.
 *   Incorporated all of the nircplt.c functions in a general way so that
 *    nircplt.c can be deleted and fitsplt.c can be run to produce plots for
 *    NIRC.
 *   Added capability to produce plots for LRIS and CCD13 with
 *    default pixel scales, etc. for those instruments.
 *   Added choice of putting in a compass for non-DSS plots.
 * 1997Aug02 - CDF - First implementation of setup file input.  Still crude.
 * 1997Oct14 - CDF - Much better setup file input.
 *   Better help available when just typing "fitsplt"
 * 1997Oct22 - CDF - Added WFPC2 compass function.
 * 1999Oct13 - CDF - Added option of doing both a screen display and a 
 *    postscript file with the same setup parameters.
 * 1999Oct16 - CDF - Added function adjust_display to allow the user to change
 *    the display values interactively during the screen display.
 * 1999Nov06 - CDF - Added a clean exit if new "quit" option is chosen in
 *    the get_instrument function.
 * 2003Mar04 - CDF - Added an option to print out a GIF file as well as a
 *    postscript file.
 * 2004May28 - CDF - Moved adjust_display function into fitsim.c
 * 2004Jun08 - CDF - Modified displaying calls to reflect new form of the
 *    display_image function.
 * 2008Jul15 - CDF - Combined most of the individual setup filling functions
 *    into the fill_setup function.
 * 2008Jul16 - CDF - Replaced PGPLOT-opening code with a simple call to
 *    open_plot_window from plotfuncs.c
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "libfits.h"
#include "cpgplot.h"
#include "structdef.h"
#include "fitsfuncs.h"
#include "fitswcs.h"
#include "fitsim.h"
#include "plotfuncs.h"

/*.......................................................................
 *
 * Function declaration
 *
 */

/*.......................................................................
 *
 * Main program
 *
 */

int main(int argc, char *argv[])
{
  int no_error=1;       /* Flag set to 0 on error */
  int autopars=ARCSEC;  /* Flag set to PIXEL if WCS read fails */
  int screenflag=1;     /* Flag set to 0 for no screen output */
  int fileflag=1;       /* Flag set to 0 for no screen output */
  int usefile=0;        /* Flag set to 1 if there is a setup file */
  int another_loop=1;   /* Flag set to 0 to stop looping */
  char file_name[MAXC]; /* The name of the FITS file to load */
  char setupname[MAXC]; /* Name of the optional setup file */
  char plotfile[MAXC];  /* Filename for postscript plot */
  char plotdev[MAXC];   /* PGPLOT display device */
  char line[MAXC];      /* General string for getting input */
  Fits *fits=NULL;      /* Fits container */
  Image *image=NULL;     /* The container of the FITS image array */
  Setup *setup=NULL;    /* Container for the image display info */

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
   * Parse the command line
   */

  strcpy(file_name,argv[1]);

  if(argc == 3 && no_error) {
    strcpy(setupname,argv[2]);
    usefile = 1;
  }
  else {
    printf("No setup file.\n");
    sprintf(setupname,"");
  }

  /*
   * See if image should be displayed to the screen
   */

  printf("\n");
  printf("Input fits file: %s\n",file_name);
  printf("\nDisplay image to the screen? [y] ");
  fgets(line,MAXC,stdin);
  if(line[0] == 'n' || line[0] == 'N')
    screenflag = 0;
  /*
   * Load data and display if desired
   */
#if 0

  if(plot_fitsfile_init(file_name,setupname,image,setup,usefile,1,screenflag)
     == ERROR) {
    fprintf(stderr,"\nERROR: Exiting fitsplt.\n");
    return 1;
  }
  if(setup == NULL)
    printf("Setup is null\n");
  if(image == NULL)
    printf("Image is null\n");
#endif


  /*
   * Load the image and print out basic image information
   */

  if(!(image = new_Image(file_name))) {
    fprintf(stderr,"\n\nERROR: fitsplt.  Could not read fits file %s\n\n",
	    file_name);
    return 1;
  }

  if(load_image_data(image) == ERROR)
    no_error = 0;

  /*
   * Create setup structure and fill with information used to set plotting
   *  parameters.
   */

  if(no_error)
    if(!(setup = fill_setup(image,setupname,usefile,1)))
      no_error = 0;

  /*
   * Open the screen display device if one is desired
   */

  if(no_error) {
    if(setup->instrument != DSS) {
      printf("\nDisplay image to the screen? [y] ");
      fgets(line,MAXC,stdin);
      if(line[0] == 'n' || line[0] == 'N')
	screenflag = 0;
      else
	open_plot_window();
    }
  }

  /*
   * Loop over displaying image until happy with results
   */

  if(screenflag && no_error) {
    printf("\nDisplaying image...\n");
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



  /*
   * Now see if an output postscript or GIF file is desired.
   * NB: postscript image automatically printed for DSS images.
   */

  if(no_error) {
    fileflag = get_plotname(image,setup,"postscript",plotfile);
    if(fileflag) {
      printf("\nWriting to file %s\n",plotfile);
      sprintf(plotdev,"%s/vps",plotfile);
      printf("\nplotdev = %s\n",plotdev);
      if(cpgopen(plotdev) <= 0)
	while(cpgopen("?") <= 0)
	  fprintf(stderr,"ERROR: Not a supported device.  Try again");
      if(no_error)
	if(display_image(image,setup))
	  no_error = 0;
    }
    fileflag = get_plotname(image,setup,"GIF",plotfile);
    if(fileflag) {
      printf("\nWriting to file %s\n",plotfile);
      sprintf(plotdev,"%s/gif",plotfile);
      if(cpgopen(plotdev) <= 0)
	no_error = 0;
      if(no_error)
	if(display_image(image,setup))
	  no_error = 0;
    }
  }

  /*
   * Clean up and exit
   */

  cpgend();
  image = del_Image(image);
  setup = del_setup(setup);

  if(no_error) {
    printf("\nProgram fitsplt.c completed.\n\n");
    return 0;
  }
  else {
    fprintf(stderr,"\nERROR. Exiting fitsplt.c\n\n");
    return 1;
  }
}


int MAIN_(void)
{
  return 0;
}
