/* catplot.c
 *
 * Usage: catplot fitsfile xycatfile ([setup_filename])
 *
 * Description:  Plots a FITS file as a greyscale plot and then marks
 *                the positions of the objects detected by SExtractor.
 *
 * To compile use the appropriate makefile in this directory.
 *
 * Revision history
 *  29May2004 CDF,  First working version
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
 * Main program
 *
 */

int main(int argc, char *argv[])
{
  int no_error=1;       /* Flag set to 0 on error */
  int usefile=0;        /* Flag set to 1 to use input setup file */
  int screenflag=1;     /* Flag set to 0 for no screen output */
  int fileflag=1;       /* Flag set to 0 for no screen output */
  int ncat;             /* Number of catalog objects */
  int another_loop=1;   /* Flag set to 0 to stop loop */
  char fitsfile[MAXC];  /* The name of the FITS file to load */
  char catfile[MAXC];   /* The name of the SExtractor catalog file */
  char setupfile[MAXC]; /* Name of the optional setup file */
  char plotfile[MAXC];  /* Filename for postscript plot */
  char plotdev[MAXC];   /* PGPLOT display device */
  Image *image;         /* The container of the FITS image array */
  Setup *setup;         /* Container for the image display info */
  Secat *catdat=NULL;   /* SExtractor catalog data */

  /*
   * Check the command line
   */

  if(argc < 3) {
    fprintf(stderr,"\ncatplot fitsfile xycatfile (setupfile)\n\n");
    return 1;
  }

  /*
   * Get the names of the files.
   */

  strcpy(fitsfile,argv[1]);
  strcpy(catfile,argv[2]);
  if(argc == 4) {
    usefile = 1;
    strcpy(setupfile,argv[3]);
  }
  else
    sprintf(setupfile,"");

  /*
   * Load the image and print out basic image information
   */

  if(!(image = new_Image(fitsfile))) {
    fprintf(stderr,"\n\nERROR: fitsplt.  Could not read fits file %s\n\n",
	    fitsfile);
    return 1;
  }

  if(load_image_data(image) == ERROR)
    no_error = 0;

  /*
   * Create setup structure and fill with information used to set plotting
   *  parameters.
   */

  if(no_error)
    if(!(setup = fill_setup(image,setupfile,usefile,1)))
      no_error = 0;

  /*
   * Loop over displaying image until happy with results
   */

  if(no_error) {
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

  /*
   * Read in SExtractor catalog file
   */

  if(no_error)
    if(!(catdat = read_secat(catfile,'#',&ncat,2)))
      no_error = 0;

  /*
   * Plot the catalog
   */

  if(no_error) {
    printf("\nPlotting SExtractor catalog positions\n");
    if(plot_secat(setup,catdat,ncat,2,3,2.0))
      no_error = 0;
  }

  cpgend();

  /*
   * Now see if an output postscript or GIF file is desired.
   * NB: postscript image automatically printed for DSS images.
   */

  if(no_error) {
    fileflag = get_plotname(image,setup,"postscript",plotfile);
    if(fileflag) {
      printf("\nWriting to file %s\n",plotfile);
      sprintf(plotdev,"%s/vps",plotfile);
      if(cpgbeg(0,plotdev,1,1) != 1)
	no_error = 0;
      if(no_error)
	if(display_image(image,setup))
	  no_error = 0;
      if(no_error) 
	if(plot_secat(setup,catdat,ncat,2,1,2.0))
	  no_error = 0;
      cpgend();
    }
  }

  if(no_error) {
    fileflag = get_plotname(image,setup,"GIF",plotfile);
    if(fileflag) {
      printf("\nWriting to file %s\n",plotfile);
      sprintf(plotdev,"%s/gif",plotfile);
      if(cpgbeg(0,plotdev,1,1) != 1)
	no_error = 0;
      if(no_error)
	if(display_image(image,setup))
	  no_error = 0;
      if(no_error)
	if(plot_secat(setup,catdat,ncat,2,3,2.0))
	  no_error = 0;
    }
  }

  /*
   * Clean up and exit
   */

  cpgend();
  image = del_Image(image);
  setup = del_setup(setup);
  catdat = del_secat(catdat);

  if(no_error) {
    printf("\nProgram catplot.c completed.\n\n");
    return 0;
  }
  else {
    fprintf(stderr,"\nERROR. Exiting catplot.c\n\n");
    return 1;
  }
}

int MAIN_(void)
{
  return 0;
}
