/* plot_phot.c
 *
 * Usage: plot_phot fitsfile catfile [setup_filename]
 *  NB: catalog format is hardwired to 17, since that is what is needed for
 *      this program.
 *
 * Description: Plots a fits image of a standard star field.  The user clicks
 *               to identify the star(s) of interest, and data on the star(s)
 *                is shown in separate plots.  In particular, the curve
 *                of growth, in magnitudes, is shown so that the user can
 *                verify the total instrumental magnitude of the standard star.
 *
 * To compile use the appropriate makefile in this directory.
 *
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "libfits.h"
#include "cpgplot.h"
#include "structdef.h"
#include "dataio.h"
#include "catlib.h"
#include "fitsfuncs.h"
#include "fitswcs.h"
#include "fitsim.h"
#include "plotfuncs.h"

/*.......................................................................
 *
 * Function declaration
 *
 */

int plot_curve_of_growth(Secat object, int nmags);

/*.......................................................................
 *
 * Main program
 *
 */

int main(int argc, char *argv[])
{
  int no_error=1;       /* Flag set to 0 on error */
  int status;           /* Return value from function */
  int ncat=0;           /* Number of objects in catalog */
  int catformat=17;     /* Catalog format */
  int autopars=ARCSEC;  /* Flag set to PIXEL if WCS read fails */
  int fileflag=1;       /* Flag set to 0 for no screen output */
  int usefile=0;        /* Flag set to 1 if there is a setup file */
  int another_loop=1;   /* Flag set to 0 to stop looping */
  char fitsfile[MAXC];  /* The name of the FITS file to load */
  char setupname[MAXC]; /* Name of the optional setup file */
  char plotdev[MAXC];   /* PGPLOT display device */
  char catfile[MAXC];   /* The name of the SExtractor catalog file */
  char line[MAXC];      /* General string for getting input */
  Pos curpos;           /* Cursor position */
  Fits *fits=NULL;      /* Fits container */
  Image *image=NULL;    /* The container of the FITS image array */
  Setup *setup=NULL;    /* Container for the image display info */
  Secat *catdat=NULL;   /* SExtractor catalog data */
  Secat bestmatch;      /* Best match to clicked position */

  /*
   * Check the command line
   */

  if(argc < 3 || argc > 4) {
    setup_help();
    return 1;
  }
  else {
    fitsplt_copyright_notice();
  }

  /*
   * Parse the command line
   */

  strcpy(fitsfile,argv[1]);
  strcpy(catfile,argv[2]);
#if 0
  if(sscanf(argv[3],"%d",&catformat) != 1)
    no_error = 0;
  else if(catformat != 17) {
    fprintf(stderr,"ERROR. Catalog must be in format 17 for this program.\n");
    no_error = 0;
  }
  else
    printf("Catalog format = %d\n",catformat);
#endif

  if(argc == 4 && no_error) {
    strcpy(setupname,argv[3]);
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
  printf("Input fits file: %s\n",fitsfile);
#if 0
  printf("\nDisplay image to the screen? [y] ");
  fgets(line,MAXC,stdin);
  if(line[0] == 'n' || line[0] == 'N')
    screenflag = 0;
#endif
  /*
   * Load data and display if desired
   */
#if 0

  if(plot_fitsfile_init(fitsfile,setupname,image,setup,usefile,1,screenflag)
     == ERROR) {
    fprintf(stderr,"\nERROR: Exiting plot_phot.\n");
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

  if(!(image = new_Image(fitsfile))) {
    fprintf(stderr,"\n\nERROR: plot_phot.  Could not read fits file %s\n\n",
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
    if(!(setup = fill_setup(image,setupname,usefile,1)))
      no_error = 0;

  /*
   * Open the screen display device if one is desired
   */

  if(no_error)
    open_plot_window();

  /*
   * Display the image
   */

  /*
   * Loop over displaying image until happy with results
   */

  if(no_error) {
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
  another_loop = 1;

  /*
   * Read in SExtractor catalog file
   */

  printf("\n");
  if(no_error)
    if(!(catdat = read_secat(catfile,'#',&ncat,catformat)))
      no_error = 0;

  /*
   * Plot the catalog
   */

  if(no_error) {
    printf("Plotting SExtractor catalog positions\n");
    if(plot_secat(setup,catdat,ncat,2,3,2.0))
      no_error = 0;
  }

  /*
   * Use mouse to select objects, and then plot curve of growth and radial
   *  plot centered on those objects.
   */

  if(no_error) {
    open_plot_window();
    open_plot_window();
    cpgslct(1);
  }


  while(another_loop && no_error) {

    /* Get cursor position */

    printf("Left-click to select object. Right-click to quit.\n");
    curpos = return_position_from_click("",1,&status);
    no_error = 1 - status;

    /* Find closest match and plot curve of growth */

    if(no_error && curpos.label[0] !='X' && curpos.label[0] !='x') {
      bestmatch = find_closest(curpos,catdat,ncat,1);
      cpgslct(2);
      cpgeras();
      cpgscf(2);
      plot_curve_of_growth(bestmatch,20);
      cpgslct(3);
      cpgeras();
      cpgscf(2);
      curpos.x = bestmatch.x;
      curpos.y = bestmatch.y;
      radplt(image,curpos,7);
      cpgslct(1);
    }
    if(curpos.label[0]=='X' || curpos.label[0]=='x')
      another_loop = 0;
  }


  /*
   * Clean up and exit
   */
#if 0
  if(no_error) {
    printf("\nHit return to exit: ");
    fgets(line,MAXC,stdin);
  }
#endif
  cpgend();
  image = del_Image(image);
  setup = del_setup(setup);
  catdat = del_secat(catdat);

  if(no_error) {
    printf("\nProgram plot_phot.c completed.\n\n");
    return 0;
  }
  else {
    fprintf(stderr,"\nERROR. Exiting plot_phot.c\n\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function: plot_curve_of_growth
 *
 * Inputs:
 *  Secat object               catalog object for which COG is to be plotted
 *  int nmags                  number of points in the curve of growth
 *
 * Output:
 *  int (SUCCESS or ERROR)
 */

int plot_curve_of_growth(Secat object, int nmags)
{
  int i;                /* Looping variable */
  float x[nmags];       /* x-vector for plotting */
  char title[MAXC];     /* Title string */

  /*
   * Create the x-axis vector
   */

  for(i=0; i<nmags; i++)
    x[i] = 2+i;

  /*
   * Plot the curve of growth
   */

  sprintf(title,"Curve of growth for object %d",object.id);
  plot_xyerr(x, object.maper, object.mapererr, nmags, "Radius (pix)",
	     "Instrumental magnitude",title,0,1);

  return SUCCESS;
}

int MAIN_(void)
{
  return 0;
}
