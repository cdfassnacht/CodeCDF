/* testplt.c
 *
 * Usage: testplt filename [setup_filename]
 *
 * Description:  Takes an image in FITS format and produces a greyscale
 *                plot.  Used to test various functions before they are
 *                incorporated into other programs. 
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
  int catformat;        /* Catalog format */
  int autopars=ARCSEC;  /* Flag set to PIXEL if WCS read fails */
  int screenflag=1;     /* Flag set to 0 for no screen output */
  int fileflag=1;       /* Flag set to 0 for no screen output */
  int usefile=0;        /* Flag set to 1 if there is a setup file */
  int another_loop=1;   /* Flag set to 0 to stop looping */
  char file_name[MAXC]; /* The name of the FITS file to load */
  char setupname[MAXC]; /* Name of the optional setup file */
  char plotfile[MAXC];  /* Filename for postscript plot */
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

  if(argc < 4 || argc > 5) {
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
  strcpy(catfile,argv[2]);
  if(sscanf(argv[3],"%d",&catformat) != 1)
    no_error = 0;
  else
    printf("Catalog format = %d\n",catformat);

  if(argc == 5 && no_error) {
    strcpy(setupname,argv[4]);
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

  if(plot_fitsfile_init(file_name,setupname,image,setup,usefile,1,screenflag)
     == ERROR) {
    fprintf(stderr,"\nERROR: Exiting testplt.\n");
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
    fprintf(stderr,"\n\nERROR: testplt.  Could not read fits file %s\n\n",
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
    open_plot_window();
#if 0
    if(setup->instrument != DSS) {
      printf("\nDisplay image to the screen? [y] ");
      fgets(line,MAXC,stdin);
      if(line[0] == 'n' || line[0] == 'N')
	screenflag = 0;
      else
	open_plot_window();
    }
#endif
  }

  /*
   * Display the image
   */

  if(screenflag && no_error) {
    printf("\nDisplaying image...\n");
    if(display_image(image,setup))
      no_error = 0;
  }

  /*
   * Read in SExtractor catalog file
   */

  printf("\n");
  if(no_error)
    if(!(catdat = read_secat(catfile,'#',&ncat,catformat)))
      no_error = 0;

  /*
   * Test mouse click return
   */

  while(another_loop && no_error) {
    curpos = return_position_from_click("",1,&status);
    no_error = 1 - status;
    if(no_error) {
      bestmatch = find_closest(curpos,catdat,ncat,1);
    }
    if(curpos.label[0]=='X' || curpos.label[0]=='x')
      another_loop = 0;
  }


  /*
   * Clean up and exit
   */

  if(no_error) {
    printf("\nHit return to exit: ");
    fgets(line,MAXC,stdin);
  }
  cpgend();
  image = del_Image(image);
  setup = del_setup(setup);
  catdat = del_secat(catdat);

  if(no_error) {
    printf("\nProgram testplt.c completed.\n\n");
    return 0;
  }
  else {
    fprintf(stderr,"\nERROR. Exiting testplt.c\n\n");
    return 1;
  }
}


int MAIN_(void)
{
  return 0;
}
