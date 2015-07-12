/* fitshead.c
 *
 * Usage: fitshead filename
 *
 * Description:  Prints out the header information for a fits file.
 *
 * To compile use the appropriate makefile in this directory.
 *
 * 10July2005 CDF, First working version
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "libfits.h"
#include "structdef.h"

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
  int i;                /* Looping variable */
  int no_error=1;       /* Flag set to 0 on error */
  int screenflag=1;     /* Flag set to 0 for no screen output */
  int fileflag=1;       /* Flag set to 0 for no screen output */
  char *headstr;     /* String for header cards */
  char file_name[MAXC]; /* The name of the FITS file to load */
  char *setupname;      /* Name of the optional setup file */
  char plotfile[MAXC];  /* Filename for postscript plot */
  char plotdev[MAXC];   /* PGPLOT display device */
  char line[MAXC];      /* General string for getting input */
  Fits *fits=NULL;      /* Fits file */
  Phdu *phdu=NULL;      /* Primary header-data-unit context descriptor */

  /*
   * Check the command line
   */

  if(argc < 2) {
    fprintf(stderr,"\n");
    fprintf(stderr,
	    "fitshead: prints out the header information for a fits file\n");
    fprintf(stderr,"\nUsage: fitshead [filename]\n\n");
    return 1;
  }

  /*
   * Get the name of the file.
   */

  strcpy(file_name,argv[1]);

  /*
   * Attempt to open the file.
   */

#if 0
  printf("\nReading header from %s.....\n\n",file_name);
#endif
  fits = new_Fits(file_name, 1, 1, 1, 0);
  if(!fits) {
    fprintf(stderr,"\n\nERROR: fitshead.  Could not read fits file\n\n");
    return 1;
  }

  /*
   * Get a pointer to the primary header-data-unit context descriptor.
   */

  phdu = (Phdu *) fits->hdu;

  /*
   * Step through the header, printing out header cards
   */

  for(i=0; i<=phdu->endline; i++) {
    headstr = rheadline(fits,fits->hdu,i);
    if(headstr != NULL) 
      printf("%s\n",headstr);
  }

  /*
   * Clean up and exit
   */

  fits = del_Fits(fits);

  if(no_error) {
    /* printf("\nProgram fitshead.c completed.\n\n"); */
    return 0;
  }
  else {
    fprintf(stderr,"\nERROR. Exiting fitshead.c\n\n");
    return 1;
  }
}


int MAIN_(void)
{
  return 0;
}
