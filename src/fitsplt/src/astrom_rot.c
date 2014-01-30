/* astrom_rot.c
 *
 * Usage: astrom_rot fitsfile astrofile xycatfile ([setup_filename])
 *
 * Description:  Plots a FITS file as a greyscale plot and then overlays
 *                the positions of a starlist (e.g., USNO stars) at a
 *                given angle.  The user can then interactively rotate the 
 *                PA of the starlist to get a first estimate of the PA
 *                of the image.
 *
 * To compile use the appropriate makefile in this directory.
 *
 * 17Jul2003 CDF,  A modification of fitsplt.c
 * v20Jul2003 CDF, Added the new find_match and print_astfile functions to 
 *                  find x,y coordinates of closest matches in xycatfile to
 *                  the rotated offsets in astrofile, and then print output
 *                  that can be used directly as input for either Marc Postman's
 *                  astrometry program (format=2) or IRAF's ccmap (format=1).
 * v11Aug2003 CDF, Added a printout to a file (usno_good.list) that can be
 *                  used as a better basis for the [source]_usno_good.[epoch]
 *                  file than the usual [source]_source_star.offset file since
 *                  the new output file only contains sources that lie within
 *                  the input fits image area.
 * v06May2004 CDF, Added a maximum pixel offset for a "good" match
 * v07Jun2004 CDF, Moved setup creation and filling into the new fill_setup
 *                  function.
 *                 Fixed displaying functions to match the new version of
 *                  display_image.
 * v2008Jul16 CDF, Small modifications to incorporate new image reading
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

#define SIZE 30

/*.......................................................................
 *
 * Function declarations
 *
 */

int pa2dpos(double centx, double centy, double pa, Secat *secat, Setup *setup);
int print_clscript(Setup *setup, Secat *astcat, char *fitsname, char *clname);
Secat *find_match(Setup *setup, Secat *astcat, Secat *xycat, int nxy,
		  int *nmatch);
int print_astfile(Setup *setup, Secat *matchcat, int nmatch, char *pfname,
		  int format);


/*.......................................................................
 *
 * Main program
 *
 */

int main(int argc, char *argv[])
{
  int no_error=1;       /* Flag set to 0 on error */
  int usefile=0;        /* Flag set to 1 to use input setup file */
  int another_loop=1;   /* Flag set to 0 to stop loop of PA changes */
  int screenflag=1;     /* Flag set to 0 for no screen output */
  int fileflag=1;       /* Flag set to 0 for no screen output */
  int nast;             /* Number of astrometric stars */
  int ncat;             /* Number of catalog objects */
  int nmatch=0;         /* Number of xy->astcat source matches */
  double centx,centy;   /* Lens position */
  double pa;            /* PA */
  char fitsfile[MAXC];  /* The name of the FITS file to load */
  char astfile[MAXC];   /* The name of the distcalc astrometric file */
  char catfile[MAXC];   /* The name of the SExtractor catalog file */
  char setupfile[MAXC]; /* Name of the optional setup file */
  char plotfile[MAXC];  /* Filename for postscript plot */
  char *plotdev;        /* PGPLOT display device */
  char line[MAXC];      /* General string for getting input */
  Image *image;         /* The container of the FITS image array */
  Setup *setup;         /* Container for the image display info */
  Secat *astdat=NULL;   /* Astrometric data */
  Secat *catdat=NULL;   /* SExtractor catalog data */
  Secat *matchcat=NULL; /* Catalog of matches */
  Labinfo *lptr;

  /*
   * Check the command line
   */

  if(argc < 4) {
    fprintf(stderr,"\nastrom_rot fitsfile astrofile xycatfile (setupfile)\n\n");
    return 1;
  }

  /*
   * Get the names of the files.
   */

  strcpy(fitsfile,argv[1]);
  strcpy(astfile,argv[2]);
  strcpy(catfile,argv[3]);
  if(argc == 5) {
    usefile = 1;
    strcpy(setupfile,argv[3]);
  }
  else
    sprintf(setupfile,"");

  /*
   * Read the image header and print out some basic info.
   */

  if(!(image = new_Image(fitsfile))) {
    fprintf(stderr,"ERROR.  Could not open %s.  Exiting program\n\n",
	    fitsfile);
    return 1;
  }

  /*
   * Actually load the image data
   */

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
  another_loop = 1;

  /*
   * Read in astrometric file which has been output from distcalc
   */

  if(no_error)
    if(!(astdat = read_distcalc(astfile,'#',&nast,1)))
      no_error = 0;
  printf("%d\n",nast);


  /*
   * Read in SExtractor catalog file
   */

  if(no_error)
    if(!(catdat = read_secat(catfile,'#',&ncat,7)))
      no_error = 0;

  /*
   * Allocate memory for circle array in setup structure
   */

  if(no_error)
    if(!(setup->circle = new_labinfo(nast)))
      no_error = 0;

  /*
   * Get lens position
   */

  if(no_error) {
    centx = catdat->x;
    centy = catdat->y;
    printf("Enter lens position [x y]: [%8.2f %8.2f] ",centx,centy);
    fgets(line,MAXC,stdin);
    if(line[0] != '\n') {
      while(sscanf(line,"%lf %lf",&centx,&centy) != 2) {
	fprintf(stderr,"ERROR: Enter lens position again: ");
	fgets(line,MAXC,stdin);
      }
    }
    setup->drawmarker = TRUE;
    setup->markx = (float) centx;
    setup->marky = (float) centy;
  }

  /*
   * Set up the circles
   */

  if(no_error) {
    setup->ncircle = nast;
  }  

  /*
   * Display the image and loop on PA.  
   */

  pa = 0.0;
  if(no_error)
    open_plot_window();

  while(another_loop && no_error) {
    if(pa2dpos(centx,centy,pa,astdat,setup))
      no_error = 0;
    if(display_image(image,setup))
      no_error = 0;
    printf("\nEnter new pa, or q to quit [%6.1f] ",pa);
    fgets(line,MAXC,stdin);
    printf("%s",line);
    if(line[0] == 'q')
      another_loop = 0;
    else if(line[0] != '\n') {
      while(sscanf(line,"%lf",&pa) != 1) {
	fprintf(stderr,"ERROR.  Enter value again: ");
	fgets(line,MAXC,stdin);
      }
      printf("pa = %f\n",pa);
    }
  }

  /*
   * Print output, with finalized PA choice, to a cl script that
   *  will be used in iraf.
   */

  if(no_error)
    if(print_clscript(setup,astdat,fitsfile,"tmp.cl"))
      no_error = 0;

  /*
   * Find matches 
   */

  if(no_error) {
    if(!(matchcat = find_match(setup,astdat,catdat,ncat,&nmatch)))
      no_error = 0;
    else
      printf("\nFound %d matches.\n",nmatch);
  }

  /*
   * Print also to an output file that can be used directly as
   *  an input file for either Mark Postman's astrometry program or
   *  IRAF's ccmap.
   */

  if(no_error)
    if(print_astfile(setup,matchcat,nmatch,"astrom.dat",1))
      no_error = 0;

  /*
   * Clean up and exit
   */

  cpgend();
  image = del_Image(image);
  setup = del_setup(setup);
  astdat = del_secat(astdat);
  catdat = del_secat(catdat);
  matchcat = del_secat(matchcat);

  if(no_error) {
    printf("\nProgram astrom_rot.c completed.\n\n");
    return 0;
  }
  else {
    fprintf(stderr,"\nERROR. Exiting astrom_rot.c\n\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function pa2dpos
 *
 * Converts a central position, a PA, and catalog offsets into rotated
 *  x and y positions.
 *
 * Inputs: double centx        x position of center
 *         double centy        y position of center
 *         double pa           PA to rotate the catalog positions by
 *         Secat *secat        catalog containing offsets
 *         Setup *setup        setup structure -- also contains the "circle"
 *                              array which is used to mark the rotated
 *                              positions.
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int pa2dpos(double centx, double centy, double pa, Secat *secat, Setup *setup)
{
  int i;             /* Looping variables */
  float r_circ;      /* Circle radius in pixels */
  double dpostheta;  /* New theta position */
  char *star;        /* Star number, without the root name */
  Pos tmppos;        /* Temporary holder for x,y position */
  Labinfo *lptr;     /* Pointer to navigate Setup->circle */
  Secat *sptr;       /* Pointer to navigate secat */

  /*
   * Set circle radius to be a fraction of image size
   */

  r_circ = 0.015 * setup->ysize;

  /*
   * Loop through, calculating new offsets given the PA
   */

  printf("centx = %f, centy = %f\n",centx,centy);
  for(i=0,sptr=secat,lptr=setup->circle; i<setup->ncircle; i++,sptr++,lptr++) {
    dpostheta = sptr->dpostheta + pa;
    rth2xy(sptr->dpos,dpostheta,&tmppos,1);
    lptr->x = centx + (tmppos.x / setup->pixscale);
    lptr->y = centy + (tmppos.y / setup->pixscale);
    lptr->size = r_circ;
    if(! (star = strchr(sptr->name,'_')))
      star = sptr->name;
    else
      star++;
    /* strcpy(lptr->text,star); */
  }

  return 0;
}

/*.......................................................................
 *
 * Function print_clscript
 *
 * Prints out a cl script that can be used in IRAF to examine all, say,
 *  USNO stars in an image.  The user then determines which USNO stars
 *  are good, i.e., not saturated and not galaxies.
 *
 * Inputs: Setup *setup        setup structure, contains (x,y) positions
 *                              of stars
 *         Secat *astcat       catalog containing (RA,Dec) positions of stars
 *         char *fitsname      name of fits file to be imexam'ed
 *         char *clname        name of output cl script
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int print_clscript(Setup *setup, Secat *astcat, char *fitsname, char *clname)
{
  int i;           /* Looping variable */
  int xlo,xhi;     /* x limits of postage stamp */
  int ylo,yhi;     /* y limits of postage stamp */
  char root[10];   /* Root name of object */
  char *star;      /* Star number, without the root name */
  Labinfo *lptr;   /* Pointer to navigate setup->circle */
  Secat *aptr;     /* Pointer to navigate astcat */
  FILE *ofp=NULL;  /* Pointer to output cl file */

  /*
   * Set root name
   */

  strncpy(root,setup->circle->text,4);
  root[4] ='\0';

  /*
   * Open output file
   */

  printf("\n");
  if(!(ofp = open_writefile(clname))) {
    fprintf(stderr,"ERROR: print_clscript\n");
    return 1;
  }

  /*
   * Print header info to output cl file
   */

  fprintf(ofp,"int fooint\n");
  fprintf(ofp,"real x,y,fwhm,junk\n");
  fprintf(ofp,"struct *coordlist\n");
  fprintf(ofp,"del(\"foo.list\")\n");

  /*
   * Loop through structures and print out 
   */

  for(i=0,lptr=setup->circle,aptr=astcat; i<setup->ncircle; i++,lptr++,aptr++) {

    /*
     * Set postage stamp limits
     */

    xlo = (int) (lptr->x - SIZE);
    xhi = (int) (lptr->x + SIZE);
    ylo = (int) (lptr->y - SIZE);
    yhi = (int) (lptr->y + SIZE);
    if(! (star = strchr(lptr->text,'_')))
      star = lptr->text;
    else
      star++;
    
    /*
     * Print out to cl script for boxes within the main image
     */

    if(xlo > 0 && xhi <= setup->xsize && ylo > 0 && yhi <= setup->ysize) {
      fprintf(ofp,"print(\"\")\n");
      fprintf(ofp,
	      "print(\"*************************************************\")\n");
      fprintf(ofp,"print(\"Displayed star is %s.\")\n",lptr->text);
      fprintf(ofp,"imexam %s[%d:%d,%d:%d] ",fitsname,xlo,xhi,ylo,yhi);
      fprintf(ofp,"keeplog+ logfil=foo.list\n");
      fprintf(ofp,"coordlist = \"foo.list\"\n");
      fprintf(ofp,"fooint = fscan(coordlist,x,y)\n");
      fprintf(ofp,"fooint = fscan(coordlist,x,y)\n");
      fprintf(ofp,"fooint = fscan(coordlist,x,y,junk,junk,junk,junk,junk,");
      fprintf(ofp,"junk,junk,junk,junk,junk,junk,junk,fwhm)\n");
      fprintf(ofp,"printf(\"%%7.2f %%7.2f %02d %02d %07.4f %+03d %02d %06.3f",
	      aptr->skypos.hr,aptr->skypos.min,aptr->skypos.sec,
	      aptr->skypos.deg,aptr->skypos.amin,aptr->skypos.asec);
      fprintf(ofp," 0.0 0.0 %s %%5.2f\",x,y,fwhm)\n",lptr->text);
      fprintf(ofp,"print(\"\")\n");
      fprintf(ofp,"printf(\"%%7.2f %%7.2f %02d %02d %07.4f %+03d %02d %06.3f",
	      aptr->skypos.hr,aptr->skypos.min,aptr->skypos.sec,
	      aptr->skypos.deg,aptr->skypos.amin,aptr->skypos.asec);
      fprintf(ofp," 0.0 0.0 %s %%5.2f\",x,y,fwhm, >> \"%s_astrom.dat\")\n",
	      lptr->text,root);
      fprintf(ofp,"print(\"\", >> \"%s_astrom.dat\")\n",root);
      fprintf(ofp,"del(\"foo.list\")\n");
    }
  }

  /*
   * Put trailer lines into cl file
   */

  fprintf(ofp,"print(\"\")\n");
  fprintf(ofp,"print(\"Finished.  Output in %s_astrom.dat\")\n",root);
  fprintf(ofp,"print(\"\")\n");
  printf("print_clscript: Finished printing to %s.\n\n",clname);

  /*
   * Clean up and exit
   */

  if(ofp)
    fclose(ofp);
  return 0;

}

/*.......................................................................
 *
 * Function find_match
 *
 * Finds matches between the possibly-rotated input astrometric catalog and 
 *  the xy file constructed by running SExtractor on the input image.
 *  The closest matching positions are returned in a Secat array.
 *
 * Inputs: Setup *setup        setup structure containing (x,y) predicted
 *                              positions of stars
 *         Secat *astcat       astrometric catalog containing (RA,Dec)
 *                              positions of stars
 *         Secat *xycat        catalog based on fits image (x,y) positions
 *         int nxy             number of members in xycat
 *         int *nmatch         number of xycat memberst that match the
 *                              predicted positions (set by this function)
 *
 * Output: Secat *matchcat     catalog containing matched (x,y) positions,
 *                              plus corresponding (RA,Dec) coordinates.
 *
 */

Secat *find_match(Setup *setup, Secat *astcat, Secat *xycat, int nxy,
		  int *nmatch)
{
  int i,j;                /* Looping variables */
  double dx,dy,dpos;      /* Offsets between (x,y) pairs */
  double mindpos;         /* Running minimum position offset */
  double matchlim=3.0;    /* Maximum pixel offset for a "good" match */
  Secat bestmatch;        /* Container for best-matching xy object */
  Secat *matchcat=NULL;   /* Container for the matches */
  Secat *aptr,*xptr;      /* Pointers to navigate astcat and xycat */
  Secat *mptr;            /* Pointer to navigate matchcat */
  Labinfo *lptr;          /* Pointer to navigate setup->circle */

  /*
   * Allocate memory for output structure
   */

  if(!(matchcat = new_secat(setup->ncircle))) {
    fprintf(stderr,"ERROR: find_match.\n");
    return NULL;
  }

  /*
   * Loop through astrometric stars
   */

  mptr = matchcat;
  for(i=0,lptr=setup->circle,aptr=astcat; i < setup->ncircle;
      i++,lptr++,aptr++) {

    /*
     * Only search for a match if the astrometric object actually falls
     *  within the observed image.
     */

    if(lptr->x >= 1 && lptr->x <= (setup->xsize - 1) && 
       lptr->y >= 1 && lptr->y <= (setup->ysize - 1)) {
    
      bestmatch.dpos = 1.0 * setup->xsize;
      for(j=0,xptr=xycat; j<nxy; j++,xptr++) {
	dx = lptr->x - xptr->x;
	dy = lptr->y - xptr->y;
	if((dpos = sqrt(dx*dx + dy*dy)) < bestmatch.dpos) {
	  bestmatch = *xptr;
	  bestmatch.dpos = dpos;
	}
      }
      printf("find_match: For %s, mindpos = %6.3f pix\n",lptr->text,
	     bestmatch.dpos);
      *mptr = bestmatch;
      mptr->skypos = aptr->skypos;
      mptr->dx = aptr->dx;
      mptr->dy = aptr->dy;
      strcpy(mptr->name,aptr->name);
      mptr++;
      (*nmatch)++;
    }
  }

  /*
   * Exit
   */

  return matchcat;
}

/*.......................................................................
 *
 * Function print_astfile
 *
 * Prints out (x,y) and (RA,Dec) pairs in a format appropriate for input
 *  to IRAF's ccmap routine.
 * Also prints out a distcalc-like output called usno_good.list which
 *  only contains the objects from the astrometric data file which lie
 *  within the boundaries of the fits image.
 *
 * Inputs:
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 * v06May2004 CDF, Added a maximum pixel offset for a "good" match
 */

int print_astfile(Setup *setup, Secat *matchcat, int nmatch, char *pfname,
		  int format)
{
  int i;            /* Looping variable */
  int no_error=1;   /* Flag set to 0 on error */
  float maxoffset;  /* Maximum offset */
  char ew[5];       /* Direction of delta_RA */
  char ns[5];       /* Direction of delta_Dec */
  char dname[MAXC]; /* Name of distcalc-like output file */
  char line[MAXC];  /* General string for reading input */
  Skypos skypos;    /* Temporary container for skypos info */
  Secat *mptr;      /* Pointer to navigate matchcat */
  FILE *ofp=NULL;   /* Pointer to output file for ccmap */
  FILE *dfp=NULL;   /* Pointer to output file with distcalc-like format */

  /*
   * Open output files
   */

  sprintf(dname,"usno_good.list");
  printf("\n");
  if(!(ofp = open_writefile(pfname))) {
    fprintf(stderr,"ERROR: print_astfile\n");
    return 1;
  }
  printf("\n");
  if(!(dfp = open_writefile(dname)))
    no_error = 0;

  /*
   * Get maximum offset
   */

  printf("\n");
  printf("Enter the maximum value of mindpos (from find_match) allowed ");
  printf("for a good match (in pixels): ");
  fgets(line,MAXC,stdin);
  while(sscanf(line,"%f",&maxoffset) != 1) {
    fprintf(stderr,"Invalid number.  Enter value again: ");
    fgets(line,MAXC,stdin);
  }

  /*
   * Loop through structures and print out 
   */

  if(no_error) {
    for(i=0,mptr=matchcat; i<nmatch; i++,mptr++) {

      /*
       * Print out objects to output file if they lie within the main image
       *  and have a dpos less than maxoffset
       */

      if(mptr->x >= 1 && mptr->x <= (setup->xsize - 1) && 
	 mptr->y >= 1 && mptr->y <= (setup->ysize - 1) &&
	 mptr->dpos < maxoffset) {

	/*
	 * Print out for astfile
	 */

	skypos = mptr->skypos;

	if(format == 2) {
	  fprintf(ofp,"%7.2f %7.2f %02d %02d %07.4f %+03d %02d %06.3f ",
		  mptr->x,mptr->y,
		  skypos.hr,skypos.min,skypos.sec,
		  skypos.deg,skypos.amin,skypos.asec);
	  fprintf(ofp,"0.0 0.0 %s %04d  %6.2f\n",mptr->name,mptr->id,
		  mptr->fwhm);
	}
	else {
	  fprintf(ofp,"%7.2f %7.2f %02d:%02d:%07.4f %+03d:%02d:%06.3f ",
		  mptr->x,mptr->y,
		  skypos.hr,skypos.min,skypos.sec,
		  skypos.deg,skypos.amin,skypos.asec);
	  fprintf(ofp,"0.0 0.0 %s %04d  %6.2f\n",mptr->name,mptr->id,
		  mptr->fwhm);
	}

	/*
	 * Determine the direction of the shift for distcalc file
	 */

	if(mptr->dx >= 0.0)
	  sprintf(ew,"E");
	else
	  sprintf(ew,"W");
	if(mptr->dy >= 0.0)
	  sprintf(ns,"N");
	else
	  sprintf(ns,"S");

	/*
	 * Print out to distcalc-like file
	 */

	fprintf(dfp,"%11s  %02d %02d %07.4f  %+03d %02d %06.3f ",
		mptr->name,skypos.hr,skypos.min,skypos.sec,
		skypos.deg,skypos.amin,skypos.asec);
	fprintf(dfp,"%10.4f %1s %10.4f %1s %10.4f\n",
		fabs(mptr->dx),ew,fabs(mptr->dy),ns,
		sqrt(mptr->dx * mptr->dx + mptr->dy * mptr->dy));
      }
    }
  }

  /*
   * Clean up and exit
   */

  if(ofp)
    fclose(ofp);
  if(dfp)
    fclose(dfp);

  if(no_error) {
    printf("print_astfile: Finished printing to %s and %s.\n",pfname,
	   dname);
    return 0;
  }
  else {
    printf("ERROR: print_astfile\n");
    return 1;
  }
}


int MAIN_(void)
{
  return 0;
}
