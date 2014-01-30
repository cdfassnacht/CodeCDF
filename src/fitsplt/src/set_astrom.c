/* set_astrom.c
 *
 * Usage: set_astrom fitsfile astrofile xycatfile ([setup_filename])
 *
 * Description:  Plots a FITS file as a greyscale plot and then overlays
 *                the positions of a starlist (e.g., USNO stars) at a
 *                given angle.  The user can then interactively change the 
 *                PA and location of the starlist.  The output is a
 *                file that can be used as input to the iraf task ccmap.
 *
 * Inputs:
 *  1. FITS file for which astrometry needs to be done.
 *  2. SExtractor file created from the FITS image.
 *  3. File containing (approximate) RA and Dec for the center of the field.
 *     This file should be in a ?? format.
 *  4. Input astrometric list.  This can be either in the distcalc
 *     (ID RA(hms) Dec(hms); format=5), USNO (RA(hms) Dec(hms); format=??),
 *     or 2MASS format (RA(degs) Dec(degs); format=4).
 *  5. Format flag
 *
 * To compile use the appropriate makefile in this directory.
 *
 * 26Jul2004 CDF,  A modification of astrom_rot.c
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

Secat *make_astcat(char *posfile, char *astfile, int *nast, int format);
void pa2dpos(double centx, double centy, double pa, Setup *setup, Secat *incat, 
	     Secat *rotcat, int ncat);
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
  int format;           /* Format of astrometric file */
  int another_loop=1;   /* Flag set to 0 to stop loop of PA changes */
  int nast;             /* Number of astrometric stars */
  int ncat;             /* Number of catalog objects */
  int nmatch=0;         /* Number of xy->astcat source matches */
  double centx,centy;   /* Central position */
  double pa;            /* PA */
  char fitsfile[MAXC];  /* The name of the FITS file to load */
  char catfile[MAXC];   /* The name of the SExtractor catalog file */
  char posfile[MAXC];   /* File containing the central position*/
  char astfile[MAXC];   /* The name of the astrometric file */
  char setupfile[MAXC]; /* Name of the optional setup file */
  char line[MAXC];      /* General string for getting input */
  Image *image;         /* The container of the FITS image array */
  Setup *setup;         /* Container for the image display info */
  Secat *astcat=NULL;   /* Astrometric data */
  Secat *catdat=NULL;   /* SExtractor catalog data */
  Secat *rotcat=NULL;   /* A working copy of catdat */
  Secat *matchcat=NULL; /* Catalog of matches */

  /*
   * Check the command line
   */

  if(argc < 6) {
    fprintf(stderr,"\nset_astrom requires at least five inputs:\n");
    fprintf(stderr,"  1. Input fits file\n");
    fprintf(stderr,"  2. SExtractor file created from input fits file\n");
    fprintf(stderr,"  3. File containing approximate RA and Dec of field ");
    fprintf(stderr,"center\n");
    fprintf(stderr,"  4. Astrometric file\n");
    fprintf(stderr,"  5. Astrometric file format\n");
    fprintf(stderr,"  6. Optional setup file\n\n");
    fprintf(stderr," Options for format:\n");
    fprintf(stderr,"  4 ==> \"RA(degrees) Dec(degrees)\" (2MASS format)\n");
    fprintf(stderr,"  5 ==> \"ID RAhr RAmin RAsec Decdeg Decmin Decsec\"");
    fprintf(stderr,"  (distcalc output format)\n");
    return 1;
  }

  /*
   * Get the names of the files.
   */

  strcpy(fitsfile,argv[1]);
  strcpy(catfile,argv[2]);
  strcpy(posfile,argv[3]);
  strcpy(astfile,argv[4]);
  if(sscanf(argv[5],"%d",&format) != 1) {
    fprintf(stderr,"ERROR: Bad value for format (input 5)\n\n");
    return 1;
  }
  if(argc == 7) {
    usefile = 1;
    strcpy(setupfile,argv[6]);
  }
  else
    sprintf(setupfile,"");

  /*
   * Read in the fits file
   */

  if(!(image = new_Image(fitsfile))) {
    fprintf(stderr,"ERROR.  Exiting program\n\n");
    return 1;
  }

  /*
   * Print out some basic image info
   */

  if(image_info(image))
    no_error = 0;

  /*
   * Create setup structure and fill with information used to set plotting
   *  parameters.
   */

  if(no_error)
    if(!(setup = fill_setup(image,setupfile,usefile,1)))
      no_error = 0;

  /*
   * Read the SExtractor file created from the input FITS file
   */

  if(no_error)
    if(!(catdat = read_secat(catfile,'#',&ncat,2)))
      no_error = 0;

  /*
   * Calculate the offsets from the central position
   */

  if(no_error)
    if(!(astcat = make_astcat(posfile,astfile,&nast,format)))
      no_error = 0;

  /*
   * Create a copy of the astrometric catalog
   */

  if(no_error)
    if(!(rotcat = new_secat(nast)))
      no_error = 0;

  /*
   * Set defaults for first round through display
   */

  pa = 0.0;
  centx = setup->centx;
  centy = setup->centy;
  printf("For first run setting central position to %7.1f %7.1f\n",
	 centx,centy);
  if(no_error)
    if(cpgopen("/xs") <= 0)
      no_error = 0;

  /*
   * Display the image and loop on PA.  
   */

  while(another_loop && no_error) {
    pa2dpos(centx,centy,pa,setup,astcat,rotcat,nast);
    if(display_image(image,setup))
      no_error = 0;
    if(plot_secat(setup,rotcat,nast,2,3,2.0))
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
#if 0
  /*
   * Print output, with finalized PA choice, to a cl script that
   *  will be used in iraf.
   */

  if(no_error)
    if(print_clscript(setup,astcat,fitsfile,"tmp.cl"))
      no_error = 0;

  /*
   * Find matches 
   */

  if(no_error) {
    if(!(matchcat = find_match(setup,astcat,catdat,ncat,&nmatch)))
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
#endif
  /*
   * Clean up and exit
   */

  cpgend();
  image = del_Image(image);
  setup = del_setup(setup);
  astcat = del_secat(astcat);
  catdat = del_secat(catdat);
  matchcat = del_secat(matchcat);

  if(no_error) {
    printf("\nProgram set_astrom.c completed.\n\n");
    return 0;
  }
  else {
    fprintf(stderr,"\nERROR. Exiting set_astrom.c\n\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function make_astcat
 *

 * Takes two input files, the first containing a central position and
 *  the second containing RA and Dec, and computes the offsets from
 *  the central position.  Then the functions places those offsets into
 *  a Secat structure, along with the input RAs and Decs, and returns
 *  the structure.
 *
 * Inputs: char *posfile       file containing central position
 *         char *astfile       file containing astrometric RA/Dec positions
 *         int *nast           number of astrometric positions/offsets (set
 *                              by this function)
 *
 * Output: Secat *astcat       astrometric catalog
 *
 */

Secat *make_astcat(char *posfile, char *astfile, int *nast, int format)
{
  int i;                /* Looping variable */
  int no_error=1;       /* Flag set to 0 on error */
  int ncent;            /* Number of lines in posfile */
  Pos *offsets=NULL;    /* Array to hold offsets */
  Pos *optr;            /* Pointer to navigate offsets */
  Skypos *skypos=NULL;  /* Array of sky positions */
  Skypos *skptr;        /* Pointer to navigate skypos */
  Secat *centpos=NULL;  /* Central position for distance calculations */
  Secat *astcat=NULL;   /* New astrometric catalog */
  Secat *sptr;          /* Pointer to navigate astcat */

  /*
   * Read in the file containing the central postion, which will be in the 
   *  format of a distcalc file.
   */

  if(!(centpos = read_distcalc(posfile,'#',&ncent,0)))
    no_error = 0;

  /*
   * Read the astrometric file 
   */

  if(no_error)
    if(!(astcat = read_secat(astfile,'#',nast,format)))
      no_error = 0;

  /*
   * Allocate memory for the skypos array
   */

  if(no_error)
    if(!(skypos = new_skypos(*nast)))
      no_error = 0;

  /*
   * Load the skypos values into the skypos array
   */

  if(no_error)
    for(i=0,sptr=astcat,skptr=skypos; i< *nast; i++,sptr++,skptr++) {
      *skptr = sptr->skypos;
      sprintf(skptr->label,"%04d",sptr->id);
    }

  /*
   * Calculate the offsets
   */

  if(no_error)
    if(!(offsets = dspos2xy(centpos->skypos,skypos,*nast)))
      no_error = 0;

  /*
   * Transfer the offsets to the astrometric catalog
   */

  for(i=0,sptr=astcat,optr=offsets; i< *nast; i++,sptr++,optr++) {
    sptr->dx = optr->x;
    sptr->dy = optr->y;
    sptr->fwhm = 10;
    xy2rth(*optr,&sptr->dpos,&sptr->dpostheta,1);
  }

  /*
   * Clean up and return
   */

  offsets = del_pos(offsets);
  skypos = del_skypos(skypos);
  centpos = del_secat(centpos);

  if(no_error)
    return astcat;
  else {
    fprintf(stderr,"ERROR: make_astcat\n");
    return del_secat(astcat);
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
 *         Secat *incat        input catalog containing offsets
 *         Secat *rotcat       output catalog with rotated positions
 *         int ncat            number of catalog members
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

void pa2dpos(double centx, double centy, double pa, Setup *setup, Secat *incat, 
	     Secat *rotcat, int ncat)
{
  int i;              /* Looping variables */
  double dpostheta;   /* New theta position */
  Pos tmppos;         /* Temporary holder for x,y position */
  Secat *iptr,*rptr;  /* Pointers to navigate incat and rotcat */

  /*
   * Loop through, calculating new offsets given the PA
   */

  for(i=0,iptr=incat,rptr=rotcat; i<ncat; i++,iptr++,rptr++) {
    *rptr = *iptr;
    dpostheta = iptr->dpostheta + pa;
    rth2xy(iptr->dpos,dpostheta,&tmppos,1);
    rptr->x = centx + tmppos.x / setup->pixscale;
    rptr->y = centy + tmppos.y / setup->pixscale;
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
#if 0
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
    strcpy(lptr->text,star);
  }

  return 0;
}
#endif
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
	if(bestmatch.dpos < matchlim)
	  break;
      }
      printf("find_match: For %s, mindpos = %6.3f pix\n",lptr->text,
	     bestmatch.dpos);
      *mptr = bestmatch;
      mptr->skypos = aptr->skypos;
      mptr->dx = aptr->dx;
      mptr->dy = aptr->dy;
      strcpy(mptr->name,lptr->text);
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
