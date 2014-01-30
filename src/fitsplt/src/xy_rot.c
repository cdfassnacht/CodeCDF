/* xy_rot.c
 *
 * Usage: xy_rot fitsfile1 catfile1 fitsfile2 catfile2 ([setup_filename])
 *
 * Description:  Plots two fits files as greyscale plots and then
 *                overlays the catalogs generated from the images.  The
 *                user then can rotate the positions obtained in the
 *                first plot until they match the positions in the second
 *                plot.
 *
 * To compile use the appropriate makefile in this directory.
 *
 * 04Jun2004 CDF,  First working version
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
#include "plotfuncs.h"

/*.......................................................................
 *
 * Function declarations
 *
 */

Secat find_closest(Pos cpos, Secat *secat, int ncat);
void markcatobj(Secat object);
void catxy2dpos(Pos cpos, Secat *secat, int ncat);
void pa2dpos(Pos cpos, double pa, Secat *incat, Secat *rotcat, int ncat);

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
  int ncat1,ncat2;      /* Number of catalog objects */
  float cx,cy;          /* Cursor position */
  double pa;            /* Position angle (N->E) */
  char cchar='0';       /* Character returned by cursor command */
  char fitsfile1[MAXC]; /* The name of the 1st FITS file to load */
  char catfile1[MAXC];  /* The name of the 1st SExtractor catalog file */
  char fitsfile2[MAXC]; /* The name of the 2nd FITS file to load */
  char catfile2[MAXC];  /* The name of the 2nd SExtractor catalog file */
  char setupfile[MAXC]; /* Name of the optional setup file */
  char plotfile[MAXC];  /* Filename for postscript plot */
  char plotdev[MAXC];   /* PGPLOT display device */
  char line[MAXC];      /* General string to read input */
  Pos cpos;             /* Cursor position */
  Image *image1;        /* The container of the 1st FITS image array */
  Image *image2;        /* The container of the 2nd FITS image array */
  Setup *setup1;        /* Container for the 1st image display info */
  Setup *setup2;        /* Container for the 2nd image display info */
  Secat *catdat1=NULL;  /* SExtractor data for 1st catalog */
  Secat *catdat2=NULL;  /* SExtractor data for 2nd catalog */
  Secat *rotcat=NULL;   /* Rotated version of first catalog */
  Secat bestmatch1;     /* Closest match in catalog 1 to the cursor position */
  Secat bestmatch2;     /* Closest match in catalog 2 to the cursor position */

  /*
   * Check the command line
   */

  if(argc < 5) {
    fprintf(stderr,"\nxy_rot fitsfile1 catfile1 fitsfile2 catfile2 ");
    fprintf(stderr,"(setupfile)\n\n");
    return 1;
  }

  /*
   * Get the names of the files.
   */

  strcpy(fitsfile1,argv[1]);
  strcpy(catfile1,argv[2]);
  strcpy(fitsfile2,argv[3]);
  strcpy(catfile2,argv[4]);
  if(argc == 6) {
    usefile = 1;
    strcpy(setupfile,argv[5]);
  }
  else
    sprintf(setupfile,"");

  /*
   * Read the 1st data array.
   */

  if(!(image1 = new_Image(fitsfile1))) {
    fprintf(stderr,"ERROR.  Exiting program. Could not open %s\n\n",
	    fitsfile1);
    return 1;
  }

  if(load_image_data(image1) == ERROR)
    no_error = 0;

  /*
   * Create setup structure and fill with information used to set plotting
   *  parameters.
   */

  if(no_error)
    if(!(setup1 = fill_setup(image1,setupfile,usefile,1)))
      no_error = 0;

  /*
   * Loop over displaying image until happy with results
   */

  if(no_error) {
    open_plot_window();
    if(display_image(image1,setup1))
      no_error = 0;
    else
      while(no_error && another_loop) {
	switch(another_loop = setup_menu(image1,setup1)) {
	case 1:
	  if(display_image(image1,setup1))
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
   * Read the 2nd data array and print out basic image info.
   */

  if(no_error)
    if(!(image2 = new_Image(fitsfile2)))
      no_error = 0;

  if(no_error)
    if(load_image_data(image2) == ERROR)
      no_error = 0;

  /*
   * Create setup structure and fill with information used to set plotting
   *  parameters.
   */

  if(no_error)
    if(!(setup2 = fill_setup(image2,setupfile,usefile,1)))
      no_error = 0;

  /*
   * Loop over displaying image until happy with results
   */

  if(no_error) {
    open_plot_window();
    if(display_image(image1,setup1))
      no_error = 0;
    else
      while(no_error && another_loop) {
	switch(another_loop = setup_menu(image1,setup1)) {
	case 1:
	  if(display_image(image1,setup1))
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
   * Read in the first SExtractor catalog file
   */

  printf("---------------------------------------------------------------\n");
  if(no_error)
    if(!(catdat1 = read_secat(catfile1,'#',&ncat1,2)))
      no_error = 0;

  /*
   * Plot the catalog
   */

  if(no_error) {
    cpgslct(1);
    printf("\nPlotting SExtractor positions for first catalog\n");
    if(plot_secat(setup1,catdat1,ncat1,2,3,2.0))
      no_error = 0;
  }

  /*
   * Read in second SExtractor catalog file
   */

  if(no_error)
    if(!(catdat2 = read_secat(catfile2,'#',&ncat2,2)))
      no_error = 0;

  /*
   * Plot the catalog
   */

  if(no_error) {
    cpgslct(2);
    printf("\nPlotting SExtractor positions for second catalog\n");
    if(plot_secat(setup2,catdat2,ncat2,2,3,2.0))
      no_error = 0;
  }

  /*
   * Get the desired position in image 1
   */

  printf("---------------------------------------------------------------\n");
  if(no_error) {
    cx = 0.0;
    cy = 0.0;
    printf("\n");
    printf("Find a compact object that is clearly detected in both images\n");
    printf("Click the mouse on the object in the FIRST image\n");
    cpgslct(1);
    if(cpgcurs(&cx,&cy,&cchar)>0) {
      printf("\nMouse clicked at: %7.2f %7.2f\n",cx,cy);
      cpos.x = cx;
      cpos.y = cy;
      bestmatch1 = find_closest(cpos,catdat1,ncat1);
      printf("Position of closest match is: %7.2f %7.2f\n",
	     bestmatch1.x,bestmatch1.y);
      markcatobj(bestmatch1);
    }
    else {
      fprintf(stderr,"*** Invalid cursor input ***\n");
      no_error = 0;
    }
  }

  /*
   * Now get position in image 2
   */

  if(no_error) {
    printf("\n Now click the mouse on the object in the SECOND image\n");
    cpgslct(2);
    if(cpgcurs(&cx,&cy,&cchar)>0) {
      printf("\nMouse clicked at: %7.2f %7.2f\n",cx,cy);
      cpos.x = cx;
      cpos.y = cy;
      bestmatch2 = find_closest(cpos,catdat2,ncat2);
      printf("Position of closest match is: %7.2f %7.2f\n\n",
	     bestmatch2.x,bestmatch2.y);
      markcatobj(bestmatch2);
    }
    else {
      fprintf(stderr,"*** Invalid cursor input ***\n");
      no_error = 0;
    }
  }

  /*
   * Calculate x,y offsets from selected position in image 1
   */

  printf("---------------------------------------------------------------\n");
  if(no_error) {
    printf("\nInital estimate of shift between images: %8.2f %8.2f\n",
	   (bestmatch2.x-bestmatch1.x),(bestmatch2.y-bestmatch1.y));
    cpos.x = bestmatch1.x;
    cpos.y = bestmatch1.y;
    catxy2dpos(cpos,catdat1,ncat1);
  }

  /*
   * Create a copy of the first catalog
   */

  if(no_error)
    if(!(rotcat = new_secat(ncat1)))
      no_error = 0;

  /*
   * Redisplay image 2 and loop on PA.  
   */

  pa = 0.0;
  cpos.x = bestmatch2.x;
  cpos.y = bestmatch2.y;
  cpgslct(2);

  printf("---------------------------------------------------------------\n");
  printf("\nRedisplaying second image, but with shifted catalog from first ");
  printf("image.\n");
  while(another_loop && no_error) {
    pa2dpos(cpos,pa,catdat1,rotcat,ncat1);
    if(display_image(image2,setup2))
      no_error = 0;
    markcatobj(bestmatch2);
    if(plot_secat(setup2,rotcat,ncat1,2,3,2.0))
      no_error = 0;
    if(plot_secat(setup2,catdat2,ncat2,2,7,1.5))
      no_error = 0;
    printf("\nGreen:  catalog 1 positions rotated by %6.1f degrees\n",pa);
    printf("Yellow: catalog 2 positions\n");
    printf("Red:    center of rotation\n");
    printf("Enter new pa, or q to quit [%6.1f] ",pa);
    fgets(line,MAXC,stdin);
    if(line[0] == 'q')
      another_loop = 0;
    else if(line[0] != '\n') {
      while(sscanf(line,"%lf",&pa) != 1) {
	fprintf(stderr,"ERROR.  Enter value again: ");
	fgets(line,MAXC,stdin);
      }
    }
  }


  cpgend();

  /*
   * Clean up and exit
   */

  cpgend();
  image1 = del_Image(image1);
  setup1 = del_setup(setup1);
  catdat1 = del_secat(catdat1);
  image2 = del_Image(image2);
  setup2 = del_setup(setup2);
  catdat2 = del_secat(catdat2);

  if(no_error) {
    printf("\nProgram xy_rot.c completed.\n\n");
    return 0;
  }
  else {
    fprintf(stderr,"\nERROR. Exiting xy_rot.c\n\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function find_closest
 *
 * Finds the object giving the closest match between the input postion and 
 *  the members of a SExtractor catalog.
 *  The closest matching positions are returned in a Secat structure.
 *
 * Inputs: Pos cpos            input position
 *         Secat *secat        catalog based on fits image (x,y) positions
 *         int ncat            number of members in secat
 *
 * Output: Secat bestmatch     best match
 *
 */

Secat find_closest(Pos cpos, Secat *secat, int ncat)
{
  int i,j;                /* Looping variables */
  double dx,dy,dpos;      /* Offsets between (x,y) pairs */
  double mindpos;         /* Running minimum position offset */
  Secat bestmatch;        /* Container for best-matching xy object */
  Secat *sptr;            /* Pointers to navigate astcat and xycat */

  /*
   * Loop through catalog
   */

  dx = cpos.x - secat->x;
  dy = cpos.y - secat->y;
  bestmatch.dpos = sqrt(dx*dx + dy*dy);
  for(i=0,sptr=secat; i < ncat; i++,sptr++) {
    dx = cpos.x - sptr->x;
    dy = cpos.y - sptr->y;
    if((dpos = sqrt(dx*dx + dy*dy)) < bestmatch.dpos) {
      bestmatch = *sptr;
      bestmatch.dpos = dpos;
    }
  }
  printf("For this position, the best match is ID=%d at",
	 bestmatch.id);
  printf(" mindpos = %6.3f pix\n",bestmatch.dpos);

  /*
   * Exit
   */

  return bestmatch;
}

/*.......................................................................
 *
 * Function markcatobj
 *
 * Marks position of object contained in the passed secat structure.
 *
 * Input:  Secat object        object to be marked
 *
 * Output: (none)
 *
 */

void markcatobj(Secat object)
{
  cpgsci(2);
  cpgsfs(2);
  cpgslw(3);
  cpgcirc(object.x,object.y,3.0*object.fwhm);
  cpgsci(1);
  cpgslw(1);
}

/*.......................................................................
 *
 * Function catxy2dpos
 *
 * Calculates offsets from some central position in (x,y) space.  Converts
 *  the offsets into (r,theta) format and puts results into the dpos and
 *  dpostheta parameters in the input catalog.
 *
 * Inputs: Pos cpos            central position
 *         Secat *secat        catalog based on fits image (x,y) positions
 *         int ncat            number of members in secat
 *
 * Output: (none)
 *
 */

void catxy2dpos(Pos cpos, Secat *secat, int ncat)
{
  int i,j;                /* Looping variables */
  double dx,dy,dpos;      /* Offsets between (x,y) pairs */
  double mindpos;         /* Running minimum position offset */
  Pos tmppos;             /* Temporary container for (x,y) offset */
  Secat *sptr;            /* Pointers to navigate astcat and xycat */

  /*
   * Loop through catalog
   */

  for(i=0,sptr=secat; i < ncat; i++,sptr++) {
    sptr->dx = sptr->x - cpos.x;
    sptr->dy = sptr->y - cpos.y;
    tmppos.x = sptr->dx;
    tmppos.y = sptr->dy;
    xy2rth(tmppos,&sptr->dpos,&sptr->dpostheta,1);
  }
}

/*.......................................................................
 *
 * Function pa2dpos
 *
 * Converts a central position, a PA, and catalog offsets into rotated
 *  x and y positions.
 *
 * Inputs: Pos cpos            x,y position of center
 *         double pa           PA to rotate the catalog positions by
 *         Secat *incat        input catalog containing offsets
 *         Secat *rotcat       output catalog with rotated positions
 *         int ncat            number of catalog members
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

void pa2dpos(Pos cpos, double pa, Secat *incat, Secat *rotcat, int ncat)
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
    rptr->x = cpos.x + tmppos.x;
    rptr->y = cpos.y + tmppos.y;
  }
}

int MAIN_(void)
{
  return 0;
}
