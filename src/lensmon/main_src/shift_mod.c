/*.......................................................................
 *
 * shift_mod.c
 *
 * Usage: shift_mod [source_rootname] [filename_extension]
 *
 * Basically a one-off program to shift difmap models which have been
 *  created with respect to one pointing center to be defined with respect
 *  to a new pointing center.
 *
 * 16Apr97 CDF
 * v02Sep97 CDF, Included dataio and structdef libraries.
 * v23May98 CDF, Put header info into output file.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include "dataio.h"
#include "structdef.h"
#include "modfuncs.h"
#include "coords.h"

/*.......................................................................
 *
 * Function declarations
 *
 */

Pos *cinfo2xy(Compinfo *cinfo, int nlines);
Compinfo *xy2cinfo(Compinfo *old, Pos *newpos, int nlines);
int get_newcent(Skypos *newcent);
int print_modfile(Compinfo *cinfo, int nlines, char *outname, Skypos newcent);

/*.......................................................................
 *
 * Main program
 *
 */

int main(int argc, char *argv[])
{
  int no_error=1;           /* Flag set to 0 on error */
  int nlines;               /* Number of lines in old model file */
  char *root;               /* Source root name */
  char *ext;                /* File name extension */
  char modname[MAXC];       /* Model file name */
  char outname[MAXC];       /* Output filename */
  Compinfo *oldmod=NULL;    /* Old model component info */
  Compinfo *newmod=NULL;    /* New model component info */
  Pos *oldxy=NULL;          /* Old x,y offsets */
  Pos *newxy=NULL;          /* New x,y offsets */
  Skypos oldcent;           /* Old model center (RA, Dec) */
  Skypos newcent;           /* New model center (RA, Dec) */
  Skypos *modpos=NULL;      /* RA and Dec of model components */
  FILE *mfp=NULL;           /* Old model file */

  /*
   * Check command line
   */

  if(argc != 3) {
    fprintf(stderr,"\nUsage: shift_mod [source_rootname] ");
    fprintf(stderr,"[filename_extension]\n\n");
    return 1;
  }

  /*
   * Open the model file
   */

  root = argv[1];
  ext = argv[2];

  get_modname(modname,root,ext,0);

  if(!(mfp = open_modfile(modname))) {
    fprintf(stderr,"ERROR: Exiting program\n");
    return 1;
  }

  /*
   * Count the number of lines in the input file
   */

  if((nlines = n_lines(mfp,'!')) < 0)
    no_error = 0;
  else
    rewind(mfp);

  /*
   * Put the model info into a Compinfo array
   */

  if(no_error)
    if(!(oldmod = fill_cinfo(mfp,nlines)))
      no_error = 0;

  /*
   * Create a list of (x,y) offsets (in arcsec) from the model center
   *  by using the (r,theta) positions in the model data file
   */

  if(no_error)
    if(!(oldxy = cinfo2xy(oldmod,nlines)))
      no_error = 0;

  /*
   * Read in the model pointing center in (RA,Dec) format
   */

  if(no_error) {
    if(mod_center(mfp,&oldcent))
      no_error = 0;
    else
      printf("Old model center: %02d %02d %08.5f %+03d %02d %08.5f\n\n",
	     oldcent.hr,oldcent.min,oldcent.sec,oldcent.deg,oldcent.amin,
	     oldcent.asec);
  }

  /*
   * Convert the model (x,y) positions to (RA,Dec) positions
   */

  if(no_error) {
    printf("Converting (x,y) --> RA and Dec.... ");
    if(!(modpos = dspos2spos(oldcent,oldxy,nlines)))
      no_error = 0;
    else
      printf("Finished\n\n");
  }

  /*
   * Get new model center
   */

  if(no_error)
    if(get_newcent(&newcent))
      no_error = 0;

  /*
   * Convert model (RA,Dec) positions plus new model center into (x,y) offsets
   */

  if(no_error) {
    printf("Converting RA and Dec --> (x,y).... ");
    if(!(newxy = dspos2xy(newcent,modpos,nlines)))
      no_error = 0;
    else
      printf("Finished\n\n");
  }

  /*
   * Convert (x,y) offsets into (r,theta) positions and put information,
   *  along with component fluxes and sizes, into a new Compinfo structure
   */

  if(no_error)
    if(!(newmod = xy2cinfo(oldmod,newxy,nlines)))
      no_error = 0;

  /*
   * Get the output filename
   */

  if(no_error) {
    if(strcmp(ext,"edt") == 0)
      sprintf(outname,"%s_shift.mod",root);
    else
      sprintf(outname,"%s_%s_shift.mod",root,ext);
    printf("Using output file name %s\n\n",outname);
  }

  /*
   * Print results to output file
   */

  if(no_error)
    if(print_modfile(newmod,nlines,outname,newcent))
      no_error = 1;

  /*
   * Clean up
   */

  if(mfp)
    fclose(mfp);

  oldmod = del_cinfo(oldmod);
  newmod = del_cinfo(newmod);
  oldxy = del_pos(oldxy);
  newxy = del_pos(newxy);
  modpos = del_skypos(modpos);

  if(no_error)
    return 0;
  else {
    fprintf(stderr,"\nERROR.  Exiting program\n\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function cinfo2xy
 *
 * Takes a Compinfo array and converts the (r,theta) positions into
 *  (x,y) positions.  The (x,y) positions are put into a Pos array, which
 *  is returned.
 *
 * Inputs: Compinfo *cinfo     model component info
 *         int nlines          number of components
 *
 * Output: Pos *pos            array of (x,y) positions. NULL on error.
 *
 */

Pos *cinfo2xy(Compinfo *cinfo, int nlines)
{
  int i;           /* Looping variable */
  Pos *pos=NULL;   /* Position array */
  Pos *pptr;       /* Pointer to navigate pos */
  Compinfo *cptr;  /* Pointer to navigate cinfo */

  /*
   * Allocate memory for Pos array
   */

  if(!(pos = new_pos(nlines,1))) {
    fprintf(stderr,"ERROR: cinfo2xy\n");
    return NULL;
  }

  /*
   * Do the conversion.
   * NB: We have to divide the model radius by 1000 since difmap records
   *      radius in mas.
   */

  printf("Converting model (r,theta) --> (x,y)... ");
  for(i=0,pptr=pos,cptr=cinfo; i<nlines; i++,pptr++,cptr++)
    rth2xy(cptr->radius/1000.0,cptr->theta,pptr,1);
  printf("Finished\n\n");

  return pos;
}

/*.......................................................................
 *
 * Function xy2cinfo
 *
 * Takes a the old Compinfo array and a new set of (x,y) component offsets
 *  and converts the offsets into (r,theta) form.  These new (r,theta)
 *  positions are combined with the component fluxes and sizes from the
 *  old Compinfo array and a new Compinfo array is made.
 *
 * Inputs: Compinfo *old       old model component info
 *         Pos *newpos         new (x,y) offsets
 *         int nlines          number of components
 *
 * Output: Compinfo *new       New Compinfo structure.  NULL on error.
 *
 */

Compinfo *xy2cinfo(Compinfo *old, Pos *newpos, int nlines)
{
  int i;               /* Looping variable */
  Pos *pptr;           /* Pointer to navigate newpos */
  Compinfo *new=NULL;  /* New Compinfo array */
  Compinfo *optr;      /* Pointer to navigate old Compinfo array */
  Compinfo *nptr;      /* Pointer to navigate new Compinfo array */

  /*
   * Allocate memory for new Compinfo array
   */

  if(!(new = new_cinfo(nlines))) {
    fprintf(stderr,"ERROR: xy2cinfo\n");
    return NULL;
  }

  /*
   * Do the conversion.
   * NB: Multiply the radius coming out of xy2rth by 1000 to convert the
   *      units to mas.
   */

  printf("Creating new model component array... ");
  for(i=0,pptr=newpos,optr=old,nptr=new; i<nlines; i++,pptr++,optr++,nptr++) {
    *nptr = *optr;
    xy2rth(*pptr,&(nptr->radius),&(nptr->theta),1);
    nptr->radius *= 1000.0;
  }
  printf("Finished\n\n");

  return new;
}

/*.......................................................................
 *
 * Function get_newcent
 *
 * Gets RA and Dec of new pointing center
 *
 * Inputs: Skypos *newcent     new pointing center (set by this function)
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int get_newcent(Skypos *newcent)
{
  char line[MAXC];  /* String used to read input */
  Skypos tmp;       /* Temporary container for new center position */

  /*
   * Initialize to crazy values
   */

  tmp.hr = 50;
  tmp.deg = 500;

  while(tmp.hr < 0 || tmp.hr > 23) {
    printf("Enter the new RA in the form hh mm ss.sssss:  ");
    gets(line);
    if(sscanf(line,"%d %d %lf",&tmp.hr,&tmp.min,&tmp.sec) != 3) {
      fprintf(stderr,"ERROR: incorrect format\n");
      tmp.hr = 50;
    }
    if(tmp.hr < 0 || tmp.hr > 23 || tmp.min < 0 || tmp.min > 59 ||
       tmp.sec < 0.0 || tmp.sec >= 60.0)  {
      fprintf(stderr,"ERROR: incorrect format\n");
      tmp.hr = 50;
    }
  }


  while(tmp.deg < -90 || tmp.deg > 90) {
    printf("Enter the new Dec in the form dd mm ss.sssss:  ");
    gets(line);
    if(sscanf(line,"%d %d %lf",&tmp.deg,&tmp.amin,&tmp.asec) != 3) {
      fprintf(stderr,"ERROR: incorrect format\n");
      tmp.deg = 500;
    }
    if(tmp.deg < -90 || tmp.deg > 90 || tmp.amin < 0 || tmp.amin > 59 ||
       tmp.asec < 0.0 || tmp.asec >= 60.0)  {
      fprintf(stderr,"ERROR: incorrect format\n");
      tmp.deg = 500;
    }
  }

  *newcent = tmp;

  return 0;
}

/*.......................................................................
 *
 * Function print_modfile
 *
 * Given a Compinfo array and an output filename, prints the model
 *  component info to an output file.
 *
 * Inputs: Compinfo *cinfo     model component info
 *         int nlines          number of components
 *         char *outname       output filename
 *         Skypos newcent      phase center of new model
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int print_modfile(Compinfo *cinfo, int nlines, char *outname, Skypos newcent)
{
  int i;           /* Looping variable */
  Compinfo *cptr;  /* Pointer for navigating cinfo */
  FILE *ofp=NULL;  /* Output file pointer */

  /*
   * Open output file for writing
   */

  if((ofp = fopen(outname,"w")) == NULL) {
    fprintf(stderr,"ERROR: print_modfile.  Cannot open %s for output\n",
	    outname);
    return 1;
  }

  /*
   * Print header info to output file
   */

  fprintf(ofp,"! Center RA: %02d %02d %08.5f, Dec: %+03d %02d %08.5f (2000.0)\n",
	 newcent.hr,newcent.min,newcent.sec,newcent.deg,newcent.amin,
	 newcent.asec);
  fprintf(ofp,"! Established model.\n");
  fprintf(ofp,"! Flux (Jy) Radius (mas)  Theta (deg)  Major (mas)  Axial ratio");
  fprintf(ofp,"   Phi (deg) T\n");

  /*
   * Print scaled model to output file
   */

  for(i=0,cptr=cinfo; i<nlines; i++,cptr++) {
    fprintf(ofp,"%#12.6g  %11.2f %8.3f",cptr->flux, 
	      cptr->radius,cptr->theta);
    if(cptr->type == 0)
      fprintf(ofp,"\n");
    else
      fprintf(ofp," %8.4f %8.6f %9.4f %d\n",
	      cptr->maj,cptr->axrat,cptr->phi,cptr->type);
  }

  if(ofp)
    fclose(ofp);

  return 0;
}
