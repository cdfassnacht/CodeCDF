/* add_offsets.c
 *
 * Usage: add_offsets input_file [output_file]

 * Note that the output filename is optional.
 * 
 * Description:  Calculates sky positions in RA and Dec from a given
 *                central position and offsets in arcsec
 *
 * 27Jun00 CDF,  A modification of distcalc.c
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "structdef.h"
#include "coords.h"
#include "dataio.h"

/*.......................................................................
 *
 * Function declarations
 *
 */

Pos *read_offset_file(char *inname, Skypos *cent, int *noffsets);
int print_spos(Skypos cent, Skypos *skypos, Pos *offsets, int noffsets,
	       char *outname);

/*.......................................................................
 *
 * Main Program
 *
 */

int main(int argc, char *argv[])
{
  int no_error = 1;     /* Flag set to 0 on error */
  int noffsets;         /* Number of offsets */
  char inname[MAXC];    /* Input filename */
  char outname[MAXC];   /* Input filename */
  Pos *offsets=NULL;    /* Array of offsets from the central position */
  Skypos cent;          /* Central position */
  Skypos *skypos=NULL;  /* Array of calculated sky positions */

  /*
   * Check input line format
   */

  if(argc<2 || argc>3) {
    fprintf(stderr,"\nUsage: add_offsets input_file [output_file]\n\n");
    fprintf(stderr,"Format of input_file:\n");
    fprintf(stderr,"  First line (central coordinate): ");
    fprintf(stderr," label rahr ramin rasec decdeg decmin decsec\n");
    fprintf(stderr,"  All following lines (offsets): label dra ddec\n");
    fprintf(stderr," NB: dra and ddec are in arcsec\n\n");
    return 1;
  }

  /*
   * Read input file
   */

  strcpy(inname,argv[1]);
  if((offsets = read_offset_file(inname,&cent,&noffsets)) == NULL ||
     noffsets == 0) {
    fprintf(stderr,"ERROR.  Problem with input file.\n");
    no_error = 0;
  }

  /*
   * Calculate the sky positions from the offsets.
   */

  if(no_error)
    if((skypos = dspos2spos(cent,offsets,noffsets)) == NULL)
      no_error = 0;

  /*
   * Print out the calculated sky positions.
   */

  if(no_error) {
    if(argc == 3) {
      strcpy(outname,argv[2]);
      if(print_spos(cent,skypos,offsets,noffsets,outname))
	no_error = 0;
    }
    else
      if(print_spos(cent,skypos,offsets,noffsets,NULL))
	no_error = 0;
  }

  /*
   * Clean up and exit
   */

  offsets = del_pos(offsets);
  skypos = del_skypos(skypos);

  if(no_error) {
    printf("\nCompleted program add_offsets.\n\n");
    return 0;
  }
  else {
    fprintf(stderr,"ERROR.  Exiting program\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function read_offset_file
 *
 * Reads in central position and offsets from an offset file.  The file
 *  should have as its first line the central position in the form
 *  "label hh mm ss.sss  dd mm ss.ss".  All of the following lines will 
 *  contain offsets in the form "label dra ddec" in arcsec.
 *
 * Inputs: char *inname        name of input file
 *         Skypos *cent        central position (set by this function)
 *         int *noffsets       number of offsets in data file (set by
 *                              this function)
 *
 * Output: Pos *offsets        array containing offsets
 *
 */

Pos *read_offset_file(char *inname, Skypos *cent, int *noffsets)
{
  int no_error=1;       /* Flag set to 0 on error */
  int nlines=0;         /* Number of data lines in input file */
  int count=0;          /* Number of lines actually read in */
  char line[MAXC];      /* General string variable for reading inputs */
  Pos *offsets=NULL;    /* Array for offsets */
  Pos *pptr;            /* Pointer to step through offsets array */
  FILE *ifp=NULL;       /* Input file pointer */

  /*
   * Open the input file
   */

  if(!(ifp = open_readfile(inname))) {
    fprintf(stderr,"ERROR: read_offset_file.\n");
    return NULL;
  }

  /*
   * Count the number of lines in the input file
   */

  if(no_error) {
    if((nlines = n_lines(ifp,'!')) == 0) {
      fprintf(stderr,"ERROR.  %s has no valid data lines.\n",inname);
      no_error = 0;
    }
    else {
      rewind(ifp);
      printf("%d data lines in %s\n\n",nlines,inname);
    }
  }

  /*
   * Allocate memory for the array of positional offsets
   */

  if(no_error) {
    if(!(offsets = new_pos(nlines-1,1)))
      no_error = 0;
    else
      pptr = offsets;
  }


  /*
   * Now step through the input file, reading in data lines.
   */

  while(fgets(line, MAXC, ifp) != NULL && no_error) {
    if(line[0] != '!') {
      count++;

      /*
       * If first data line, read in central position.
       */

      if(count == 1) {
	if(sscanf(line,"%s %d %d %lf %d %d %lf",
		  cent->label,&cent->hr,&cent->min,&cent->sec,
		  &cent->deg,&cent->amin,&cent->asec) != 7) {
	  fprintf(stderr,"ERROR: read_offset_file.\n");
	  fprintf(stderr,"   First data line of input file must be of ");
	  fprintf(stderr,"the form\n");
	  fprintf(stderr,"   label rahr ramin rasec decdeg decmin decsec\n");
	  no_error = 0;
	}
      }

      /*
       * For all other data lines read in positional offsets.
       */

      else if(line[0] != '\n') {
	if(sscanf(line,"%s %lf %lf",pptr->label,&pptr->x,&pptr->y) != 3) {
	  fprintf(stderr,"ERROR: read_offset_file.\n");
	  fprintf(stderr,"   Positional offset data lines must be of ");
	  fprintf(stderr,"the form\n");
	  fprintf(stderr,"     label dra ddec,\n");
	  fprintf(stderr,"   where dra and ddec are in arcsec.\n");
	  no_error = 0;
	}
	else
	  pptr++;
      }
    }
  }

  /*
   * Clean up and exit
   */

  if(no_error) {
    *noffsets = nlines-1;
    return offsets;
  }
  else {
    fprintf(stderr,"ERROR: read_offset_file\n");
    return(del_pos(offsets));
  }
}

/*.......................................................................
 *
 * Function print_spos
 *
 * Prints out the sky positions calculated from the given central position
 *  and the list of offsets.
 * NB:  This produces output of the same format as distcalc.c.
 *
 * Inputs: Skypos cent         central position (RA, Dec)
 *         Skypos *skypos      calculated positions (RA, Dec)
 *         Pos *offsets        (x,y) offsets
 *         int noffsets        number of offsets
 *         char *outname       output file name (NULL if no output file)
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int print_spos(Skypos cent, Skypos *skypos, Pos *offsets, int noffsets,
	       char *outname)
{
  int i;          /* Looping variable */
  char ew[5];     /* Direction of delta_RA */
  char ns[5];     /* Direction of delta_Dec */
  Skypos *sptr;   /* Pointer used to navigate skypos */
  Pos *optr;      /* Pointer used to navigate offsets */
  FILE *ofp=NULL; /* Output file pointer */

  /*
   * Open the output file if one is desired
   */

  if(outname)
    if(!(ofp = open_writefile(outname))) {
      fprintf(stderr,"\nUnable to open output file.  Sending output to ");
      fprintf(stderr,"STDOUT.\n");
    }

  /*
   * Print out the central position
   */

  if(ofp) {
    fprintf(ofp,"#\n");
    fprintf(ofp,"#\n");
    fprintf(ofp,"#                CENTRAL SOURCE:  %s\n",cent.label);
    fprintf(ofp,"#\n");
    fprintf(ofp,"#Central source:  %02d %02d %07.4f  %+03d %02d %07.4f\n",
	    cent.hr,cent.min,cent.sec,cent.deg,cent.amin,cent.asec);
    fprintf(ofp,"#\n");
    fprintf(ofp,"#                                         ");
    fprintf(ofp,"Shift FROM central source\n");
    fprintf(ofp,"#                                         ");
    fprintf(ofp,"TO listed source (arcsec)\n");
    fprintf(ofp,"#Source          RA              Dec      ");
    fprintf(ofp,"    d_RA          d_Dec   Total Sep.\n");
    fprintf(ofp,"#----------  -------------  -------------  ");
    fprintf(ofp,"-----------  -----------  ---------\n");
  }
  else {
    printf("#\n");
    printf("#\n");
    if(strcmp(cent.label,"") != 0) {
      printf("#                CENTRAL SOURCE:  %s\n",cent.label);
    }
    printf("#\n");
    printf("#\n");
    printf("#Central source:  %02d %02d %07.4f  %+03d %02d %07.4f\n",
	   cent.hr,cent.min,cent.sec,cent.deg,cent.amin,cent.asec);
    printf("#\n");
    printf("#                                         ");
    printf("Shift FROM central source\n");
    printf("#                                         ");
    printf("TO listed source (arcsec)\n");
    printf("#Source          RA              Dec      ");
    printf("    d_RA          d_Dec   Total Sep.\n");
    printf("#----------  -------------  -------------  ");
    printf("-----------  -----------  ---------\n");
  }

  for(i=0,sptr=skypos,optr=offsets; i<noffsets; i++,sptr++,optr++) {

    /*
     * Determine the direction of the shift
     */

    if(optr->x >= 0.0)
      sprintf(ew,"E");
    else
      sprintf(ew,"W");
    if(optr->y >= 0.0)
      sprintf(ns,"N");
    else
      sprintf(ns,"S");

    /*
     * Print output to file or to STDOUT
     */

    if(ofp) {
      fprintf(ofp,"%11s  %02d %02d %07.4f  %+03d %02d %06.3f ",
	      optr->label,sptr->hr,sptr->min,sptr->sec,
	      sptr->deg,sptr->amin,sptr->asec);
      fprintf(ofp,"%10.4f %1s %10.4f %1s %10.4f\n",
	      fabs(optr->x),ew,fabs(optr->y),ns,
	      sqrt(optr->x * optr->x + optr->y * optr->y));
    }
    else {
      printf("%11s  %02d %02d %07.4f  %+03d %02d %06.3f ",
	      optr->label,sptr->hr,sptr->min,sptr->sec,
	      sptr->deg,sptr->amin,sptr->asec);
      printf("%10.4f %1s %10.4f %1s %10.4f\n",
	      fabs(optr->x),ew,fabs(optr->y),ns,
	      sqrt(optr->x * optr->x + optr->y * optr->y));
    }
  }

  if(ofp) {
    fprintf(ofp,"#\n");
    fprintf(ofp,
	    "#-----------------------------------------------------------");
    fprintf(ofp,"-----------------\n");
  }
  else {
    printf("#\n");
    printf("#-----------------------------------------------------------");
    printf("-----------------\n");
  }

  /*
   * Clean up and exit.
   */

  if(ofp)
    fclose(ofp);
    
  return 0;
}
 
