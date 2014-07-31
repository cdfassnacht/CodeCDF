/* distcalc.c
 *
 * Usage: distcalc -i
 * 	  distcalc -c input_file [output_file]
 * 	  distcalc -p input_file [output_file]
 * Note that the output filename is optional for both the -c and -p options.
 * 
 * OPTIONS:
 *  -i   Interactive mode.
 *  -c   Compute a set of offsets from a given CENTRAL position in
 * 	 batch mode.  For this option, the first line of the input
 * 	 file contains the position of the central source.  All
 * 	 other lines in the input file contain positions of the
 * 	 sources for which offsets from the central source will.
 * 	 be computed.
 * 	 Input file format: label ra_hr ra_min ra_sec dec_deg dec_amin dec_asec
 *  -p   Compute offsets between PAIRS of positions in batch mode.
 * 	 For this file, each line of the input file contains the
 * 	 two positions to be compared.
 * 	 Input file format: label ra_1 dec_1 ra_2 dec_2
 * 	  where ra_i has the format: ra_hr ra_min ra_sec
 * 	  and dec_i has the format:  dec_deg dec_amin dec_asec.
 *
 * Description:  Takes positions entered as RA and Dec and finds the distance 
 *                between them in arcsec.
 *
 * 18Jul95 CDF
 * v17Nov95 CDF, Changed precision of output.
 * v10Jun96 CDF, Added possibility of file input
 * v07Aug97 CDF, Massive re-write to include the much more accurate
 *                offset calculations contained in coords.c.
 *               Added documentation.
 *               Added error checking.
 *               Better I/O.
 * v21Nov97 CDF, Improved output format.
 * v09Apr99 CDF, Another much needed improvement of the output format.
 * v03May99 CDF, Changed command line format.
 *               Fixed errors in interactive input.
 * v16Oct99 CDF, Added a "pairwise offset" option for batch mode in
 *                addition to the "offset from center" option that
 *                was previously the only possibility.
 *               Changed command-line invocation slightly.
 *               Changed interactive input format for more convenience.
 * v12Jan00 CDF, Fixed memory leakage problem in the all three distance
 *                calculation functions.
 * v24Jan00 CDF, Changed output format to include source name for -c
 *                option.
 * v26Jan00 CDF, Fixed bug in output filename acquisition.
 * v04Jan01 CDF, Fixed bug in pair-wise calculation output labeling.
 * v29Jul03 CDF, Moved the print_offsets function to dataio.c
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "structdef.h"
#include "coords.h"
#include "dataio.h"

/*.......................................................................
 *
 * Enumeration
 *
 */

enum {
  CENT,
  PAIR
};



/*.......................................................................
 *
 * Function declarations
 *
 */

int calc_interactive();
int calc_cent(char *inname, char *outname);
int calc_pair(char *inname, char *outname);
void distcalc_help();

/*.......................................................................
 *
 * Main Program
 *
 */

int main(int argc, char *argv[])
{
  int no_error = 1;     /* Flag set to 0 on error */

  /*
   * Check input line format
   */

  switch(argc) {
  case 1:
    distcalc_help();
    return 1;
    break;
  case 2:
    if(strcmp(argv[1],"-i") == 0) {
      /* Interactive input */
      if(calc_interactive())
	no_error = 0;
    }
    else {
      distcalc_help();
    }
    break;
  case 3:


    /* File-based input with NO designated output filename */
    if(strcmp(argv[1],"-c") == 0) {
      if(calc_cent(argv[2],NULL))
	no_error = 0;
    }
    else if(strcmp(argv[1],"-p") == 0) {
      if(calc_pair(argv[2],NULL))
	no_error = 0;
    }
    else if(strcmp(argv[1],"-i") == 0) {
      if(calc_interactive())
	no_error = 0;
    }
    else {
      distcalc_help();
    }
    break;
  case 4:
    /* File-based input with a designated output filename */
    if(strcmp(argv[1],"-c") == 0) {
      if(calc_cent(argv[2],argv[3]))
	no_error = 0;
    }
    else if(strcmp(argv[1],"-p") == 0) {
      if(calc_pair(argv[2],argv[3]))
	no_error = 0;
    }
    else if(strcmp(argv[1],"-i") == 0) {
      if(calc_interactive())
	no_error = 0;
    }
    else {
      distcalc_help();
    }
    break;
  default:
    fprintf(stderr,"\n***WARNING: Too many arguments.***\n\n");
    distcalc_help();
    return 1;
    break;
  }

  /*
   * Clean up and exit
   */


  if(no_error) {
    printf("\nCompleted program distcalc.\n\n");
    return 0;
  }
  else {
    fprintf(stderr,"ERROR.  Exiting program\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function calc_interactive
 *
 * Takes positions entered by the user and then calls print_offsets.
 *
 * Inputs: none
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 */

int calc_interactive()
{
  int no_error = 1;     /* Flag set to 0 on error */
  char line[MAXC];      /* General string for reading input */
  Pos *offsets=NULL;    /* Output list of position offsets */
  Skypos cent;          /* First sky position entered */
  Skypos *skypos=NULL;  /* List of all other sky positions */

  /*
   * Allocate memory for second position
   */

  if(!(skypos = new_skypos(1))) {
    fprintf(stderr,"ERROR: Exiting program.\n");
    return 1;
  }

  /*
   * Fill in central RA and Dec
   */

  printf("Enter the RA and Dec of the first object ");
  printf("(hh mm ss.sss dd mm ss.ss)\n");
  printf("Position 1: ");
  fgets(line,MAXC,stdin);
  while(sscanf(line,"%d %d %lf %d %d %lf",&cent.hr,&cent.min,&cent.sec,
	       &cent.deg,&cent.amin,&cent.asec) != 6) {
    fprintf(stderr,"ERROR: bad input\n");
    fprintf(stderr,"Enter coordinates again:  ");
    fgets(line,MAXC,stdin);
  }

  /*
   * Fill in second position RA and Dec
   */

  printf("Enter the RA and Dec of the second object ");
  printf("(hh mm ss.sss dd mm ss.ss)\n");
  printf("Position 2: ");
  fgets(line,MAXC,stdin);
  while(sscanf(line,"%d %d %lf %d %d %lf",&skypos->hr,&skypos->min,
	       &skypos->sec,&skypos->deg,&skypos->amin,&skypos->asec) != 6) {
    fprintf(stderr,"ERROR: bad input\n");
    fprintf(stderr,"Enter coordinates again:  ");
    fgets(line,MAXC,stdin);
  }

  /*
   * Calculate the separation and put into a Pos array
   */

  if(!(offsets = dspos2xy(cent,skypos,1)))
    no_error = 0;

  /*
   * Clear label string
   */

  sprintf(cent.label,"");

  /*
   * Print out results
   */

  if(print_offsets(cent,skypos,offsets,1,NULL))
    no_error = 0;

  /*
   * Clean up and exit
   */

  offsets = del_pos(offsets);
  skypos = del_skypos(skypos);

  if(no_error)
    return 0;
  else {
    fprintf(stderr,"ERROR: calc_interactive.\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function calc_cent
 *
 * Reads in lines from a file to calculate a set of offsets from a central
 *  position.
 *
 * Inputs: char *inname        input filename
 *         char *outname       output filename (NULL if none designated on
 *                              command line)
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 */

int calc_cent(char *inname, char *outname)
{
  int no_error = 1;     /* Flag set to 0 on error */
  int nlines;           /* Number of lines in the input file */
  char newoutn[MAXC];   /* New output name on error */
  char line[MAXC];      /* General string for reading input */
  Pos *offsets=NULL;    /* Output list of position offsets */
  Skypos cent;          /* First sky position entered */
  Skypos *skypos=NULL;  /* List of all other sky positions */
  Skypos *sptr;         /* Pointer used to navigate skypos */
  FILE *ifp=NULL;       /* Input file pointer */
  FILE *ofp=NULL;       /* Output file pointer */
 
  /*
   * Open the input file
   */

  if(!(ifp = open_readfile(inname))) {
    fprintf(stderr,"ERROR: calc_cent.\n");
    return 1;
  }

  /*
   * Open the output file if one is desired
   */

  if(outname == NULL) {
    printf("\nEnter name of output file if one is desired or\n");
    printf("just hit return for no output file: [no output file] ");
    fgets(line,MAXC,stdin);
    if(line[0] != '\n') {
      while(sscanf(line,"%s",newoutn) != 1) {
	fprintf(stderr,"ERROR: Bad input.  Enter filename again:  ");
	fgets(line,MAXC,stdin);
      }
      if(!(ofp = open_writefile(newoutn))) {
	fprintf(stderr,"ERROR: calc_cent\n");
	return 1;
      }
    }
  }
  else {
    if(!(ofp = open_writefile(outname))) {
      fprintf(stderr,"ERROR: calc_cent\n");
      return 1;
    }
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
      printf("\nRead %d lines from %s\n",nlines,inname);
    }
  }

  /*
   * Get the central position
   */

  while(fgets(line, MAXC, ifp) != NULL && no_error) {
    if(line[0] != '!') {
      if(sscanf(line,"%s %d %d %lf %d %d %lf",cent.label,&cent.hr,&cent.min,
		&cent.sec,&cent.deg,&cent.amin,&cent.asec) != 7) {
	fprintf(stderr,"ERROR.  Bad input file format.\n");
	fprintf(stderr,"This program requires the following format:\n");
	fprintf(stderr,"  label rahr ramin rasec decdeg decmin decsec.\n");
	no_error = 0;
      }
      else {
	break;
      }
    }
  }

  /*
   * Allocate container for all other positions
   */

  if(no_error)
    if(!(skypos = new_skypos(nlines-1)))
      no_error = 0;

  /*
   * Read in the other positions
   */

  if(no_error) {
    sptr = skypos;
    while(fgets(line, MAXC, ifp) != NULL && no_error) {
      if(line[0] != '!') {
	if(sscanf(line,"%s %d %d %lf %d %d %lf",sptr->label,
		  &sptr->hr,&sptr->min,&sptr->sec,
		  &sptr->deg,&sptr->amin,&sptr->asec) != 7) {
	  fprintf(stderr,"ERROR.  Bad input file format.\n");
	  fprintf(stderr,"This program requires the following format:\n");
	  fprintf(stderr,"  label rahr ramin rasec decdeg decmin decsec.\n");
	  no_error = 0;
	}
	else
	  sptr++;
      }
    }
  }

  /*
   * Calculate the positional offsets
   */

  if(no_error)
    if(!(offsets = dspos2xy(cent,skypos,nlines-1)))
      no_error = 0;

  /*
   * Print out the offsets
   */

  if(no_error)
    if(print_offsets(cent,skypos,offsets,nlines-1,ofp))
      no_error = 0;

  /*
   * Clean up and exit
   */

  if(ifp)
    fclose(ifp);
  if(ofp)
    fclose(ofp);
  offsets = del_pos(offsets);
  skypos = del_skypos(skypos);

  if(no_error)
    return 0;
  else {
    fprintf(stderr,"ERROR: calc_cent\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function calc_pair
 *
 * Reads in lines from a file to calculate a set of offsets in a pairwise
 *  fashion.  That is, each line of the input file contains a pair of
 *  coordinates.  This function calculates the offset between the pair
 *  of coordinates for each line in the input file separately.
 *
 * Inputs: char *inname        input filename
 *         char *outname       output filename (NULL if none designated on
 *                              command line)
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 */

int calc_pair(char *inname, char *outname)
{
  int no_error = 1;     /* Flag set to 0 on error */
  int nlines;           /* Number of lines in the input file */
  char newoutn[MAXC];   /* New output name on error */
  char line[MAXC];      /* General string for reading input */
  Pos *offsets=NULL;    /* Output list of position offsets */
  Skypos cent;          /* First sky position entered */
  Skypos *skypos=NULL;  /* List of all other sky positions */
  FILE *ifp=NULL;       /* Input file pointer */
  FILE *ofp=NULL;       /* Output file pointer */
 
  /*
   * Open the input file
   */

  if(!(ifp = open_readfile(inname))) {
    fprintf(stderr,"ERROR: calc_cent.\n");
    return 1;
  }

  /*
   * Open the output file if one is desired
   */

  if(outname == NULL) {
    printf("\nEnter name of output file if one is desired or\n");
    printf("just hit return for no output file: [no output file] ");
    fgets(line,MAXC,stdin);
    if(line[0] != '\n') {
      while(sscanf(line,"%s",newoutn) != 1) {
	fprintf(stderr,"ERROR: Bad input.  Enter filename again:  ");
	fgets(line,MAXC,stdin);
      }
      if(!(ofp = open_writefile(newoutn))) {
	fprintf(stderr,"ERROR: calc_pair\n");
	return 1;
      }
    }
  }
  else {
    if(!(ofp = open_writefile(outname))) {
      fprintf(stderr,"ERROR: calc_pair\n");
      return 1;
    }
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
      printf("\nRead %d lines from %s\n",nlines,inname);
    }
  }

  /*
   * Allocate memory for the one-coordinate container, skypos.
   */

  if(no_error)
    if(!(skypos = new_skypos(1)))
      no_error = 0;

  /*
   * Now cycle through the input file, reading in a pair of coordinates
   *  for each line and calculating the offset between them.
   */

  while(fgets(line, MAXC, ifp) != NULL && no_error) {
    if(line[0] != '!') {
      if(sscanf(line,"%s %d %d %lf %d %d %lf %d %d %lf %d %d %lf",
		cent.label,&cent.hr,&cent.min,&cent.sec,
		&cent.deg,&cent.amin,&cent.asec,
		&skypos->hr,&skypos->min,&skypos->sec,
		&skypos->deg,&skypos->amin,&skypos->asec) != 13) {
	fprintf(stderr,"ERROR.  Bad input file format.\n");
	fprintf(stderr,"This program requires the following format:\n");
	fprintf(stderr," label rahr1 ramin1 rasec1 decdeg1 decmin1 decsec1");
	fprintf(stderr," rahr2 ramin2 rasec2 decdeg2 decmin2 decsec2.\n");
	no_error = 0;
      }

      /*
       * Calculate the positional offsets
       */

      if(no_error) {
	strcpy(skypos->label,cent.label);
	if(!(offsets = dspos2xy(cent,skypos,1)))
	  no_error = 0;
      }

      /*
       * Print out the offsets
       */

      if(no_error)
	if(print_offsets(cent,skypos,offsets,1,ofp))
	  no_error = 0;
    }
  }

  /*
   * Clean up and exit
   */

  if(ifp)
    fclose(ifp);
  if(ofp)
    fclose(ofp);
  offsets = del_pos(offsets);
  skypos = del_skypos(skypos);

  if(no_error)
    return 0;
  else {
    fprintf(stderr,"ERROR: calc_pair\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function distcalc_help
 *
 * Prints out helpful information when the program is called with no 
 *  arguments.
 *
 * Inputs: none
 *
 * Output: none
 *
 */

void distcalc_help()
{
  printf("\n*******************************************************\n\n");
  printf("\nProgram: distcalc -- ");
  printf("Calculates offsets between pairs of coordinates.\n\n");
  printf("Usage: distcalc -i\n");
  printf("       distcalc -c input_file [output_file]\n");
  printf("       distcalc -p input_file [output_file]\n\n");
  printf("Note that the output filename is optional for both the -c and -p ");
  printf("options.\n\n");
  printf("OPTIONS:\n");
  printf(" -i   Interactive mode.\n\n");
  printf(" -c   Compute a set of offsets from a given CENTRAL position in\n");
  printf("      batch mode.  For this option, the first line of the input\n");
  printf("      file contains the position of the central source.  All\n");
  printf("      other lines in the input file contain positions of the\n");
  printf("      sources for which offsets from the central source will.\n");
  printf("      be computed.\n");
  printf("      Input file format: label ra_hr ra_min ra_sec dec_deg dec_amin ");
  printf("dec_asec\n\n");
  printf(" -p   Compute offsets between PAIRS of positions in batch mode.\n");
  printf("      For this file, each line of the input file contains the\n");
  printf("      two positions to be compared.\n");
  printf("      Input file format: label ra_1 dec_1 ra_2 dec_2\n");
  printf("       where ra_i has the format: ra_hr ra_min ra_sec\n");
  printf("       and dec_i has the format:  dec_deg dec_amin dec_asec.\n\n");
  printf("*******************************************************\n\n");
}
