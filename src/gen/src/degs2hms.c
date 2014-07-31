/* degs2hms.c
 *
 * Usage: degs2hms -i
 *        degs2hms -b [input_file] {[output_file]}
 *
 * Description: Takes a set of coordinates in decimal degree format and
 *               converts them to RA,Dec format (hms).  
 *
 * 17Nov2004, CDF  A inversion of hms2degs.c
 *
 */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "structdef.h"
#include "coords.h"
#include "dataio.h"

/*.......................................................................
 *
 * Function declarations
 *
 */

int degs2hms_interactive();
int degs2hms_batch(char *inname, char *outname);
void degs2hms_help();

/*.......................................................................
 *
 * Main program
 *
 */


int main(int argc, char *argv[])
{
  int no_error=1;  /* Flag set to 0 on error */

  /*
   * Check input line format
   */

  switch(argc) {
  case 1:
    degs2hms_help();
    return 1;
    break;
  case 2:
    if(strcmp(argv[1],"-i") == 0) {
      /* Interactive input */
      if(degs2hms_interactive())
	no_error = 0;
    }
    else {
      degs2hms_help();
    }
    break;
  case 3:
    /* File-based input with NO designated output filename */
    if(strcmp(argv[1],"-b") == 0) {
      if(degs2hms_batch(argv[2],NULL))
	no_error = 0;
    }
    else if(strcmp(argv[1],"-i") == 0) {
      if(degs2hms_interactive())
	no_error = 0;
    }
    else {
      degs2hms_help();
    }
    break;
  case 4:
    /* File-based input with a designated output filename */
    if(strcmp(argv[1],"-b") == 0) {
      if(degs2hms_batch(argv[2],argv[3]))
	no_error = 0;
    }
    else if(strcmp(argv[1],"-i") == 0) {
      if(degs2hms_interactive())
	no_error = 0;
    }
    else {
      degs2hms_help();
    }
    break;
  default:
    fprintf(stderr,"\n***WARNING: Too many arguments.***\n\n");
    degs2hms_help();
    return 1;
    break;
  }

  return 0;
}

/*.......................................................................
 *
 * Function degs2hms_interactive
 *
 * Takes positions entered by the user and prints output.
 *
 * Inputs: none
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 */

int degs2hms_interactive()
{
  int go_on=1;      /* Flag set to 0 when user wants to quit */
  double alphadeg;  /* RA in decimal format */
  double deltadeg;  /* Dec in decimal format */
  char line[MAXC];  /* General string for reading input */
  Skypos skypos;    /* Sky position in hh mm ss.sss dd mm ss.sss format */

  /*
   * Loop until the user wants to quit
   */

  while(go_on) {

    /*
     * Fill in RA and Dec
     */

    printf("Enter the alpha and delta to be converted.\n");
    printf(" (format is ddd.ddddd  sddd.dddd):  ");
    fgets(line,MAXC,stdin);
    while(sscanf(line,"%lf %lf",&alphadeg,&deltadeg) != 2) {
      fprintf(stderr,"ERROR: bad input\n");
      fprintf(stderr,"Enter coordinates again:  ");
      fgets(line,MAXC,stdin);
    }

    /*
     * Convert to hms format
     */

    deg2spos(alphadeg,deltadeg,&skypos);

    /*
     * Print out results
     */

    printf("\n%9.5f %9.5f ---> %02d:%02d:%07.4f %+03d:%02d:%06.3f\n",
	   alphadeg,deltadeg,
	   skypos.hr,skypos.min,skypos.sec,
	   skypos.deg,skypos.amin,skypos.asec);

    /*
     * Check if further runs are required
     */

    printf("\nDo another conversion? [y] ");
    fgets(line,MAXC,stdin);
    if(line[0] == 'n' || line[0] == 'N')
      go_on = 0;
  }

  printf("\nFinished with degs2hms.c\n\n");
  return 0;
}


/*.......................................................................
 *
 * Function degs2hms_batch
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

int degs2hms_batch(char *inname, char *outname)
{
  int no_error = 1;    /* Flag set to 0 on error */
  int nlines;          /* Number of lines in the input file */
  double alphadeg;     /* RA in decimal format */
  double deltadeg;     /* Dec in decimal format */
  char line[MAXC];     /* General string for reading input */
  char newname[MAXC];  /* New name for output file */
  Skypos skypos;       /* First sky position entered */
  FILE *ifp=NULL;      /* Optional input file */
  FILE *ofp=NULL;      /* Optional output file */

 
  /*
   * Open the input file
   */

  if(!(ifp = open_readfile(inname))) {
    fprintf(stderr,"ERROR: calc_skypos.\n");
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
      while(sscanf(line,"%s",newname) != 1) {
	fprintf(stderr,"ERROR: Bad input.  Enter filename again:  ");
	fgets(line,MAXC,stdin);
      }
      if(!(ofp = open_writefile(newname)))
	no_error = 0;
    }
  }
  else {
    if(!(ofp = open_writefile(outname)))
      no_error = 0;
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
      printf("\nFound %d lines in %s\n",nlines,inname);
      rewind(ifp);
    }
  }

  /*
   * Print out header
   */

  if(ofp) {
    fprintf(ofp,"#  Object          RA           Dec   \n");
    fprintf(ofp,"#----------   ------------ -------------\n");
  }
  else {
    printf("\n  Object        RA             Dec    \n");
    printf("----------   ------------- -------------\n");
  }

  /*
   * Now cycle through the input file, reading in a pair of coordinates
   *  for each line and calculating the offset between them.
   */

  while(fgets(line, MAXC, ifp) != NULL && no_error) {
    if(line[0] != '!') {
      if(sscanf(line,"%s %lf %lf",skypos.label,&alphadeg,&deltadeg) != 3) {
	fprintf(stderr,"ERROR.  Bad input file format.\n");
	fprintf(stderr,"This program requires the following format:\n");
	fprintf(stderr," label alpha delta.\n");
	no_error = 0;
      }

      if(no_error) {

	/*
	 * Convert to RA, Dec
	 */

	deg2spos(alphadeg,deltadeg,&skypos);

	/*
	 * Print out the offsets
	 */

	if(ofp)
	  fprintf(ofp,"%-12s %02d %02d %07.4f %+03d %02d %06.3f\n",skypos.label,
		  skypos.hr,skypos.min,skypos.sec,
		  skypos.deg,skypos.amin,skypos.asec);
	else
	  printf("%-12s %02d %02d %07.4f %+03d %02d %06.3f\n",skypos.label,
		 skypos.hr,skypos.min,skypos.sec,
		 skypos.deg,skypos.amin,skypos.asec);
      }
    }
  }

  /*
   * Clean up and exit
   */

  if(ifp)
    fclose(ifp);
  if(ofp)
    fclose(ofp);

  if(no_error) {
    printf("\nExiting program degs2hms.c\n\n");
    return 0;
  }
  else {
    fprintf(stderr,"ERROR: degs2hms_batch\n");
    return 1;
  }
}



/*.......................................................................
 *
 * Function degs2hms_help
 *
 * Prints out helpful information when the program is called with no 
 *  arguments.
 *
 * Inputs: none
 *
 * Output: none
 *
 */

void degs2hms_help()
{
  printf("\n*******************************************************\n\n");
  printf("\nProgram: degs2hms -- ");
  printf("Converts positions in decimal degree format into RA,Dec format\n\n");
  printf("Usage: degs2hms -i\n");
  printf("       degs2hms -b input_file [output_file]\n");
  printf("Note that the output filename is optional for the -b option.\n\n");
  printf("OPTIONS:\n");
  printf(" -i   Interactive mode.\n\n");
  printf(" -b   Convert positions in a batch mode from an input file.\n");
  printf("      For this option, each line of the input file contains a\n");
  printf("      position in alpha, delta (decimal degrees) format as follows:");
  printf("\n\n");
  printf("      label alpha \n\n");
  printf("*******************************************************\n\n");
}
