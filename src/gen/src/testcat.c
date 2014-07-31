/* testcat.c
 * 
 * Program to test functions before they are incorporated into other programs
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "structdef.h"
#include "dataio.h"
#include "catlib.h"

int main(int argc, char *argv[])
{
  int i;                /* Looping variable */
  int no_error=1;       /* Flag set to 0 on error */
  int keep_reading=1;   /* Flag used for reading header */
  int nheader=0;       /* Number of header lines */
  int nlines=0;
  int ncols=0;
  int col;
  int ncent;            /* Number of lines in posfile (should be 1) */
  int ncat;             /* Number of lines in the input catalog */
  int format=6;         /* Format of input/output files */
  float test1,test2,test3;
  char comment='#';     /* SExtractor comment character */
  char posfile[MAXC];   /* Filename for input position file */
  char catfile[MAXC];   /* Filename for input catalog */
  char outfile[MAXC];   /* Filename for distcalc-like output file */
  char junk[MAXC];
  char keyword[MAXC];
  char line[MAXC];
  char strformat[MAXC]; /* Format string */
  Pos *offsets=NULL;    /* Array to hold offsets */
  Skypos *skypos=NULL;  /* Array of sky positions */
  Skypos *skptr;        /* Pointer to navigate skypos */
  Secat *centpos=NULL;  /* Central position for distance calculations */
  Secat *secat=NULL;    /* Data array from catalog */
  Secat *sptr;          /* Pointer to navigate secat */
  FILE *ifp=NULL;       /* Input file pointer */
  FILE *ofp=NULL;       /* Output file pointer */

  /*
   * Get the inputs from the command line.
   */

  strcpy(catfile,argv[1]);

  /*
   * Open input file
   */

  if(!(ifp = open_readfile(catfile))) {
    fprintf(stderr,"ERROR: read_secat.\n");
    return ERROR;
  }

  /*
   * Get number of lines in input file
   */

  if((nlines = n_lines(ifp,comment)) == 0) {
    fprintf(stderr,"ERROR: read_secat.  No valid data in input file.\n");
    no_error = 0;
  }
  else {
    rewind(ifp);
  }

  /*
   * Read in header information
   */

  printf("\nChecking for header information, (expected to begin ");
  printf("on the first line)\n");
  while(no_error && keep_reading && fgets(line,MAXC,ifp) != NULL) {
    if(line[0] == comment) {
      printf("%s",line);
      nheader++;
      if(sscanf(line,"%s %d %s",junk,&col,keyword) != 3) {
	fprintf(stderr,"ERROR: read_secat.  Bad format in header for %s\n",
		catfile);
	no_error = 0;
      }
    }
    else
      keep_reading = 0;
  }
  printf("\nFound %d header lines.\n",nheader);
  ncols = n_cols(line,comment,1);
  sscanf(line,"%f",&test1);
  sscanf(line,"%f",&test2);
  printf("%f %f\n",test1,test2);
#if 0
  sprintf(strformat,"%%d %%f %%f");
  sscanf(line,strformat,&col,&test2,&test3);
  printf("%d %f %f\n",col,test2,test3);
#endif
  printf("\n");
  

  /*
   * Read the input catalog
   */

  if(no_error)
    if(!(secat = read_secat(catfile,'#',&ncat,format)))
      no_error = 0;

  /*
   * Clean up and exit
   */

  secat = del_secat(secat);
  if(ofp)
    fclose(ofp);

  if(no_error) {
    printf("\nProgram testcat finished.\n\n");
    return 0;
  }
  else {
    fprintf(stderr,"\nERROR.  Exiting testcat.\n\n");
    return 1;
  }

  return 0;

}
