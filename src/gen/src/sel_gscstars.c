/*
 * sel_gscstars.c
 *
 * Description: Takes the output of the GSC star catalog program
 *   found on the ESO GSC web server (CDF home page under Astronomy)
 *   are closest to the requested field center.  This program is
 *   useful for preparing for P60 spectroscopy observing runs because
 *   the guider camera on the P60 can only see objects down to
 *   13th or 14th mag, so blind offsetting needs to be done.
 *
 * Usage: sel_gscstars [infile] [outfile]
 *    where [infile] is a file containing a list of files to be processed
 *    and [outfile] is the name of the desired output file.
 *
 * Output: a list of 3 appropriate GSC stars per input file.
 *
 * 18Nov97 CDF
 */

#include <stdio.h>
#include <stdlib.h>
#include "structdef.h"
#include "dataio.h"

#define SKIPBEG 22    /* Number of lines to skip at beginning of files */
#define SKIPEND 23    /* Number of lines to skip at end of files */
#define FAINTMAG 100.0 /* Faintest magnitude for an acceptable GSC star */

/*.......................................................................
 *
 * Function declarations
 *
 */

int read_gsc_file(char *gname, FILE *ofp);

/*.......................................................................
 *
 * Main program
 *
 */

int main(int argc, char *argv[])
{
  int nl1;         /* Number of lines in filelist file */
  int nstars;      /* Number of good stars found in GSC file */
  char *flname;    /* Filelist filename */
  char *outname;   /* Output filename */
  char line[MAX];  /* General string for reading input */
  char gname[MAX]; /* Names of GSC star files */
  FILE *ifp;       /* Pointer for input filelist file */
  FILE *ofp;       /* Pointer for output file */

  /*
   * Check command line format
   */

  if(argc != 3) {
    fprintf(stderr,"\nUsage: sel_gscstars [infile] [outfile]\n");
    fprintf(stderr," where [infile] is a file containing a list of files ");
    fprintf(stderr,"to be processed\n");
    fprintf(stderr," and [outfile] is the name of the desired output file.");
    fprintf(stderr,"\n\n");
    return 1;
  }

  /*
   * Open filelist file
   */

  flname = argv[1];
  while((ifp = fopen(flname,"r")) == NULL) {
    fprintf(stderr,"ERROR: %s does not exist.  Enter new filename:  ",
	    flname);
    gets(flname);
  }

  /*
   * Open output file
   */

  outname = argv[2];
  while((ofp = fopen(outname,"w")) == NULL) {
    fprintf(stderr,"ERROR: %s does not exist.  Enter new filename:  ",
	    outname);
    gets(outname);
  }

  /*
   * Count number of lines in filelist file
   */

  if((nl1 = n_lines(ifp,'!')) == 0) {
    fprintf(stderr,"\nERROR: No lines in input file\n");
    return 1;
  }
  else
    rewind(ifp);

  /*
   * Open each GSC star file in turn, selecting 3 stars per file
   */

  while(fgets(line,MAX,ifp) != NULL) {
    if(sscanf(line,"%s",gname) != 1) {
      fprintf(stderr,"ERROR: Problem reading %s\n",flname);
      return 1;
    }
    if((nstars = read_gsc_file(gname,ofp)) != 3)
      fprintf(stderr,"Warning: Found only %d stars for %s\n",nstars,gname);
  }

  /*
   * Clean up and close
   */

  if(ifp)
    fclose(ifp);
  if(ofp)
    fclose(ofp);

  return 0;
}
 
/*.......................................................................
 *
 * Function read_gsc_file
 *
 * Reads in data from a GSC star file and selects 3 stars for output
 *
 *
 * Inputs: char *gname         name of GSC star file
 *         FILE *ofp           output file pointer
 *
 * Output: int nstars          number of stars found that meet the
 *                              selection criteria
 *
 */

int read_gsc_file(char *gname, FILE *ofp)
{
  int count=0;     /* Number of lines read from GSC star files */
  int nlgsc;       /* Number of lines in GSC star files */
  int nstars=0;    /* Number of good GSC stars found */
  int pa;          /* Position angle of star from field center */
  int hr,min;      /* Position of GSC star */
  float sec;
  int deg,amin;
  float asec;
  float mag;       /* Star magnitude */
  float dist;      /* Distance of star from field center */
  char line[MAX];  /* General string for reading input */
  char junk[MAX];  /* String for reading in useless info */
  FILE *gscfp;     /* Pointer for GSC star files */

  /*
   * Open file
   */

  if((gscfp = fopen(gname,"r")) == NULL) {
    fprintf(stderr,"ERROR: read_gsc_file.  Cannot open %s\n",gname);
    return -1;
  }

  /*
   * Read in number of lines in file
   */

  if((nlgsc = n_lines(gscfp,'!')) <= 0) {
    fprintf(stderr,"ERROR: read_gsc_file. No lines in %s\n",gname);
    return -1;
  }
  else
    rewind(gscfp);

  /*
   * Read in star info from file, skipping SKIPBEG lines at the
   *  beginning and SKIPEND lines at the end, which do not contain
   *  star info
   */

  while(nstars < 3 && fgets(line,MAX,gscfp) != NULL) {
    count ++;
    if(count == (SKIPBEG+1)) {
      if(sscanf(line,"</b>   %s %s %d %d %f %d %d %f %f %s %f %d",
		junk,junk,&hr,&min,&sec,&deg,&amin,&asec,&mag,junk,
		&dist,&pa) != 12) {
	fprintf(stderr,"ERROR: read_gsc_file.  Error in format.\n");
	return -1;
      }
      if(mag < FAINTMAG) {
	nstars++;
	fprintf(ofp,"%sG%d  %02d %02d %05.2f  %+03d %02d %05.2f  %5.2f  ",
		gname,nstars,hr,min,sec,deg,amin,asec,mag);
	fprintf(ofp,"%5.2f %3d\n",dist,pa);
      }
    }
    if(count > (SKIPBEG+1) && (nlgsc-count) >= SKIPEND) {
      if(sscanf(line,"%s %s %d %d %f %d %d %f %f %s %f %d",
		junk,junk,&hr,&min,&sec,&deg,&amin,&asec,&mag,junk,
		&dist,&pa) != 12) {
	fprintf(stderr,"ERROR: read_gsc_file.  Error in format.\n");
	return -1;
      }
      if(mag < FAINTMAG) {
	nstars++;
	fprintf(ofp,"%sG%d  %02d %02d %05.2f  %+03d %02d %05.2f  %5.2f  ",
		gname,nstars,hr,min,sec,deg,amin,asec,mag);
	fprintf(ofp,"%5.2f %3d\n",dist,pa);
      }
    }
  }

  /*
   * Close input file
   */

  if(gscfp)
    fclose(gscfp);

  return nstars;
}

