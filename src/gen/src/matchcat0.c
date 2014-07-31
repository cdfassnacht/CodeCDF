/*
 * matchcat.c
 *
 * Usage: matchcat [source root name] [phot system]
 *        The root name is the name used to ID the lens files, 
 *         e.g. '0712' for the CLASS B0712+472 files.
 *        The phot system is either (for now) 'gri' or 'BRI'.
 *        The expected input files are called either 
 *         [source_root]_g.cat, [source_root]_r.cat, [source_root]_i.cat
 *                                   OR
 *         [source_root]_B.cat, [source_root]_R.cat, [source_root]_I.cat
 *
 * This program compares three ascii catalogs that are output from the
 *  SExtractor program.  The catalogs are assumed to be in the g, r, and i
 *  or B, R, and I bands.  This program finds matches between sources in 
 *  the catalogs, based on their positions.
 *
 * 13Nov99 CDF
 * v22Nov99 CDF, Revised to do all three catalogs (g, r, and i) at once,
 *                rather than two at a time.
 * v25Feb00 CDF, Now sort the r-band catalog in x before doing the
 *                matching.  This should make the output file easier to read.
 * v26Feb00 CDF, Read in object flags as well, and add a flag cutoff (>=16)
 *                to reject seriously corrupted magnitudes.
 * v28Feb00 CDF, Now check for SExtractor error value (99.00) in aperture
 *                magnitude in match catalogs.
 * v01Jun00 CDF, Added a source ID to the catalog.
 * v08Jun00 CDF, Split catalog matching into a function that is called for
 *                both the g and i-band matching.
 *               New flag cutoff for g and i-band catalogs is COMPLIM.
 *               Split the output catalog generation into a function.
 *               Added a purge of the r-band (master) catalog to delete
 *                members with fitflag >= MASTERLIM.
 * v10Jan01 CDF, Added option of doing BRI catalog in addition to the
 *                gri catalog.
 *               
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "structdef.h"
#include "dataio.h"

#define MASTERLIM 16
#define COMPLIM 64

/*.......................................................................
 *
 * Structure declarations
 *
 */

typedef struct {
  double x;          /* First coord */
  double y;          /* Second coord */
  float maper;       /* Aperture magnitude */
  float mbest;       /* Total magnitude by "best" method */
  float merr;        /* Error on mbest */
  float fwhm;        /* FWHM of image */
  float class;       /* Star/gal classifier (0.0 --> 1.0 (most starlike)) */
  int fitflag;       /* Flag returned by SExtractor */
  int matchflag;     /* Flag indicating a match has been found */
  float sep;         /* Positional separation between catalogs */
} Secat;             /* Structure for SExtractor output info */

/*.......................................................................
 *
 * Function declarations
 *
 */

Secat *new_secat(int size);
Secat *del_secat(Secat *secat);
Secat *read_secat(char *inname, char comment, int *nlines);
Secat *purge_cat(Secat *in_cat, int nincat, char *catname, int *npurged, 
		 int purgeflag);
int xcmp(const void *v1, const void *v2);
int run_match(Secat *masterval, Secat *compcat, int ncompcat, 
	      Secat *matchptr, double match_thresh, int switchcolor);
int print_matchcat(char *root, char *ext1, char *ext2, char *ext3, 
		   Secat *mastercat, int nmaster, Secat *match12,
		   int nmatch1, Secat *match23, int nmatch2, 
		   int nmatchall);

/*.......................................................................
 *
 * Main program
 *
 */

int main(int argc, char *argv[])
{
  int i;                    /* Looping variable  */
  int no_error=1;           /* Flag set to 0 on error */
  int ncat1,ncat2,ncat3;    /* Number of lines in the three catalogs */
  int nmaster;              /* Number of members of master catalog */
  int nmatch1=0;            /* Number of matches found btwn cat1 and cat2 */
  int nmatch2=0;            /* Number of matches found btwn cat2 and cat3 */
  int mflag1,mflag2;        /* Flags set when matches are found */
  int nmatchall=0;          /* Number of sources in all three catalogs */
  double match_thresh=3.0;  /* Threshold for a valid match */
  char *root;               /* source "root name" for file IDs" */
  char name1[MAXC];         /* Filename for catalog 1 */
  char name2[MAXC];         /* Filename for catalog 2 */
  char name3[MAXC];         /* Filename for catalog 3 */
  char ext1[MAXC];          /* Band name for catalog 1 */
  char ext2[MAXC];          /* Band name for catalog 2 */
  char ext3[MAXC];          /* Band name for catalog 3 */
  char line[MAXC];          /* General string for getting input */
  Secat *cat1=NULL;         /* Data array from catalog 1 */
  Secat *cat2=NULL;         /* Data array from catalog 2 */
  Secat *cat3=NULL;         /* Data array from catalog 3 */
  Secat *mastercat=NULL;    /* Master catalog (purged version of cat2) */
  Secat *match12=NULL;      /* Array containing matches between cats 1 and 2 */
  Secat *match23=NULL;      /* Array containing matches between cats 2 and 3 */
  Secat *mstrptr;           /* Pointer to navigate mastercat */
  Secat *mptr1;             /* Pointer to navigate match12 */
  Secat *mptr2;             /* Pointer to navigate match23 */

  /*
   * Check the command line invocation
   */

  if(argc != 3) {
    fprintf(stderr,"\nUsage: matchcat [source root name] [phot system]\n\n");
    fprintf(stderr,
	    "The root name is the name used to ID the lens files and the\n");
    fprintf(stderr,
	    "  photometric system is either \'gri\' or \'BRI\'.\n\n");
    fprintf(stderr,"For example,"); 
    fprintf(stderr,
	    " if the system were CLASS B0712+472, the root name might be\n");
    fprintf(stderr,
	    " \'0712\' and the phot system might be \'gri\'.\n");
    fprintf(stderr,"In that case, the input files would  be called");
    fprintf(stderr," 0712_g.cat, 0712_r.cat,\n and 0712_i.cat.\n\n");
    return 1;
  }
  printf("\n");

  /*
   * Get the two input filenames from the command line
   */

  root = argv[1];

  /*
   * Create the filenames from the source root name -- this is hardwired
   *  for now.
   */

  if(strcmp(argv[2],"BRI") == 0) {
    sprintf(ext1,"B");
    sprintf(ext2,"R");
    sprintf(ext3,"I");
  }
  else {
    sprintf(ext1,"g");
    sprintf(ext2,"r");
    sprintf(ext3,"i");
  }

  sprintf(name1,"%s_%s.cat",root,ext1);
  sprintf(name2,"%s_%s.cat",root,ext2);
  sprintf(name3,"%s_%s.cat",root,ext3);

  /*
   * Fill the data structures
   */

  if(!(cat1 = read_secat(name1,'#',&ncat1,1)))
    no_error = 0;
  if(!(cat2 = read_secat(name2,'#',&ncat2,1)))
    no_error = 0;
  if(!(cat3 = read_secat(name3,'#',&ncat3,1)))
    no_error = 0;

  /*
   * Allocate memory for the match arrays.  Make them the same
   *  size as cat2, which will be used as the master catalog
   */

  if(!(match12 = new_secat(ncat2)))
    no_error = 0;

  if(!(match23 = new_secat(ncat2)))
    no_error = 0;


  /*
   * Purge catalog 2 of entries with SExtractor flags greater than 
   *  MASTERLIM to create the "master" catalog.
   */

  if(no_error)
    if(!(mastercat = purge_cat(cat2,ncat2,name2,&nmaster,MASTERLIM)))
      no_error = 0;

  /*
   * Before matching, sort the "master" catalog in order of
   *  increasing x, for ease of reading output.
   */

  if(no_error) {
    printf("\nSorting the master catalog in order of increasing x....");
    qsort(mastercat,nmaster,sizeof(mastercat[0]),xcmp);
    printf(" Done.\n");
  }


  /*
   * Set threshold levels
   */

  printf("\nSet maximum separation for a match in pixels: [%4.1f] ",
	 match_thresh);
  fgets(line,MAXC,stdin);
  if(line[0] != '\n') {
    while(sscanf(line,"%lf",&match_thresh) != 1 || match_thresh <= 0.0) {
      fprintf(stderr,"ERROR.  Bad input.  Enter max separation again:  ");
      fgets(line,MAXC,stdin);
    }
  }

  /*
   * Fill match arrays by comparing catalogs
   */

  if(no_error) {

    /*
     *  Loop through the arrays computing positional offsets.
     *  Use catalog 2 as the "master" catalog since the r image is usually
     *   deeper and cleaner than the g and i images.
     */

    for(i=0,mstrptr=mastercat,mptr1=match12,mptr2=match23; i<nmaster; 
	i++,mstrptr++,mptr1++,mptr2++) {

      /*
       * First compare to catalog 1
       */

      if((mflag1 = run_match(mstrptr,cat1,ncat1,mptr1,match_thresh,1)) == 1)
	nmatch1++;

      /*
       * Now compare to catalog 3
       */

      if((mflag2 = run_match(mstrptr,cat3,ncat3,mptr2,match_thresh,-1)) == 1)
	nmatch2++;

      /*
       * Increment nmatchall if matches found in both catalogs.
       */

      if(mflag1 && mflag2)
	nmatchall++;
    }
  }

  /*
   * Write output file
   */

  if(no_error)
    if(print_matchcat(root,ext1,ext2,ext3,mastercat,nmaster,match12,nmatch1,
		      match23,nmatch2,nmatchall))
       no_error = 0;

  /*
   * Clean up and exit
   */

  cat1 = del_secat(cat1);
  cat2 = del_secat(cat2);
  cat3 = del_secat(cat3);
  mastercat = del_secat(mastercat);
  match12 = del_secat(match12);
  match23 = del_secat(match23);

  if(no_error) {
    printf("\nProgram matchcat finished.\n\n");
    return 0;
  }
  else {
    fprintf(stderr,"\nERROR.  Exiting matchcat.\n\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function new_secat
 *
 * Allocates dynamic memory for a pointer array of Secat structures
 *
 * Input:  int size            size of the array
 *
 * Output: Secat *newinfo    pointer to the new array.  NULL if error
 *
 */

Secat *new_secat(int size)
{
  int i;
  Secat *newinfo;
  Secat *sptr;
 
  newinfo = (Secat *) malloc(sizeof(Secat) * size);
  if(!newinfo) {
    fprintf(stderr,"new_secat: \n");
    fprintf(stderr,"Insufficient memory for Secat array.\n");
    return NULL;
  }

  /*
   * Initialize values in the array
   */

  for(i=0,sptr=newinfo; i<size; i++,sptr++) {
    sptr->x = sptr->y;
    sptr->maper = sptr->mbest = -99.0;
    sptr->matchflag = 0;
  }

  return newinfo;
}

/*.......................................................................
 *
 * Function del_secat
 *
 * Frees up memory allocated to Secat array
 *
 * Input:  Secat *secat    array to be freed
 *
 * Output: NULL
 */

Secat *del_secat(Secat *secat)
{
  if(secat)
    free(secat);
 
  return NULL;
}

/*.......................................................................,
 *
 * Function read_secat
 *
 * Reads a the output catalog from a SExtractor run and puts the results
 *  into a Secat array (defined above), which has its memory allocation 
 *  performed in the function.
 * NB:  This function expects input in the form of 16 columns of floats,
 *       of which it reads in columns 1, 2, 3, 5, 6, 14, and 16
 *
 * Inputs: char *inname        name of input file
 *         char comment        comment character
 *         int *nlines         number of lines in input file -- set by
 *                              this function.
 *         
 * Output: Secat *newdata filled array
 */

Secat *read_secat(char *inname, char comment, int *nlines)
{
  int no_error=1;           /* Flag set to 0 on error */
  int ncols;                /* Number of columns in input file */
  float junk;               /* Variable used to read in unrequired data */
  char line[MAXC];          /* General input string */
  Secat *newdata=NULL;      /* Filled secat array */
  Secat *sptr;              /* Pointer to navigate secat */
  FILE *ifp=NULL;           /* Input file pointer */

  /*
   * Open output file
   */

  while(!(ifp = fopen(inname,"r"))) {
    fprintf(stderr,"ERROR: read_secat.  Cannot open %s.\n",inname);
    fprintf(stderr,"Enter new filename:  ");
    fgets(line,MAXC,stdin);
    if(sscanf(line,"%s",inname) != 1) {
      fprintf(stderr,"ERROR: read_secat.  Bad input -- ");
      fprintf(stderr,"setting filename = foo.\n");
      sprintf(inname,"foo");
    }
  }

  /*
   * Get number of lines in input file
   */

  if((*nlines = n_lines(ifp,comment)) == 0) {
    fprintf(stderr,"ERROR: read_secat.  No valid data in input file.\n");
    no_error = 0;
  }
  else {
    rewind(ifp);
  }

  /*
   * Allocate memory for secat and point sptr at the beginning of
   *  the array.
   */

  if(no_error) {
    if(!(newdata = new_secat(*nlines)))
      no_error = 0;
    else
      sptr = newdata;
  }

  /*
   * Read in data
   */

  while(no_error && fgets(line,MAXC,ifp) != NULL) {
    if(line[0] != comment) {
      if((ncols = sscanf(line,
		 "%lf %lf %f %f %f %f %f %f %f %f %f %f %f %f %d %f",
		 &sptr->x,&sptr->y,&sptr->maper,&junk,&sptr->mbest,
		 &sptr->merr,&junk,&junk,&junk,&junk,&junk,&junk,&junk,
		 &sptr->fwhm,&sptr->fitflag,&sptr->class)) != 16) {
	fprintf(stderr,"ERROR: read_secat.  Bad input format in %s.\n",
		inname);
	fprintf(stderr," Data must contain 16 columns.\n");
	no_error = 0;
      }
      else {
	sptr++;
      }
    }
  }

  /*
   * Clean up and exit
   */

  if(ifp)
    fclose(ifp);

  if(no_error) {
    printf("read_secat: %s has %d columns and %d lines\n",inname,
	   ncols,*nlines);
    return newdata;
  }
  else {
    fprintf(stderr,"ERROR: read_secat.\n");
    return del_secat(newdata);
  }
}


/*.......................................................................
 *
 * Function purge_cat
 *
 * Purges a catalog of members with SExtractor flags greater than the
 *  passed purgeflag value to create a new catalog.
 *
 * Inputs: Secat *in_cat       input catalog
 *         int nincat          number of members of the input catalog
 *         char *catname       input catalog name
 *         int *npurged        number of members of the output catalog (set
 *                              by this function)
 *         int purgeflag       flag cutoff value.
 *
 * Output: Secat *out_cat      purged output catalog
 *
 */

Secat *purge_cat(Secat *in_cat, int nincat, char *catname, int *npurged, 
		 int purgeflag)
{
  int i;                 /* Looping variable */
  Secat *out_cat=NULL;   /* Purged output catalog */
  Secat *iptr,*optr;     /* Pointers to navigate the input and output cats */

  /*
   * Initialize
   */

  *npurged = 0;
  printf("\npurge_cat: Purging %s catalog based on flag value....",catname);

  /*
   * Allocate memory for output catalog
   */

  if(!(out_cat = new_secat(nincat))) {
    printf("\n\nERROR: purge_cat\n");
    return NULL;
  }

  /*
   * Loop through the input catalog and select all members with 
   *  fitflag < purgeflag to go into the output catalog.
   */

  optr = out_cat;
  for(i=0,iptr=in_cat; i<nincat; i++,iptr++) {
    if(iptr->fitflag < purgeflag) {
      *optr = *iptr;
      (*npurged)++;
      optr++;
    }
  }

  /*
   * Return
   */

  printf("Done\n");
  printf("purge_cat: New catalog contains %d out of %d original members.\n",
	 *npurged,nincat);
  return out_cat;
}

/*.......................................................................
 *
 * Function xcmp
 *
 * Compares the "x" fields of two Secat structures and returns 1 if
 *  the first x value is greater than the second, 0 if they're equal, and
 *  -1 if the first is less than the second.  This function is called
 *  by qsort.
 *
 * Inputs: void *v1            first structure (cast to void)
 *         void *v2            second structure (cast to void)
 *
 * Output: int (-1,0,1)        as described above
 *
 */

int xcmp(const void *v1, const void *v2)
{
  Secat *s1 = (Secat *) v1;  /* Secat casting of v1 */
  Secat *s2 = (Secat *) v2;  /* Secat casting of v2 */

  /*
   * Do the comparison
   */

  if(s1->x > s2->x)
    return 1;
  else if(s1->x == s2->x)
    return 0;
  else
    return -1;
}

/*.......................................................................
 *
 * Function run_match
 *
 * Searches through a "comparison" catalog for a match to an entry in the
 *  "master" catalog.  The matching is based on the position of the
 *  detected source.  If a match is found, the "match" catalog is updated
 *  to include the entry.
 *
 * Inputs: Secat *masterval    value of entry in master catalog
 *         Secat *compcat      comparison catalog
 *         int ncompcat        number of entries in comparison catalog
 *         Secat *matchptr     pointer to current position in match catalog
 *         double match_thresh threshold for a valid positional match
 *         int switchcolor     constant set to +1 if color is defined as
 *                              (comparison-master) and set to -1 if color
 *                              is defined as (master-comparison)
 *
 * Output: int matchflag       flag set to 1 if a match is found
 *
 */

int run_match(Secat *masterval, Secat *compcat, int ncompcat, 
	      Secat *matchptr, double match_thresh, int switchcolor)
{
  int j;              /* Looping variable */
  int matchflag=0;    /* Flag set to 1 if a match is found */
  double dx,dy;       /* Difference between x and y positions */
  double sep;         /* Position separation */
  double minsep;      /* Minimum separation between catalog positions */
  Secat closest;      /* Closest match between catalog sources */
  Secat *cptr;        /* Pointer to navigate compcat */
  
  /*
   * Run through the comparison catalog looking for the closest 
   *  positional match.  All info from this match is put into the
   *  "closest" structure.
   */

  minsep = 50.0 * match_thresh;
  closest = *masterval;
  matchflag = 0;
  for(j=0,cptr=compcat; j<ncompcat; j++,cptr++) {
    if((dx = cptr->x - masterval->x) < match_thresh &&
       (dy = cptr->y - masterval->y) < match_thresh) {
      sep = sqrt(dx*dx + dy*dy);
      if(sep < minsep) {
	minsep = sep;
	closest = *cptr;
	closest.sep = sep;
      }
    }
  }

  /*
   * If we've got a match, put color and compcat magnitude into the
   *  match catalog.  A valid match is defined by a small separation,
   *  (separation < match_thresh), a valid fitflag (< COMPLIM), and a valid 
   *  aperture magnitude (not equal to 99.00, which is SExtractor's 
   *  error code).
   */

  if(minsep < match_thresh && closest.fitflag < COMPLIM &&
     fabs(closest.maper - 99.0) > 1.0) {
    matchflag = 1;
    *matchptr = *masterval;
    matchptr->matchflag = 1;
    matchptr->sep = closest.sep;

    /*
     * Put relevant info from the comparison catalog into the match 
     *  catalog.  The comparison magnitude goes into the mbest container 
     *  and the color (= switchcolor * (m_aper_compcat - m_aper_master)) 
     *  goes into the maper container.
     */

    matchptr->mbest = closest.mbest;
    matchptr->merr = closest.merr;
    matchptr->maper = switchcolor * (closest.maper - masterval->maper);
    matchptr->fwhm = closest.fwhm;
    matchptr->fitflag = closest.fitflag;
  }

  /*
   * Otherwise, if there is no match, set matchflag to 0
   */

  else {
    matchptr->matchflag = 0;
    matchptr->maper = -99.0;
    matchptr->mbest = -99.0;
    matchptr->merr = 99.0;
    matchptr->fitflag = -1;
  }

  return matchflag;

}

/*.......................................................................
 *
 * Function print_matchcat
 *
 * Prints the matched catalog to an output file.
 *
 * Inputs: char *root          lens system root name
 *         char *ext1          photometric band for catalog 1
 *         char *ext2          photometric band for catalog 2
 *         char *ext3          photometric band for catalog 3
 *         Secat *mastercat    master catalog
 *         int nmaster         number in master catalog
 *         Secat *match12      catalog of matches between catalogs 1 and 2
 *         int nmatch1         number of members in match12
 *         Secat *match23      catalog of matches between catalogs 2 and 3
 *         int nmatch2         number of members in match23
 *         int nmatchall       number of sources found in all three catalogs
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int print_matchcat(char *root, char *ext1, char *ext2, char *ext3, 
		   Secat *mastercat, int nmaster, Secat *match12,
		   int nmatch1, Secat *match23, int nmatch2, 
		   int nmatchall)
{
  int i;                    /* Looping variable */
  int no_error=1;           /* Flag set to 0 on error */
  char outname[MAXC];       /* Output filename */
  Secat *mstrptr;           /* Pointer to navigate mastercat */
  Secat *mptr1;             /* Pointer to navigate match12 */
  Secat *mptr2;             /* Pointer to navigate match23 */
  FILE *ofp=NULL;           /* Pointer to output file */

  /*
   * Open output file
   */

  sprintf(outname,"%s_%s%s%s.cat",root,ext1,ext2,ext3);
  ofp = open_writefile(outname);

  /*
   * Print out results
   */

  if(no_error) {
    printf("\nNumber of valid matches:\n");
    printf("   %d in %s and %s catalogs\n",nmatch1,ext1,ext2);
    printf("   %d in %s and %s catalogs\n",nmatch2,ext2,ext3);
    printf("   %d in all three catalogs\n",nmatchall);
    fprintf(ofp,"#ID     %s_x     %s_y    %sflag  %s_tot   %s_err  %s_aper class ",
	    ext2,ext2,ext2,ext2,ext2,ext2);
    fprintf(ofp,"%sflag    %s     %s_err  (%s-%s)  %sflag    %s     %s_err  (%s-%s)  ",
	    ext1,ext1,ext1,ext1,ext2,ext3,ext3,ext3,ext2,ext3);
    fprintf(ofp,"fwhm_%s fwhm_%s fwhm_%s\n",ext1,ext2,ext3);
    fprintf(ofp,"#--- -------- -------- ----- ------- ------ ------- ----- ");
    fprintf(ofp,"----- ------- ------ ------- ----- ------- ------ ------- ");
    fprintf(ofp,"------ ------ ------\n");
    for(i=0,mstrptr=mastercat,mptr1=match12,mptr2=match23; i<nmaster; 
	i++,mstrptr++,mptr1++,mptr2++) {
      fprintf(ofp,
	      "%04d %8.2f %8.2f %5d %7.3f %6.3f %7.3f %5.2f ",
	      i+1,mstrptr->x,mstrptr->y,mstrptr->fitflag,mstrptr->mbest,
	      mstrptr->merr,mstrptr->maper,mstrptr->class);
      fprintf(ofp,"%5d %7.3f %6.3f %7.3f ",
	      mptr1->fitflag,mptr1->mbest,mptr1->merr,mptr1->maper);
      fprintf(ofp,"%5d %7.3f %6.3f %7.3f %6.2f %6.2f %6.2f\n",
	      mptr2->fitflag,mptr2->mbest,mptr2->merr,mptr2->maper,
	      mptr1->fwhm,mstrptr->fwhm,mptr2->fwhm);
    }
  }

  /*
   * Clean up and exit
   */

  if(ofp)
    fclose(ofp);

  return 0;
}
