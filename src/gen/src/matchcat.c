/*
 * matchcat.c
 *
 * Usage: matchcat [source root name] [epoch_name] [xyz]
 *        The root name is the name used to ID the lens files, 
 *          e.g. '0712' for the CLASS B0712+472 files.
 *        The epoch name is the epoch of the observations, e.g. jan01, comb,
 *         etc.
 *        'xyz' is a 3-character (for now) photometric system, e.g. 'gri' 
 *          or 'BRI'.  The second band (e.g. r or R) is assumed to be the
 *          master catalog.
 *        The expected input files are called 
 *         [source_root]_[x].cat, [source_root]_[y].cat, 
 *          and [source_root]_[z].cat
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
 * v16Oct01 CDF, Added epoch to command-line input.
 *               Changed read_secat, run_match to include new aperture error 
 *                column.
 * v28Jun02 CDF, Added new members of Secat structure and modified the
 *                read_secat function accordingly.
 * v01Jul02 CDF, Changed sorting of output catalog to be in increasing
 *                distance from the lens rather than in increasing x.
 *                Thus, changed xcmp function to dposcmp function.
 *               Moved acquisition of lensx and lensy to the find_lens
 *                function.
 * v02Jul02 CDF, Changed print_matchcat to include dpos field in output.
 * v04Jul02 CDF, Added aperture color information to output, so changed
 *                run_match and print_matchcat functions.
 * v12Jul03 CDF, Moved new_secat and del_secat into structdef.c/structdef.h
 *               Moved read_secat into dataio.c/dataio.h
 * v20Jul03 CDF, Added printout of RA and Dec of catalog sources.
 * v29Jul03 CDF, Moved purge_cat, find_lens, and dposcmp functions into
 *                the new catlib.c/catlib.h library files.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "structdef.h"
#include "dataio.h"
#include "catlib.h"


/*.......................................................................
 *
 * Function declarations
 *
 */

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
  int nband=0;              /* Number of bands */
  int ncat1,ncat2,ncat3;    /* Number of lines in the three catalogs */
  int nmaster;              /* Number of members of master catalog */
  int lensindex;            /* Catalog index of closest match to lens pos */
  int nmatch1=0;            /* Number of matches found btwn cat1 and cat2 */
  int nmatch2=0;            /* Number of matches found btwn cat2 and cat3 */
  int mflag1,mflag2;        /* Flags set when matches are found */
  int nmatchall=0;          /* Number of sources in all three catalogs */
  double match_thresh=3.0;  /* Threshold for a valid match */
  char *root;               /* source "root name" for file IDs" */
  char name1[MAXC];         /* Filename for catalog 1 */
  char name2[MAXC];         /* Filename for catalog 2 */
  char name3[MAXC];         /* Filename for catalog 3 */
  char epoch[MAXC];         /* Observing epoch */
  char ext[MAXC];           /* Array to hold band names */
  char mext;                /* Master catalog band name */
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

  if(argc < 4) {
    fprintf(stderr,"\nUsage: matchcat [source root name] [epoch] [xyz]\n\n");
    fprintf(stderr,
	    "The root name is the name used to ID the lens files.\n");
    fprintf(stderr,
	    "The epoch is the observing epoch, e.g., jan01, comb, etc.\n");
    fprintf(stderr,"\'xyz\' is a 3-character photometric system, e.g. ");
    fprintf(stderr,"\'gri\' or 'BRI'.\n");
    fprintf(stderr,"The second band (e.g. r or R) is assumed to be the ");
    fprintf(stderr,"master catalog.\n\n");
    fprintf(stderr,"For example,"); 
    fprintf(stderr,
	    " if the system were CLASS B0712+472, the root name might be\n");
    fprintf(stderr,
	    " \'0712\', the epoch might be jan01 and the phot system might ");
    fprintf(stderr,"be \'BRI\'.\n");
    fprintf(stderr,"In that case, the input files would  be called:");
    fprintf(stderr," 0712_jan01_B.cat, 0712_jan01_R.cat, and");
    fprintf(stderr," 0712_jan01_I.cat.\n\n");
    return 1;
  }
  printf("\n");

  /*
   * Get the root name from the command line.
   */

  root = argv[1];
#if 0
  /*
   * Get the master band name
   */

  if(no_error) {
    mext = *argv[3];
    printf("\nMaster catalog is %c\n",mext);
  }
#endif
  /*
   * Get the bands for the catalogs from the third input string
   */

  if((nband = strlen(argv[3])) != 3)
    no_error = 0;
  else {
    printf("\nCatalogs in %d bands:\n",nband);
    for(i=0; i<nband; i++) {
      ext[i] = *(argv[3]+i);
      printf("  %c\n",ext[i]);
    }
    printf("\n");
  }
  /*
   * Create the filenames from the source root name
   */

  if(no_error) {
    sprintf(ext1,"%c",ext[0]);
    sprintf(ext2,"%c",ext[1]);
    sprintf(ext3,"%c",ext[2]);
    sprintf(name1,"%s_%s_%s.cat",root,argv[2],ext1);
    sprintf(name2,"%s_%s_%s.cat",root,argv[2],ext2);
    sprintf(name3,"%s_%s_%s.cat",root,argv[2],ext3);

  /*
   * Fill the data structures
   */

    if(!(cat1 = read_secat(name1,'#',&ncat1,1)))
      no_error = 0;
    if(!(cat2 = read_secat(name2,'#',&ncat2,1)))
      no_error = 0;
    if(!(cat3 = read_secat(name3,'#',&ncat3,1)))
      no_error = 0;
  }

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
   * Find lens system in master catalog.
   */

  if(no_error)
    if(find_lens(mastercat,nmaster,&lensindex))
      no_error = 0;

  /*
   * Sort the "master" catalog in order of
   *  increasing distance from the lens.
   */

  if(no_error) {
    printf("\nSorting the master catalog in order of increasing distance ");
    printf("from lens...");
    qsort(mastercat,nmaster,sizeof(mastercat[0]),dposcmp);
    printf(" Done.\n");
  }

  /*
   * Set threshold levels
   */

  if(no_error) {
    printf("\nSet maximum separation for a match in pixels: [%4.1f] ",
	   match_thresh);
    fgets(line,MAXC,stdin);
    if(line[0] != '\n') {
      while(sscanf(line,"%lf",&match_thresh) != 1 || match_thresh <= 0.0) {
	fprintf(stderr,"ERROR.  Bad input.  Enter max separation again:  ");
	fgets(line,MAXC,stdin);
      }
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

  if(no_error) {
    sprintf(name1,"%s_%s",root,argv[2]);
    if(print_matchcat(name1,ext1,ext2,ext3,mastercat,nmaster,match12,nmatch1,
		      match23,nmatch2,nmatchall))
       no_error = 0;
  }

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
   *  isophotal magnitude (not equal to 99.00, which is SExtractor's 
   *  error code).
   */

  if(minsep < match_thresh && closest.fitflag < COMPLIM &&
     fabs(closest.miso - 99.0) > 1.0) {
    matchflag = 1;
    *matchptr = *masterval;
    matchptr->matchflag = 1;
    matchptr->sep = closest.sep;

    /*
     * Put relevant info from the comparison catalog into the match 
     *  catalog.  The comparison magnitude goes into the mtot container 
     *  and the aperture colors (= switchcolor * (m*_compcat - m*_master)) 
     *  goes into the miso and ma* containers.
     */

    matchptr->mtot = closest.mtot;
    matchptr->merr = closest.merr;
    matchptr->misoerr = closest.misoerr;
    matchptr->miso = switchcolor * (closest.miso - masterval->miso);
    matchptr->ma1 = switchcolor * (closest.ma1 - masterval->ma1);
    matchptr->ma2 = switchcolor * (closest.ma2 - masterval->ma2);
    matchptr->ma3 = switchcolor * (closest.ma3 - masterval->ma3);
    matchptr->fwhm = closest.fwhm;
    matchptr->fitflag = closest.fitflag;
  }

  /*
   * Otherwise, if there is no match, set matchflag to 0
   */

  else {
    matchptr->matchflag = 0;
    matchptr->miso = -99.0;
    matchptr->ma1 = -99.0;
    matchptr->ma2 = -99.0;
    matchptr->ma3 = -99.0;
    matchptr->mtot = -99.0;
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
 * v20Jul03 CDF, Added printout of RA and Dec of catalog sources.
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
    fprintf(ofp,"#ID  dpos    %s_x     %s_y    class ",ext2,ext2);
    fprintf(ofp,"%sflg  %s_tot %s_err %s_iso  %s_ap2  fwhm_%s ",
	    ext2,ext2,ext2,ext2,ext2,ext2);
    fprintf(ofp,"%sflg  %s_tot %s_err (%s-%s)iso (%s-%s)ap2 fwhm_%s ",
	    ext1,ext1,ext1,ext1,ext2,ext1,ext2,ext1);
    fprintf(ofp,"%sflg  %s_tot %s_err (%s-%s)iso (%s-%s)ap2 fwhm_%s ",
	    ext3,ext3,ext3,ext2,ext3,ext2,ext3,ext3);
    fprintf(ofp," RA2000   Dec2000\n");
    fprintf(ofp,"#--- ---- -------- -------- ----- ");
    fprintf(ofp,"---- ------ ----- ------ ------ ------ ");
    fprintf(ofp,"---- ------ ----- -------- -------- ------ ");
    fprintf(ofp,"---- ------ ----- -------- -------- ------\n");
    for(i=0,mstrptr=mastercat,mptr1=match12,mptr2=match23; i<nmaster; 
	i++,mstrptr++,mptr1++,mptr2++) {
      fprintf(ofp,
	      "%04d %4.0f %8.2f %8.2f %5.2f ",
	      i+1,mstrptr->dpos,mstrptr->x,mstrptr->y,mstrptr->class);
      fprintf(ofp,"%3d  %6.2f %5.2f %6.2f %6.2f %6.2f ",
	      mstrptr->fitflag,mstrptr->mtot,mstrptr->merr,mstrptr->miso,
	      mstrptr->ma2,mstrptr->fwhm);
      fprintf(ofp,"%3d  %6.2f %5.2f  %6.2f   %6.2f  %6.2f ",
	      mptr1->fitflag,mptr1->mtot,mptr1->merr,mptr1->miso,
	      mptr1->ma2,mptr1->fwhm);
      fprintf(ofp,"%3d  %6.2f %5.2f  %6.2f   %6.2f  %6.2f  ",
	      mptr2->fitflag,mptr2->mtot,mptr2->merr,mptr2->miso,
	      mptr2->ma2,mptr2->fwhm);
      fprintf(ofp,"%02d %02d %07.4f %+03d %02d %06.3f\n",
	      mptr2->skypos.hr,mptr2->skypos.min,mptr2->skypos.sec,
	      mptr2->skypos.deg,mptr2->skypos.amin,mptr2->skypos.asec);
    }
  }

  /*
   * Clean up and exit
   */

  if(ofp)
    fclose(ofp);

  return 0;
}
