/*
 * catlib.c
 *
 * This is a library of functions to process catalog (Secat) structures.
 *  purge_cat       - deletes members of a catalog with flag values
 *                     greater than a given limit
 *  find_closest    - finds the closest catalog member to the given (x,y)
                       position
 *  find_lens       - finds the closest source in the catalog to a given 
 *                     position
 *  dposcmp         - compares the dpos members of two Secat structures --
 *                     used in sorting catalogs.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "structdef.h"
#include "dataio.h"
#include "catlib.h"

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
  printf("\npurge_cat: Purging %s catalog based on flag value....\n",catname);
  printf("purge_cat: Rejecting sources with flags >= %d\n",purgeflag);

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

  printf("   ....Done\n");
  printf("purge_cat: New catalog contains %d out of %d original members.\n",
	 *npurged,nincat);
  return out_cat;
}

/*.......................................................................
 *
 * Function find_lens
 *
 * Finds the "true" lens position in the catalog given the rough lens
 *  position.
 *
 * Inputs: Secat *cat          catalog in which to search
 *         int ncat            number of catalog entries
 *         int *lensindex      catalog index of closest match to lens
 *                              position (set by this function)
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error.
 *
 * v01Aug2003 CDF, Added a returned value of lensindex, corresponding
 *                  to the catalog index value for the lens.
 */

int find_lens(Secat *cat, int ncat, int *lensindex)
{
  int i;              /* Looping variable */
  double lensx,lensy; /* Lens position */
  double dx,dy,dpos;  /* Offsets from lens position */
  double dposmin;     /* Temporary minimum offset from lens position */
  char line[MAXC];    /* General string for getting input. */
  Secat tmppos;       /* Temporary best match to lens position */
  Secat *sptr;        /* Pointer to navigate cat */

  /*
   * Get lens position
   */

  printf("\nEnter position of lens system [x y]: ");
  fgets(line,MAXC,stdin);
  while(sscanf(line,"%lf %lf",&lensx,&lensy) != 2) {
    fprintf(stderr,"ERROR.  Bad input.  Enter lens position again:  ");
    fgets(line,MAXC,stdin);
  }

  /*
   * Initialize
   */

  dx = cat->x - lensx;
  dy = cat->y - lensy;
  dposmin = sqrt(dx*dx + dy*dy);
  tmppos = *cat;
  *lensindex = 0;

  /*
   * Do the search by looping through the catalog
   */

  for(i=0,sptr=cat; i<ncat; i++,sptr++) {
    dx = sptr->x - lensx;
    dy = sptr->y - lensy;
    if((dpos = sqrt(dx*dx + dy*dy)) < dposmin) {
      dposmin = dpos;
      tmppos = *sptr;
      *lensindex = i;
    }
  }

  /*
   * Query the user about the best match
   */

  printf("\nfind_lens: Closest object to lens position is %6.2f pixels ",
	 dposmin);
  printf("away.\n");
  printf("Is this acceptable? [y] ");
  fgets(line,MAXC,stdin);
  if(line[0] == 'n' || line[0] == 'N') {
    printf("\n*** WARNING -- Lens not found in master catalog.\n");
    printf("Keeping rough position for lens (%7.2f %7.2f).\n",lensx,lensy);
    *lensindex = 0;
  }
  else {
    lensx = tmppos.x;
    lensy = tmppos.y;
    printf("\nfind_lens: Setting lens position to (%7.2f %7.2f).\n",
	   lensx,lensy);
  }

  /*
   * Put displacements from actual lens position into the catalog
   */

  for(i=0,sptr=cat; i<ncat; i++,sptr++) {
    sptr->dx = sptr->x - lensx;
    sptr->dy = sptr->y - lensy;
    sptr->dpos = sqrt(sptr->dx * sptr->dx + sptr->dy * sptr->dy);
  }

  return 0;
}


/*.......................................................................
 *
 * Function find_closest
 *
 * Finds the object giving the closest match between the input (x,y) postion
 *  and  the members of a SExtractor catalog.
 *  The closest matching positions are returned in a Secat structure.
 *
 * Inputs: Pos cpos            input position
 *         Secat *secat        catalog based on fits image (x,y) positions
 *         int ncat            number of members in secat
 *         int verbose         flag set to 1 for verbose output
 *
 * Output: Secat bestmatch     best match
 *
 */

Secat find_closest(Pos cpos, Secat *secat, int ncat, int verbose)
{
  int i,j;                /* Looping variables */
  double dx,dy,dpos;      /* Offsets between (x,y) pairs */
  double mindpos;         /* Running minimum position offset */
  Secat bestmatch;        /* Container for best-matching xy object */
  Secat *sptr;            /* Pointers to navigate astcat and xycat */

  /*
   * Loop through catalog
   */

  bestmatch = *secat;
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

  if(verbose) {
    printf("For this position, the best match is ID=%d at",
	   bestmatch.id);
    printf(" mindpos = %6.3f pix\n",bestmatch.dpos);
  }

  /*
   * Exit
   */

  return bestmatch;
}

/*.......................................................................
 *
 * Function dposcmp
 *
 * Compares the dpos field for two Secat structures and returns 1
 *  if the first is greater than the second, 0 if they're equal, and
 *  -1 if the first is less than the second.  This function is called
 *  by qsort.
 *
 * Inputs: void *v1            first structure (cast to void)
 *         void *v2            second structure (cast to void)
 *
 * Output: int (-1,0,1)        as described above
 *
 */

int dposcmp(const void *v1, const void *v2)
{
  Secat *s1 = (Secat *) v1;  /* Secat casting of v1 */
  Secat *s2 = (Secat *) v2;  /* Secat casting of v2 */

  /*
   * Do the comparison
   */

  if(s1->dpos > s2->dpos)
    return 1;
  else if(s1->dpos == s2->dpos)
    return 0;
  else
    return -1;
}

/*.......................................................................
 *
 * Function dposcmp_sdss
 *
 * Compares the dpos field for two SDSScat structures and returns 1
 *  if the first is greater than the second, 0 if they're equal, and
 *  -1 if the first is less than the second.  This function is called
 *  by qsort.
 *
 * Inputs: void *v1            first structure (cast to void)
 *         void *v2            second structure (cast to void)
 *
 * Output: int (-1,0,1)        as described above
 *
 */

int dposcmp_sdss(const void *v1, const void *v2)
{
  SDSScat *s1 = (SDSScat *) v1;  /* SDSScat casting of v1 */
  SDSScat *s2 = (SDSScat *) v2;  /* SDSScat casting of v2 */

  /*
   * Do the comparison
   */

  if(s1->dpos > s2->dpos)
    return 1;
  else if(s1->dpos == s2->dpos)
    return 0;
  else
    return -1;
}

/*.......................................................................
 *
 * Function dcmp
 *
 * Compares two doubles and returns 1
 *  if the first is greater than the second, 0 if they're equal, and
 *  -1 if the first is less than the second.  This function is called
 *  by qsort.
 *
 * Inputs: void *v1            first structure (cast to void)
 *         void *v2            second structure (cast to void)
 *
 * Output: int (-1,0,1)        as described above
 *
 */

int dcmp(const void *v1, const void *v2)
{
  double *s1 = (double *) v1;  /* double casting of v1 */
  double *s2 = (double *) v2;  /* double casting of v2 */

  /*
   * Do the comparison
   */

  if(*s1 > *s2)
    return 1;
  else if(*s1 == *s2)
    return 0;
  else
    return -1;
}

