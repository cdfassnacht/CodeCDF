/*.......................................................................
 *
 * coords.c
 *
 * A series of functions dealing with conversions of sky coordinates
 *  to (x,y) measured in arcsec.
 *
 * v30Jul2001 CDF, added the print_spos function.
 * v08Aug2001 CDF, fixed a bug in the rad2spos function, which 
 *                  was not correctly dealing with the 24/0 RA boundary.
 * v12Jul2003 CDF, Moved new_skypos and del_skypos to structdef.c
 * v29Jul2003 CDF, Moved actual calculation of offsets from a pair of
 *                  positions in radians into the new rad2offset function.
 * v21Jun2007 CDF, Added new ddeg2deg function.
 * v26Jun2007 CDF, Added read_skypos function.
 * v08Jul2007 CDF, Moved main functionality of dspos2spos and ddeg2deg
 *                  into new offset2rad function.
 *                 Added new ddeg2xy function.
 * v04Jul2008 CDF, Added two new functions:
 *                  deg2rad,       to convert decimal degrees to radians
 *                  secat2offset:  to loop through a secat array and calculate
 *                                  offsets from a central position.
 * v24Jan2009 CDF, Fixed a bug in secat2offset
 * v20Mar2013 CDF, Fixed a bug in rad2offset
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "structdef.h"
#include "coords.h"

/*.......................................................................
 *
 * Function print_spos
 *
 * Prints out a sky position in a standard format.
 *
 * Inputs: FILE *fp            pointer to file in which info will be
 *                              printed.  NB: use stdout for standard
 *                              output.
 *         Skypos pos          sky position
 *         
 *
 * Output: (none)
 *
 */

void print_skypos(FILE *fp, Skypos pos)
{
  fprintf(fp,"%02d %02d %07.4f  %+03d %02d %06.3f ",pos.hr,pos.min,pos.sec,
	  pos.deg,pos.amin,pos.asec);
}

/*.......................................................................
 *
 * Function read_skypos
 *
 * Reads in a sky position in either hms or decimal degrees format and
 *  returns the answer as a Pos structure, with its x component set to
 *  alpha and its y component set to delta.
 *
 * Inputs:
 *  (none)                     Run interactively
 *
 * Output:
 *  Pos posdegs                RA and Dec in decimal degrees
 *
 */

Pos read_skypos()
{
  int informat=2;   /* Format of input, either hms or decimal degrees */
  char line[CMAXC]; /* General string for reading input */
  Pos posdegs;      /* Output RA and Dec in decimal degrees */
  Skypos spos;      /* hms-format of input, if needed */

  /*
   * Get input format
   */

  printf("Choose input format for position:\n");
  printf("   1. Decimal degrees\n");
  printf("   2. hms format (default)\n");
  printf("Enter format [%d] ",informat);
  fgets(line,CMAXC,stdin);
  if(line[0] != '\n') {
    while(sscanf(line,"%d",&informat) != 1) {
      fprintf(stderr,"\nERROR: Bad format.  Enter format again");
      fgets(line,CMAXC,stdin);
    }
  }
  if(informat == 1) {
    printf("\nEnter coordinate in decimal degrees format:");
    printf("  ddd.ddddd pddd.dddddd\n");
    printf("   (where p is the sign of the Dec, either + or -)\n");
    fgets(line,CMAXC,stdin);
    while(sscanf(line,"%lf %lf",&posdegs.x,&posdegs.y) != 2) {
      fprintf(stderr,"ERROR: Input does not match input format.  Try again:  ");
      fgets(line,CMAXC,stdin);
    }
  }
  else {
    if(informat != 2)
      printf("Invalid format: Assuming hms format\n");
    printf("\nEnter coordinate in hms format: hh mm ss.sss pdd mm ss.sss\n");
    printf("   (where p is the sign of the Dec, either + or -)\n");
    fgets(line,CMAXC,stdin);
    while(sscanf(line,"%d %d %lf %d %d %lf",
		 &spos.hr,&spos.min,&spos.sec,&spos.deg,&spos.amin,&spos.asec)
	  != 6) {
      fprintf(stderr,"ERROR: Input does not match hms format.  Try again:  ");
      fgets(line,CMAXC,stdin);
    }
    /*
     * Convert to decimal degrees
     */
    spos2deg(spos,&posdegs.x,&posdegs.y);
  }

  return posdegs;
}

/*.......................................................................
 *
 * Function rth2xy
 *
 * Converts a postion in r and theta (in degrees) to a position in (x,y)
 *
 * Inputs: float r             radial distance from center IN ARCSEC
 *         float theta         position angle IN DEGREES
 *         Pos *pos            position in (x,y) (set by this function)
 *         int isastro         flag set to 1 if using astronomical conventions,
 *                              i.e., with the angle measured from N through E,
 *                              and with the E-W axis flipped.
 *
 * Output: None
 *
 */

void rth2xy(double r, double theta, Pos *pos, int isastro)
{
  /*
   * Convert theta to radians
   */

  theta = theta * PI / 180.0;

  /*
   * Take care of the astronomical case
   */

  if(isastro) {
    pos->x = r * sin(theta);
    pos->y = r * cos(theta);
  }

  /*
   * Take care of the mathematical case
   */

  else {
    pos->x = r * cos(theta);
    pos->y = r * sin(theta);
  }
}

/*.......................................................................
 *
 * Function xy2rth
 *
 * Converts a postion in (x,y) to a position in r and theta 
 *
 * Inputs: Pos pos             position in (x,y) 
 *         double *r           radial distance from center IN ARCSEC (set
 *                               by this function)
 *         double *theta       position angle IN DEGREES (set by function)
 *         int isastro         flag set to 1 if using astronomical conventions,
 *                              i.e., with the angle measured from N through E,
 *                              and with the E-W axis flipped.
 *
 * Output: None
 *
 */

void xy2rth(Pos pos, double *r, double *theta, int isastro)
{
  *r = sqrt(pos.x * pos.x + pos.y * pos.y);

  /*
   * Take care of the astronomical case
   */

  if(isastro)
    *theta = atan2(pos.x,pos.y);

  /*
   * Take care of the mathematical case
   */

  else
    *theta = atan2(pos.y,pos.x);

  /*
   * Convert theta to degrees
   */

  *theta = *theta * 180.0 / PI;
}

/*.......................................................................
 *
 * Function deg2rad
 *
 * Converts a sky position (RA,Dec) in decimal degrees to a position in 
 *  radians (ralpha,rdelta)
 *
 * Inputs: 
 *   double dalpha          RA in  decimal degrees
 *   double ddelta          Dec in decimal degrees
 *   double *ralpha         RA in radians (set by this function)
 *   double *rdelta         Dec in radians (set by this function)
 *
 * Output: None
 *
 */

void deg2rad(double dalpha, double ddelta, double *ralpha, double *rdelta)
{
  *ralpha = PI * dalpha / 180.0;
  *rdelta = PI * ddelta / 180.0;
}

/*.......................................................................
 *
 * Function spos2rad
 *
 * Converts a sky position (RA,Dec) to a position in radians (alpha,delta)
 *
 * Inputs: Skypos spos         sky position
 *         double *alpha       RA in radians (set by this function)
 *         double *delta       Dec in radians (set by this function)
 *
 * Output: None
 *
 */

void spos2rad(Skypos spos, double *alpha, double *delta)
{
  *alpha = 15.0 * PI * (spos.hr + spos.min/60.0 + spos.sec/3600.0)/180.0;
  if(spos.deg < 0)
    *delta = PI * (spos.deg - spos.amin/60.0 - spos.asec/3600.0)/180.0;
  else
    *delta = PI * (spos.deg + spos.amin/60.0 + spos.asec/3600.0)/180.0;
}

/*.......................................................................
 *
 * Function rad2spos
 *
 * Converts a position in radians (alpha,dec) into a sky position (RA,Dec)
 *
 * Inputs: double alpha        RA in radians
 *         double delta        Dec in radians
 *         Skypos *spos        sky position (set by this function)
 *
 * Output: None:
 *
 * v08Aug01 CDF, fixed a bug which was not correctly dealing with the 
 *                24/0 RA boundary.
 *
 */

void rad2spos(double alpha, double delta, Skypos *spos)
{
  double alphdeg;     /* alpha in degrees */
  double deltdeg;     /* delta in degrees */

  /*
   * Convert from radians to degrees.
   */

  alphdeg = 180.0 * alpha / (15*PI);
  deltdeg = 180.0 * delta / PI;

  /*
   * Check to make sure that alpha is in the proper range.
   */

  if(alphdeg < 0.0)
    alphdeg += 24.0;
  if(alphdeg >= 24.0)
    alphdeg -= 24.0;

  /*
   * Convert alpha to hms format.
   */

  spos->hr = floor(alphdeg);
  spos->min = floor(60.0 * (alphdeg - spos->hr));
  spos->sec = 3600.0 * (alphdeg - spos->hr - (spos->min/60.0));

  /*
   * Convert delta to hms format.
   */

  if(deltdeg > 0) {
    spos->deg = floor(deltdeg);
    spos->amin = floor(60.0 * (deltdeg - spos->deg));
    spos->asec = 3600.0 * (deltdeg - spos->deg - (spos->amin/60.0));
  }
  else {
    spos->deg = ceil(deltdeg);
    spos->amin = floor(60.0 * (spos->deg - deltdeg));
    spos->asec = 3600.0 * (spos->deg - deltdeg - (spos->amin/60.0));
  }
}

/*.......................................................................
 *
 * Function spos2deg
 *
 * Converts a sky position (RA,Dec) to a position in decimal degrees
 *  (alphadeg,deltadeg)
 *
 * Inputs: Skypos spos         sky position
 *         double *alphadeg    RA in degrees (set by this function)
 *         double *deltadeg    Dec in degrees (set by this function)
 *
 * Output: None
 *
 */

void spos2deg(Skypos spos, double *alphadeg, double *deltadeg)
{
  *alphadeg = 15.0 * (spos.hr + spos.min/60.0 + spos.sec/3600.0);
  if(spos.deg < 0)
    *deltadeg = spos.deg - spos.amin/60.0 - spos.asec/3600.0;
  else
    *deltadeg = spos.deg + spos.amin/60.0 + spos.asec/3600.0;
}

/*.......................................................................
 *
 * Function deg2spos
 *
 * Converts a position in decimal degrees (alpha,dec) into a sky position 
 *  (RA,Dec)
 *
 * Inputs: double alphdeg      RA in decimal degrees
 *         double deltdeg      Dec in decimal degrees
 *         Skypos *spos        sky position (set by this function)
 *
 * Output: None:
 *
 */

void deg2spos(double alphdeg, double deltdeg, Skypos *spos)
{
  alphdeg /= 15.0;
  spos->hr = floor(alphdeg);
  spos->min = floor(60.0 * (alphdeg - spos->hr));
  spos->sec = 3600.0 * (alphdeg - spos->hr - (spos->min/60.0));

  if(deltdeg > 0) {
    spos->deg = floor(deltdeg);
    spos->amin = floor(60.0 * (deltdeg - spos->deg));
    spos->asec = 3600.0 * (deltdeg - spos->deg - (spos->amin/60.0));
  }
  else {
    spos->deg = ceil(deltdeg);
    spos->amin = floor(60.0*(spos->deg - deltdeg));
    spos->asec = 3600.0 * (spos->deg - deltdeg - (spos->amin/60.0));
  }
}

/*.......................................................................
 *
 * Function dspos2xy
 *
 * Takes skypos variables, i.e. R.A. and Dec, for a central position and a 
 *  array of other positions and computes the x and y offsets between them 
 *  in arcsec.  This function uses the formulas in AIPS Memo 27, with the 
 *  SIN geometry.
 *
 * Inputs: Skypos cent         central position position
 *         Skypos *spos        array of other postions
 *         int npos            number of positions in spos array
 *
 * Output: Pos *pos            array of (x,y) offsets.  NULL on error
 *
 * v29Jul03 CDF, Moved actual calculation of offsets into the new
 *                rad2offset function.
 */

Pos *dspos2xy(Skypos cent, Skypos *spos, int npos)
{
  int i;                  /* Looping variable */
  double alpha0,delta0;   /* RA and Dec of central position in radians */
  double alpha2,delta2;   /* RA and Dec of offset positions in radians */
  Pos *pos;               /* Array of x,y offsets */
  Pos *pptr;              /* Pointer to navigate pos */
  Skypos *sptr;           /* Pointer to navigate spos */

  /*
   * Convert central position to radians
   */

  spos2rad(cent,&alpha0,&delta0);

  /*
   * Allocate memory for Pos array
   */

  if(!(pos = new_pos(npos,1))) {
    fprintf(stderr,"ERROR: dspos2xy\n");
    return NULL;
  }

  /*
   * Do the calculations
   */

  for(i=0,pptr=pos,sptr=spos; i<npos; i++,pptr++,sptr++) {

    /*
     * Convert offset position to radians
     */

    spos2rad(*sptr,&alpha2,&delta2);

    /*
     * Compute offset
     */

    *pptr = rad2offset(alpha0,delta0,alpha2,delta2);

  }

  return pos;
}

/*.......................................................................
 *
 * Function ddeg2xy
 *
 * Takes two positions (RA and Dec, in decimal degrees) and computes the 
 *  x and y offsets between them in arcsec.  
 * This function uses the formulas in AIPS Memo 27, with the SIN geometry.
 *
 * Inputs: 
 *  double alphadeg0           central RA in decimal degrees
 *  double deltadeg0           central Dec in decimal degrees
 *  double alphadeg            other RA in decimal degrees
 *  double deltadeg            other Dec in decimal degrees
 *
 * Output: Pos pos             (x,y) offset in arcsec
 *
 */

Pos ddeg2xy(double alphadeg0, double deltadeg0, double alphadeg,
	    double deltadeg)
{
  double alpha0,delta0; /* Central position in decimal degrees */
  double alpha,delta;   /* Other position in decimal degrees */
  Pos pos;              /* Offset between the two positions, in arcsec */

  /*
   * Convert positions to radians
   */

  alpha0 = PI * alphadeg0 / 180.0;
  delta0 = PI * deltadeg0 / 180.0;
  alpha = PI * alphadeg / 180.0;
  delta = PI * deltadeg / 180.0;

  /*
   * Compute offset
   */

  pos = rad2offset(alpha0,delta0,alpha,delta);

  return pos;
}

/*.......................................................................
 *
 * Function rad2offset
 *
 * Takes two input positions, expressed in radians, and returns the
 *  offset between them in arcseconds.  The offset is computed using
 *  AIPS Memo 27 with the SIN geometry
 *
 * Inputs: double alpha1       position 1 right ascension in radians
 *         double delta1       position 1 declination in radians
 *         double alpha2       position 2 right ascension in radians
 *         double delta2       position 2 declination in radians
 *
 * Output: Pos offset          offset in arcseconds
 *
 */

Pos rad2offset(double alpha1, double delta1, double alpha2, double delta2)
{
  double L,M;    /* RA and Dec shifts */
  Pos offset;    /* Offset in arcsec */

  /*
   * Initialize, just in case
   */

  offset.x = 0.0;
  offset.y = 0.0;

  /*
   * Use the formulae from Memo 27
   */

  L = cos(delta2)*sin(alpha2 - alpha1);
  M = sin(delta2)*cos(delta1) - cos(delta2)*sin(delta1)*cos(alpha2-alpha1);

  /*
   * Convert L and M to arcsec and put them into Pos container
   */

  offset.x = L * 180.0 * 3600.0 / PI;
  offset.y = M * 180.0 * 3600.0 / PI;

  return offset;
}

/*.......................................................................
 *
 * Function offset2rad
 *
 * Takes a central position in radians (alpha0,delta0) and an offset
 *  in arcsec (x,y) and converts it into a second position in radians
 *  (alpha,delta) that results from adding the offsets to the first
 *  postion.
 *
 * Inputs:
 *  double alpha0              central RA in radians
 *  double delta0              central Dec in radians
 *  Pos offset                 offset in arcseconds
 *  double *alpha              output RA in radians
 *  double *delta              output Dec in radians
 *
 * Output:
 *  int (0 or 1)               0 ==> success, 1 ==> error
 *
 */

int offset2rad(double alpha0, double delta0, Pos offset, double *alpha,
	       double *delta)
{
  double L,M;         /* Offset positions in radians */
  double L2M2;        /* Defined as sqrt(1 - L^2 - M^2) */

  /*
   * Convert offset to radians
   */

  L = offset.x / 206265.0;
  M = offset.y / 206265.0;

  if((L*L + M*M) > 1.0) {
    fprintf(stderr,"ERROR: offset2rad. Shifts too large\n");
    return 1;
  }
  else
    L2M2 = sqrt(1 - L*L - M*M);

  /*
   * Use Memo 27 formulae to get alpha and delta (positions in radians)
   */

  *alpha = alpha0 + atan(L / (cos(delta0)*L2M2 - M*sin(delta0)));
  *delta = asin(M*cos(delta0) + sin(delta0)*L2M2);

  return 0;
}

/*.......................................................................
 *
 * Function dspos2spos
 *
 * Converts a fixed RA and Dec position plus a series of x,y offsets from
 *  that fixed position to a series of RA and Dec positions.  This
 *  function uses the formulae in AIPS Memo 27 with the SIN geometry
 *
 * Inputs: Skypos cent         central position
 *         Pos *pos            offsets from the central position
 *         int npos            number of offsets
 *
 * Output: Skypos *spos        RAs and Decs off the offset positions,
 *                              NULL on error
 */

Skypos *dspos2spos(Skypos cent, Pos *pos, int npos)
{
  int i;              /* Looping variable */
  double alpha0;      /* RA of center in radians */
  double alpha;       /* RA of offset position in radians */
  double delta0;      /* Dec of center in radians */
  double delta;       /* Dec of offset position in radians */
  Skypos *spos=NULL;  /* Array of RA and Dec positions */
  Skypos *sptr;       /* Pointer to navigate spos array */
  Pos *pptr;          /* Pointer to navigate pos array */

  /*
   * Convert center position to radians
   */

  spos2rad(cent,&alpha0,&delta0);

  /*
   * Allocate memory for Skypos array
   */

  if(!(spos = new_skypos(npos))) {
    fprintf(stderr,"ERROR:  dspos2spos\n");
    return NULL;
  }

  /*
   * Make the conversions
   */

  for(i=0,pptr=pos,sptr=spos; i<npos; i++,pptr++,sptr++) {

    /*
     * Convert central position and offset to new position in radians
     */

    if(offset2rad(alpha0,delta0,*pptr,&alpha,&delta)) {
      fprintf(stderr,"ERROR: dspos2spos.\n");
      return del_skypos(spos);
    }

    /*
     * Convert alpha and delta to standard RA/Dec format
     */

    rad2spos(alpha,delta,sptr);
  }

  return spos;
}

/*.......................................................................
 *
 * Function ddeg2deg
 *
 * Converts a fixed RA and Dec position (in degrees) plus a series of x,y 
 *  offsets (in arcsec) from that fixed position to a series of RA and Dec
 *  positions.  This function uses the formulae in AIPS Memo 27 with the SIN 
 *  geometry
 *
 * Inputs: Pos cent            central position in decimal degrees
 *         Pos *pos            offsets from the central position in arcsec
 *         int npos            number of offsets
 *
 * Output: Pos *adpos          RAs and Decs off the offset positions, in
 *                              decimal degrees.  NULL on error
 */

Pos *ddeg2deg(Pos cent, Pos *pos, int npos)
{
  int i;              /* Looping variable */
  double L,M;         /* Offset positions in radians */
  double L2M2;        /* Defined as sqrt(1 - L^2 - M^2) */
  double alpha0;      /* RA of center in radians */
  double alpha;       /* RA of offset position in radians */
  double delta0;      /* Dec of center in radians */
  double delta;       /* Dec of offset position in radians */
  Pos *adpos=NULL;    /* Array of RA and Dec positions in decimal degrees */
  Pos *aptr;          /* Pointer to navigate adpos array */
  Pos *pptr;          /* Pointer to navigate pos array */

  /*
   * Convert center position to radians
   */

  alpha0 = cent.x * PI / 180.0;
  delta0 = cent.y * PI / 180.0;

  /*
   * Allocate memory for Skypos array
   */

  if(!(adpos = new_pos(npos,1))) {
    fprintf(stderr,"ERROR:  ddeg2deg\n");
    return NULL;
  }

  /*
   * Make the conversions
   */

  for(i=0,pptr=pos,aptr=adpos; i<npos; i++,pptr++,aptr++) {

    /*
     * Convert central position and offset to new position in radians
     */

    if(offset2rad(alpha0,delta0,*pptr,&alpha,&delta)) {
      fprintf(stderr,"ERROR: dspos2spos.\n");
      return del_pos(adpos);
    }

    /*
     * Convert alpha and delta to decimal degrees
     */

    aptr->x = alpha * 180.0 / PI;
    aptr->y = delta * 180.0 / PI;
  }

  return adpos;
}

/*......................................................................
 *
 * Function mod_center
 *
 * Takes a difmap model file and returns the position of the phase
 *  center in (RA,Dec) format.
 *
 * Inputs: FILE *ifp           model file pointer
 *         Skypos *center      RA and Dec of the model center (set by
 *                              this function)
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> failure
 *
 */

int mod_center(FILE *ifp, Skypos *center)
{
  int count=0;       /* Number of central positions read in */
  char line[CMAXC];  /* General string for reading input */
  char word1[CMAXC]; /* First word on an input line */
  Skypos tmp;        /* Temporary container for center position */

  /*
   * Put the file pointer at the beginning of the file
   */

  rewind(ifp);

  /*
   * Read lines from the model file until the line with the model
   *  center is found.
   */

  while(fgets(line,CMAXC,ifp) != NULL && count == 0) {
    if(line[0] == '!') {
      if(sscanf(line,"! %s",word1) != 1) {
	fprintf(stderr,"ERROR: mod_center.  Bad input format.\n");
	return 1;
      }
      if(strcmp(word1,"Center") == 0) {
	if(sscanf(line,"! Center RA: %d %d %lf,  Dec: %d %d %lf",
		  &tmp.hr,&tmp.min,&tmp.sec,&tmp.deg,&tmp.amin,&tmp.asec)
	   != 6)  {
	  fprintf(stderr,"ERROR: mod_center.  Bad input format.\n");
	  return 1;
	}
	count++;
      }
    }
  }

  if(count == 0) {
    fprintf(stderr,"ERROR: mod_center. No model center line in input file\n");
    return 1;
  }
  else {
    *center = tmp;
    return 0;
  }
}

/*.......................................................................
 *
 * Function secat2offset
 *
 * Given a file containing a central position and a Secat array containing
 *  RA and Dec positions (either in decimal degrees or in hms format),
 *  calculate the offsets between the central position and each member
 *  of the Secat array and store the results in the dx, dy, and dpos
 *  members of the Secat array.
 *
 * Inputs:
 *  Secat *secat               Secat array
 *  int ncat                   number of members of the array
 *  int posformat              format of RA/Dec positions in secat, either
 *                              HMS or DEGS
 *  Secat centpos              Secat container for central position
 *  int centformat             format of central position, either HMS or DEGS
 *
 * Output:
 *  int (SUCCESS or ERROR)
 *
 * 2009Jan24 CDF, Fixed a bug 
 */

int secat2offset(Secat *secat, int ncat, int posformat, Secat centpos,
		 int centformat)
{
  int i;                   /* Looping variable */
  int no_error=1;          /* Flag set to 0 on error */
  int ncent;               /* Number of lines in posfile (should be 1) */
  double alpha0,delta0;    /* Central position in radians */
  double alpha_i,delta_i;  /* Catalog object position in radians */
  Pos tmpoffset;           /* Temporary holder for offsets */
  Secat *sptr;             /* Pointer to navigate secat */

  /*
   * Convert the central position to radians
   */

  if(centformat == HMS)
    spos2rad(centpos.skypos,&alpha0,&delta0);
  else
    deg2rad(centpos.alpha,centpos.delta,&alpha0,&delta0);

  /*
   * Loop through the secat array and first convert position to radians
   *  and then call rad2offset to get offset for each object
   */

  for(i=0,sptr=secat; i<ncat; i++,sptr++) {

    /*
     * Convert input position to radians
     */

    if(posformat == HMS)
      spos2rad(sptr->skypos,&alpha_i,&delta_i);
    else
      deg2rad(sptr->alpha,sptr->delta,&alpha_i,&delta_i);

    /*
     * Compute offset
     */

    tmpoffset = rad2offset(alpha0,delta0,alpha_i,delta_i);

    /*
     * Put offset into secat structure
     */
    
    sptr->dx = tmpoffset.x;
    sptr->dy = tmpoffset.y;
    sptr->dpos = sqrt(tmpoffset.x * tmpoffset.x + tmpoffset.y * tmpoffset.y);
  }

  /*
   * Return
   */

  if(no_error)
    return SUCCESS;
  else {
    fprintf(stderr,"ERROR: secat2offset\n");
    return ERROR;
  }
}
