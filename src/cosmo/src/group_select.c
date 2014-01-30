/*
 * group_select.c
 *
 * This program uses the method from Wilman et al. to determine group
 *  membership from a redshift distribution, based on the method in
 *  Wilman et al. (2005).
 *
 * 09Apr2005 CDF,  First version
 * 02Feb2007 CDF,  Got rid of old distance calculations and replaced with
 *                  newer calc_cosdist library function
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "structdef.h"
#include "dataio.h"
#include "cosmo.h"

typedef struct {
  int id;            /* Galaxy ID */
  double x;          /* x position in arcsec offset from lens */
  double y;          /* y position in arcsec offset from lens */
  double z;          /* Redshift */
  double mag;        /* Magnitude of object */
  double dpos;       /* Offset from centroid */
  int dataflag;      /* Flag */
  char text[MAXC];   /* Label */
} Zcat;          /* Structure for redshift catalog */

Zcat *new_zcat(int size);
Zcat *del_zcat(Zcat *zcat);
Zcat *read_zcat(char *inname, char comment, int *nlines);
void search_radii(double sigma_obs, double z, double Hz, double da, 
		  double *dzmax, double *drmax, double *dthmax);
int select_members(Zcat *zcat, int ncat, Zcat zmed, double dzmax,
		    double dthmax, Zcat *gcat, int *ngroup);
Zcat find_median(Zcat *gcat, int ngroup);
double sigma_gapper(Zcat *gcat, int ngroup);
int doubcmp(const void *v1, const void *v2);



/*.......................................................................
 *
 * Main program
 *
 */

int main(int argc, char *argv[])
{
  int i;                 /* Looping variable */
  int no_error=1;        /* Flag set to 0 on error */
  int count;             /* Number of times through membership loop */
  int contin=1;          /* Flag set to 0 to break loop */
  int contin2;           /* Flag set to 0 to break loop */
  int iddiff;            /* Difference between id's.  If == 0, break loop */
  int nlines;            /* Number of lines in redshift file */
  int ngroup;            /* Number of members in group */
  double z;              /* Group redshift */
  double dz=0.0;         /* Error on lens redshift */
  double sigma_obs;      /* Observed velocity dispersion */
  double r200;           /* Approximate r_200 for given redshift and sigma_v */
  double dzmax;          /* Maximum value of delta(z) */
  double drmax;          /* Maximum value of delta(r) */
  double dthmax;         /* Maximum value of delta(theta) */
  char line[MAX];        /* General string for reading input */
  char zfile[MAXC];      /* File containing redshifts */
  Cosmo cosmo;           /* Cosmological world model */
  Cosdist cosdist;       /* Structure for all of the distance measures */
  Zcat zmed;             /* Centroid */
  Zcat *zcat=NULL;       /* Redshift catalog */
  Zcat *gcat=NULL;       /* Group catalog */
  Zcat *gcat2=NULL;      /* Previous version of the group catalog */
  Zcat *gptr,*gptr2;     /* Pointers to navigate Zcats */

  /*
   * Check the command line invocation
   */

  if(argc < 2) {
    fprintf(stderr,"\nUsage: group_select [zcat] ");
    fprintf(stderr,"\n\n");
    fprintf(stderr," zcat is the file containing the redshifts\n\n");
    return 1;
  }
  printf("\n");

  /*
   * Read in data file
   */

  strcpy(zfile,argv[1]);
  if(!(zcat = read_zcat(zfile,'#',&nlines)))
    no_error = 0;

  /*
   * Make a group catalog of the same size
   */

  if(!(gcat = new_zcat(nlines)))
    no_error = 0;
  if(!(gcat2 = new_zcat(nlines)))
    no_error = 0;

  /*
   * Initialize world model
   */

  cosmo.omega_m = 0.3;
  cosmo.omega_de = 0.7;
  cosmo.w = -1.0;

  /*
   * Get the cosmological parameters
   */

  printf("\n");
  get_cosmo(&cosmo);
 
  /*
   * Loop through groups, if desired
   */

  while(contin && no_error) {

    /*
     * Get redshift
     */

    if(get_valerr(&z,&dz,"estimated group redshift") == 0) {
      fprintf(stderr,"\nERROR.  Exiting group_select.\n\n");
      return 1;
    }
    zmed.z = z;
    printf("Redshift = %f\n",zmed.z);

    /*
     * Compute the distance measures
     */

    cosdist = calc_cosdist(0,z,cosmo);

    /*
     * Compute approximate r_200
     */

    r200 = sigma_obs / ((1 + z) * 11.5 * cosdist.hz);

    /*
     * Print out quantities of interest
     */

    printf("\n--------------------------------------------------\n\n");
    printf("Inputs:\n");
    printf("  z             = %f\n",z);
    printf("  Omega_m       = %f\n",cosmo.omega_m);
    printf("  Omega_Lambda  = %f\n",cosmo.omega_de);
    if(cosmo.omega_de != 0.0)
      printf("  w             = %f\n\n",cosmo.w);
    printf("Calculated quantities:\n");
    printf("   H(z)                      = %5.1f h km/s/Mpc\n",cosdist.hz);
    printf("   Angular diameter distance = %4.0f h^{-1} Mpc\n",
	   cosdist.d_a/MPC2CM);
    printf("   Luminosity distance       = %4.0f h^{-1} Mpc\n",
	   cosdist.d_l/MPC2CM);
    printf("   1.0 h^{-1} Mpc            = %5.2f arcmin\n",
	   RAD2ASEC * MPC2CM/(cosdist.d_a * 60.0));
    printf("   1.0 h^{-1} _comoving_ Mpc = %5.2f arcmin\n",
	   RAD2ASEC * MPC2CM/(cosdist.d_m * 60.0));
    printf("   r_200                     = %4.2f h^{-1} Mpc\n",r200);
    printf("\n--------------------------------------------------\n\n");

    /*
     * Initialize velocity dispersion and centroid
     */

    sigma_obs = 500.0;
    zmed.x = 0.0;
    zmed.y = 0.0;
    count = 1;
    contin2 = 1;

    /*
     * Iterate until convergence 
     */

    while(contin2 && no_error) {

      /*
       * Compute search radii based on redshift and sigma_v
       */

      search_radii(sigma_obs,zmed.z,cosdist.hz,cosdist.d_a,
		   &dzmax,&drmax,&dthmax);

      /*
       * Print out quantities of interest
       */

      printf("\n--------------------------------------------------\n\n");
      printf("Inputs:\n");
      printf("  sigma_obs     = %5.0f\n",sigma_obs);
      printf("  sigma_rest    = %5.0f\n",sigma_obs/(1.0+z));
      printf("  centroid      = (%4.0f,%4.0f)\n",zmed.x,zmed.y);
      printf("Calculated quantities:\n");
      printf("   delta(z)_max              = %6.3f\n",dzmax);
      printf("   delta(r)_max              = %f h^{-1} Mpc\n",drmax);
      printf("   delta(theta_max)          = %7.2f arcsec\n",dthmax);
      printf("\n--------------------------------------------------\n\n");

      /*
       * Select only sources that are within delta(z)_max and 
       *  delta(theta)_max of the centroid
       */

      if(select_members(zcat,nlines,zmed,dzmax,dthmax,gcat,&ngroup) == 1)
	no_error = 0;

      /*
       * Find the velocity dispersion of the group
       */

      if(no_error) {
	if((sigma_obs = sigma_gapper(gcat,ngroup)) < 0.0)
	  no_error = 0;
	else {
	  printf("Group observed velocity dispersion (gapper) = %4.0f km/s\n",
		 sigma_obs);
	  printf("Group rest-frame velocity dispersion (gapper) = %4.0f km/s\n",
		 sigma_obs/(1.0 + z));
	}
      }

      /*
       * Calculate new group centroid
       */

      if(no_error) {
	zmed = find_median(gcat,ngroup);
	if(zmed.z < 0)
	  no_error = 0;
      }

      /*
       * Check to see if membership is stable
       */

      if(no_error) {
	iddiff = 0;
	if(count == 1) {
	  for(i=0,gptr=gcat,gptr2=gcat2; i<ngroup; i++,gptr++,gptr2++)
	    *gptr2 = *gptr;
	}
	else {
	  for(i=0,gptr=gcat,gptr2=gcat2; i<ngroup; i++,gptr++,gptr2++)
	    iddiff += fabs(gptr2->id - gptr2->id);
	  if(iddiff == 0)
	    contin2 = 0;
	  else {
	    printf("Group membership not stable. iddiff = %d\n",iddiff);
	    for(i=0,gptr=gcat,gptr2=gcat2; i<ngroup; i++,gptr++,gptr2++)
	      *gptr2 = *gptr;
	  }
	}
      }

      count++;
    }

    /*
     * Continue?
     */

    printf("\nAnother (group)? (y/n) [y] ");
    fgets(line,MAX,stdin);
    if(line[0] == 'N' || line[0] == 'n')
      contin = 0;
  }

  /*
   * Clean up and exit
   */

  zcat = del_zcat(zcat);
  gcat = del_zcat(gcat);
  gcat2 = del_zcat(gcat2);
  

  printf("\n");
  return 0;
}

/*.......................................................................,
 *
 * Function read_zcat
 *
 * Reads a input file and puts results into a Zcat array,
 *  which has its memory allocation performed in the function.
 * NB:  This function requires input in the form:
 *     ID  x  y  z  mag  
 *
 * Inputs: char *inname        name of input file
 *         char comment        comment character
 *         int *nlines         number of lines in input file -- set by
 *                              this function.
 *         
 * Output: zcat *newdata filled array
 */

Zcat *read_zcat(char *inname, char comment, int *nlines)
{
  int no_error=1;           /* Flag set to 0 on error */
  int ncols;                /* Number of columns in input file */
  char line[MAXC];          /* General input string */
  Zcat *newdata=NULL;       /* Filled zcat array */
  Zcat *dptr;               /* Pointer to navigate zcat */
  FILE *ifp=NULL;           /* Input file pointer */

  /*
   * Open input file
   */

  if(!(ifp = open_readfile(inname))) {
    fprintf(stderr,"ERROR: read_zcat.\n");
    return NULL;
  }

  /*
   * Get number of lines in input file
   */

  if((*nlines = n_lines(ifp,comment)) == 0) {
    fprintf(stderr,"ERROR: read_zcat.  No valid data in input file.\n");
    no_error = 0;
  }
  else {
    rewind(ifp);
  }

  /*
   * Allocate memory for zcat and point dptr at the beginning of
   *  the array.
   */

  if(no_error) {
    if(!(newdata = new_zcat(*nlines)))
      no_error = 0;
    else
      dptr = newdata;
  }

  /*
   * Read in data
   */

  while(no_error && fgets(line,MAXC,ifp) != NULL) {
    if(line[0] != comment) {
      switch(sscanf(line,"%d %lf %lf %lf %lf",&dptr->id,&dptr->x,&dptr->y,
		    &dptr->z,&dptr->mag)) {
      case 5:
	ncols = 5;
	dptr++;
	break;
      default:
	fprintf(stderr,"ERROR: read_zcat.  Bad input format in %s.\n",
		inname);
	fprintf(stderr," Data must be in format %%d %%f %%f %%f %%f\n");
	fprintf(stderr," i.e., integer, float, float, float, float\n");
	no_error = 0;
      }
    }
  }

  /*
   * Clean up and exit
   */

  if(ifp)
    fclose(ifp);

  if(no_error) {
    printf("\nread_zcat: %s has %d columns and %d lines\n",inname,
	   ncols,*nlines);
    return newdata;
  }
  else {
    fprintf(stderr,"ERROR: read_zcat.\n");
    return del_zcat(newdata);
  }
}

/*.......................................................................
 *
 * Function new_zcat
 *
 * Allocates dynamic memory for a pointer array of Zcat structures
 *
 * Input:  int size            size of the array
 *
 * Output: Zcat *newinfo    pointer to the new array.  NULL if error
 *
 */

Zcat *new_zcat(int size)
{
  int i;
  Zcat *newinfo;
  Zcat *dptr;
 
  newinfo = (Zcat *) malloc(sizeof(Zcat) * size);
  if(!newinfo) {
    fprintf(stderr,"new_zcat: \n");
    fprintf(stderr,"Insufficient memory for Zcat array.\n");
    return NULL;
  }

  /*
   * Initialize values in the array
   */

  for(i=0,dptr=newinfo; i<size; i++,dptr++) {
    dptr->x = dptr->y = dptr->z = 0.0;
    dptr->dataflag = 0;
  }

  return newinfo;
}

/*.......................................................................
 *
 * Function del_zcat
 *
 * Frees up memory allocated to Zcat array
 *
 * Input:  Zcat *zcat    array to be freed
 *
 * Output: NULL
 */

Zcat *del_zcat(Zcat *zcat)
{
  if(zcat)
    free(zcat);
 
  return NULL;
}

/*.......................................................................
 *
 * Function search_radii
 *
 * Computes search radii based on redshift and velocity dispersion
 *
 * Inputs: double sigma_obs    observed velocity dispersion
 *         double z            mean group redshift
 *         double Hz           H(z)
 *         double *dzmax       redshit search radius (set by this function)
 *         double *drmax       physical search radius (set by this function)
 *         double *dthmax       angular search radius (set by this function)
 *
 * Output: (none)
 */

void search_radii(double sigma_obs, double z, double Hz, double da,
		  double *dzmax, double *drmax, double *dthmax)
{
  double b=3.5;          /* Factor used in delta(r)_max */

  /*
   * Set delta(z)_max
   */

  *dzmax = 2.0 * sigma_obs / 3.0e5;

  /*
   * Set delta(r)_max
   */

  *drmax = 2.0 * sigma_obs / (b * (1 + z) * Hz);

  /*
   * Set delta(theta)_max
   */

  *dthmax = RAD2ASEC * *drmax / da;
}

/*.......................................................................
 *
 * Function select_members
 *
 * Selects group members based on dzmax and dthmax
 *
 * Inputs: Zcat *zcat          input catalog
 *         int ncat            number of members
 *         Zcat zmed           centroid of position and redshift
 *         double dzmax        maximum offset in redshift space
 *         double dthmax       maximum offset in position
 *         Zcat *gcat          group catalog (refilled by this function)
 *         int *ngroup         number of group members (set by this function)
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int select_members(Zcat *zcat, int ncat, Zcat zmed, double dzmax,
		   double dthmax, Zcat *gcat, int *ngroup)
{
  int i;             /* Looping variable */
  double dz;         /* Redshift offset */
  double dth;        /* Angular offset */
  Zcat *zptr,*gptr;  /* Pointers to zcat and gcat */


  /*
   * Find members
   */

  *ngroup = 0;
  gptr = gcat;

  for(i=0,zptr=zcat; i<ncat; i++,zptr++) {
    dz = fabs(zptr->z - zmed.z);
    dth = sqrt((zptr->x - zmed.x)*(zptr->x - zmed.x) +
	       (zptr->y - zmed.y)*(zptr->y - zmed.y));
    if(dz<dzmax && dth<dthmax) {
      *gptr = *zptr;
      gptr->dpos = dth;
      *ngroup = *ngroup + 1;
      gptr++;
    }
  }
  printf("ngroup = %d\n",*ngroup);

  /*
   * Print out results
   */

  for(i=0,gptr=gcat; i<*ngroup; i++,gptr++) {
    printf("%5d %7.2f %7.2f %7.2f %6.4f  %5.2f\n",gptr->id,gptr->x,gptr->y,
	   gptr->dpos,gptr->z,gptr->mag);
  }
  printf("\n");

  return 0;

}

/*.......................................................................
 *
 * Function sigma_gapper
 *
 * Uses the Beers et al. gapper algorithm to compute a velocity
 *  dispersion for the group.  This algorithm gives
 *                          sqrt(pi)   n-1
 *  sigma(v)_obs = 1.135 c --------- * SUM (w_i * g_i)
 *                         n * (n-1)   i=1
 *
 *  where n is the number of group members, w_i = i * (n - 1)
 *  and g_i = z_(i+1) - z_i
 *
 * NB: The redshifts MUST BE SORTED before the algorithm is used
 *
 * Inputs: Zcat *gcat            group catalog output from select_members
 *         int ngroup            number of members in the group
 *
 * Output: double sigma          velocity dispersion
 *
 */

double sigma_gapper(Zcat *gcat, int ngroup)
{
  int i;                /* Looping variable */
  double magfac;        /* Constant factor in front of sum */
  double sigma=0.0;     /* Velocity dispersion */
  double *zlist;        /* List of redshifts - will be sorted */
  double *zptr;         /* Pointer to navigate zlist */
  Zcat *gptr;           /* Pointer to navigate gcat */

  /*
   * Allocate memory for the redshift list
   */

  if(!(zlist = new_doubarray(ngroup))) {
    fprintf(stderr,"ERROR: sigma_gapper.\n");
    return -1.0;
  }

  /*
   * Transfer redshifts and sort
   */

  for(i=0,zptr=zlist,gptr=gcat; i<ngroup; i++,zptr++,gptr++)
    *zptr = gptr->z;

  qsort(zlist,ngroup,sizeof(double),doubcmp);

  /*
   * Compute the gapper statistic
   */

  magfac = 3.0e5 * 1.135 * sqrt(PI) / (ngroup * (ngroup - 1));
  for(i=0,zptr=zlist; i<ngroup-1; i++,zptr++) {
    sigma += i * (ngroup - i) * (*(zptr+1) - *zptr);
  }
  sigma *= magfac;

  /*
   * Clean up and exit
   */

  zlist = del_doubarray(zlist);
  return sigma;
}

/*.......................................................................
 *
 * Function find_median
 *
 * Finds the median redshift and position for a group
 *
 * Inputs: Zcat *gcat          group catalog
 *         int ngroup          number of group members
 *
 * Output: Zcat medpos         container for median values
 */

Zcat find_median(Zcat *gcat, int ngroup)
{
  int i;               /* Looping variable */
  int medpos;          /* Position of median in sorted array */
  int avgflag;         /* Flag set to 1 if averaging is needed */
  double *temp=NULL;   /* Temporary container for redshift, x, and y */
  double *dptr;        /* Pointer to navigate temp */
  Zcat zmed;           /* Median position */
  Zcat *gptr;          /* Pointer to navigate gcat */

  /*
   * Allocate memory for temporary array
   */

  if(!(temp = new_doubarray(ngroup))) {
    fprintf(stderr,"ERROR: find_median\n");
    zmed.z = -999.0;
    return zmed;
  }

  /*
   * Find index location of median
   */

  if((ngroup/2.0) - floor(ngroup/2.0) == 0.0) {
    avgflag = 1;
    medpos = floor(ngroup/2) - 1;
  }
  else {
    avgflag = 0;
    medpos = floor(ngroup/2.0);
  }

  /*
   * Median z
   */

  for(i=0,gptr=gcat,dptr=temp; i<ngroup; i++,gptr++,dptr++)
    *dptr = gptr->x;
  qsort(temp,ngroup,sizeof(double),doubcmp);
  if(avgflag)
    zmed.x = (*(temp+medpos) + *(temp+medpos+1)) / 2.0;
  else
    zmed.x = *(temp+medpos);
  
  /*
   * Median z
   */

  for(i=0,gptr=gcat,dptr=temp; i<ngroup; i++,gptr++,dptr++)
    *dptr = gptr->y;
  qsort(temp,ngroup,sizeof(double),doubcmp);
  if(avgflag)
    zmed.y = (*(temp+medpos) + *(temp+medpos+1)) / 2.0;
  else
    zmed.y = *(temp+medpos);
  
  /*
   * Median z
   */

  for(i=0,gptr=gcat,dptr=temp; i<ngroup; i++,gptr++,dptr++)
    *dptr = gptr->z;
  qsort(temp,ngroup,sizeof(double),doubcmp);
  if(avgflag)
    zmed.z = (*(temp+medpos) + *(temp+medpos+1)) / 2.0;
  else
    zmed.z = *(temp+medpos);

  /*
   * Give results
   */

  printf("find_median: Median position = (%5.0f,%5.0f)\n",zmed.x,zmed.y);
  printf("find_median: Median redshift = %6.4f\n",zmed.z);

  /*
   * Clean up and exit
   */

  temp = del_doubarray(temp);
  return zmed;
}

/*.......................................................................
 *
 * Function doubcmp
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

int doubcmp(const void *v1, const void *v2)
{
  double *d1 = (double *) v1;  /* double casting of v1 */
  double *d2 = (double *) v2;  /* double casting of v2 */

  /*
   * Do the comparison
   */

  if(*d1 > *d2)
    return 1;
  else if(*d1 == *d2)
    return 0;
  else
    return -1;
}
