/*
 * find_d.c
 *
 * This program takes the redshifts of a lens and background source and
 *  computes D = D_l D_s / D_ls for a variety of world models.
 *
 * 11Nov2004 CDF, Based very closely on phys_units.c
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "cosmo.h"

#define MAX 80
#define DTH0 2.845397e-6
#define PI 3.141592653589793

/*.......................................................................
 *
 * Function declarations
 *
 */

void calc_d(float zl, float zs, float dzl, float dzs, float theta, 
	    float dtheta, double omega0, double lambda, float *d);
void calc_errs(float zl, float zs, float dzl, float dzs, float theta, 
	       float dtheta, double omega0, double lambda, float dl,
	       float ds, float dls, float d, float me);

/*.......................................................................
 *
 * Main program
 *
 */

int main(int argc, char *argv[])
{
  int i;              /* Loopig variable */
  int readok=0;       /* Flag set to 1 when data successfully read */
  int contin=1;       /* Flag set to 0 to break loop */
  float zl;           /* Lens redshift */
  float zs;           /* Source redshift */
  float dzl=0.0;      /* Error on lens redshift */
  float dzs=0.0;      /* Error on source redshift */
  float d;            /* D = D_l D_s / D_ls */
  double omega0=0.3;  /* Omega_0 */
  double lambda=0.7;  /* Omega_Lambda */
  float theta;        /* Angular separation in arcsec */
  float dtheta;       /* Error in angular separation in arcsec */
  char line[MAX];     /* General string for reading input */
  FILE *ofp=NULL;     /* Output file pointer */

  /*
   * Get redshifts and omega
   */

  readok = 0;
  printf("Enter lens redshift and optional error on redshift:  ");
  while(readok == 0) {
    fgets(line,MAX,stdin);
    switch(sscanf(line,"%f %f",&zl,&dzl)) {
    case 1:
      if(zl < 0.0) {
	fprintf(stderr,"ERROR: improper input.  Enter value(s) again:  ");
	break;
      }
      else {
	readok = 1;
	break;
      }
    case 2:
      if(zl < 0.0 || dzl < 0.0) {
	fprintf(stderr,"ERROR: improper input.  Enter value(s) again:  ");
	break;
      }
      else {
	readok = 1;
	break;
      }
    default:
      fprintf(stderr,"ERROR: improper input.  Enter value(s) again:  ");
    }
  }

  readok = 0;
  printf("Enter source redshift and optional error on redshift:  ");
  while(readok == 0) {
    fgets(line,MAX,stdin);
    switch(sscanf(line,"%f %f",&zs,&dzs)) {
    case 1:
      if(zs < 0.0) {
	fprintf(stderr,"ERROR: improper input.  Enter value(s) again:  ");
	break;
      }
      else {
	readok = 1;
	break;
      }
    case 2:
      if(zs < 0.0 || dzs < 0.0) {
	fprintf(stderr,"ERROR: improper input.  Enter value(s) again:  ");
	break;
      }
      else {
	readok = 1;
	break;
      }
    default:
      fprintf(stderr,"ERROR: improper input.  Enter value(s) again:  ");
    }
  }

  readok = 0;
  printf("Enter max image separation (arcsec) and optional error:  ");
  while(readok == 0) {
    fgets(line,MAX,stdin);
    switch(sscanf(line,"%f %f",&theta,&dtheta)) {
    case 1:
      if(theta < 0.0) {
	fprintf(stderr,"ERROR: improper input.  Enter value(s) again:  ");
	break;
      }
      else {
	readok = 1;
	break;
      }
    case 2:
      if(theta < 0.0 || dtheta < 0.0) {
	fprintf(stderr,"ERROR: improper input.  Enter value(s) again:  ");
	break;
      }
      else {
	readok = 1;
	break;
      }
    default:
      fprintf(stderr,"ERROR: improper input.  Enter value(s) again:  ");
    }
  }

  if((ofp = fopen("d.out","w")) == NULL) {
    fprintf(stderr,"ERROR: Cannot open d.out.\n");
    return 1;
  }
  printf("d.out is now open\n");

  /*
   * Loop through world models
   */

  for(i=0; i<101; i++) {

    /*
     * Run for a flat cosmology
     */

    omega0 = i / 100.0;
    lambda = 1.0 - omega0;

    /*
     * Calculate parameters
     */

    calc_d(zl,zs,dzl,dzs,theta,dtheta,omega0,lambda,&d);
    fprintf(ofp,"%4.2f %4.2f %f\n",omega0,lambda,d);

  }

  fclose(ofp);
  return 0;
}

/*.......................................................................
 *
 * Function calc_params
 *
 * Calculates the physical parameters (angular diameter distances,
 *  Einstein masses, and velocity dispersions) for a lens system
 *  in a given cosmology.  If the redshift errors are non-zero,
 *  call calc_errs to calculate uncertainties.
 *
 * Inputs: float zl            lens redshift
 *         float zs            source redshift
 *         float dzl           error on zl
 *         float dzs           error on zs
 *         float theta         image separation in arcsec
 *         float dtheta        uncertainty in theta
 *         double omega0        Omega_M
 *         double lambda        Omega_Lambda
 *
 * Output: none
 *
 */

void calc_d(float zl, float zs, float dzl, float dzs, float theta, 
	    float dtheta, double omega0, double lambda, float *d)
{
  float dl;           /* Ang. diam. dist. to lens */
  float ds;           /* Ang. diam. dist. to source */
  float dls;          /* Ang. diam. dist. btwn source and lens */
  float me;           /* Mass enclosed in Einstein ring */
  float r_phys;       /* Physical radius of Einstein ring */
  float v_circ;       /* Implied circular velocity, given M_E */
  float sigbn;        /* Velocity dispersion from Blandford & Narayan eqn. */

  /*
   * Compute the angular diameter distances in Mpc
   */

  dl = ang_dist(0.0,zl,omega0,lambda) / MPC2CM;
  ds = ang_dist(0.0,zs,omega0,lambda) / MPC2CM;
  dls = ang_dist(zl,zs,omega0,lambda) / MPC2CM;

  /*
   * Compute D  -- put D in Gpc
   */

  *d = dl * ds / (dls * 1000.0);

}

/*.......................................................................
 *
 * Function calc_errs
 *
 * Calculates the uncertainties in physical parameters derived in
 *  calc_params.
 *
 * Inputs: float zl            lens redshift
 *         float zs            source redshift
 *         float dzl           error on zl
 *         float dzs           error on zs
 *         float theta         image separation in arcsec
 *         float dtheta        uncertainty in theta
 *         double omega0        Omega_M
 *         double lambda        Omega_Lambda
 *         float dl            lens angular diameter distance
 *         float dl            source angular diameter distance
 *         float dl            lens-source angular diameter distance
 *         float d             ratio of distances
 *         float me            Einstein ring mass
 *
 * Output: none
 *
 */

void calc_errs(float zl, float zs, float dzl, float dzs, float theta, 
	       float dtheta, double omega0, double lambda, float dl,
	       float ds, float dls, float d, float me)
{
  float ddl,dds,ddls; /* Errors on D_l, D_s and D_ls */
  float dd;           /* Error on D */
  float dme;          /* Error on M_E */

  /*
   * Do the calculations for angular diameter distances -- answers in Mpc
   */

  ddl = ang_dist_err(0.0,zl,0.0,dzl) / MPC2CM;
  dds = ang_dist_err(0.0,zs,0.0,dzs) / MPC2CM;
  ddls = ang_dist_err(zl,zs,dzl,dzs) / MPC2CM;

  /*
   * Error on D -- answer in Gpc
   */

  dd = d * sqrt(ddl*ddl/(dl*dl) + dds*dds/(ds*ds) + ddls*ddls/(dls*dls));

  /*
   * Error on M_E, answer in solar masses
   */

  dme = me * sqrt(dd*dd/(d*d) + 4.0*dtheta*dtheta/(theta*theta));

  /*
   * Print out answers
   */

  printf("\nUncertainties on the above quantities are:\n");
  printf("   dD_l  = %7.4f h^{-1} Mpc\n",ddl);
  printf("   dD_s  = %7.4f h^{-1} Mpc\n",dds);
  printf("   dD_ls = %7.4f h^{-1} Mpc\n",ddls);
  printf("   dD    = %7.4f h^{-1} Gpc\n",dd);
  printf("   dM_E  = %9.3e h^{-1} M_sun\n\n",dme);

}
