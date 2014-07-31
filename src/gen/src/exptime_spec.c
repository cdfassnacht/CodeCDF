/*
 * exptime_spec.c
 *
 * Calculates an estimate of the exposure time needed to achieve a
 *  desired S/N ratio for a spectroscopic observation.
 *
 * Usage: 
 *
 * 16Sep98 CDF
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define PI 3.141592653589793
#define MAXC 80
#define C 3.0e10
#define H 6.627e-27

/*.......................................................................
 *
 * Enumeration for telescope
 *
 */

enum {
  LRIS,
  DBSP,
  P60,
  OTHER
};

/*.......................................................................
 *
 * Function declarations
 *
 */

void get_info(float *fvar, char *prompt);

/*.......................................................................
 *
 * Main program
 *
 */

int main(int argc, char *argv[])
{
  int tel;            /* Telescope/spectrograph combination */
  int npix;           /* Number of pixels per resolution element */
  float slitwidth;    /* Width of slit in arcsec */
  float pixscale=0.0; /* Pixel scale in the spatial direction */
  float rdnoise=0.0;  /* Read noise of CCD */
  float r;            /* Resolution (lambda / Delta lambda) */
  float sn;           /* Desired S/N ratio */
  float dtel=0.0;     /* Diameter of the telescope in m */
  float atel;         /* Area of the telescope in cm^2 */
  float effic=0.0;    /* Overall efficiency of sky + tel + optics + CCD */
  float rmag;         /* R or r band magnitude */
  float sky=21.0;     /* Sky background in mag/arcsec^2 */
  float fsource;      /* Flux from source in erg/cm^2/sec/Hz */
  float fsky;         /* Flux from sky in erg/cm^2/sec/Hz/arcsec^2 */
  float c,x;          /* Temporary variables used in calculations */
  float texp;         /* Estimated exposure time */
  char line[MAXC];    /* General string for reading input */

  /*
   * Get telescope -- just used to set default values
   */

  tel = OTHER;
  printf("\nTelescope/spectrograph choices ");
  printf("(used to set default values only)\n");
  printf("  %d. Keck/LRIS\n",LRIS);
  printf("  %d. P200/DBSP\n",DBSP);
  printf("  %d  P60/Echelle\n",P60);
  printf("  %d. Other\n",OTHER);
  printf(" ----------------\n");
  printf("Enter choice [%d] ",OTHER);
  fgets(line,MAXC,stdin);
  if(line[0] != '\n') {
    while(sscanf(line,"%d",&tel) != 1 || tel < LRIS || tel > OTHER) {
      fprintf(stderr,"ERROR: Bad input.  Enter choice again:  ");
      fgets(line,MAXC,stdin);
    }
  }

  /*
   * Set default values based on telescope/spectrograph
   */

  switch(tel) {
  case LRIS:
    dtel = 10.0;
    atel = PI * 500.0 * 500.0;
    effic = 0.2;
    sky = 22.0;
    rdnoise = 6.5;
    pixscale = 0.215;
    break;
  case DBSP:
    dtel = 5.0;
    atel = PI * 250.0 * 250.0;
    effic = 0.05;
    sky = 21;
    rdnoise = 8.0;
    pixscale = 0.468;
    break;
  default:
  }

  /*
   * Set other general defaults
   */

  sn = 20.0;
  r = 1200.0;
  slitwidth = 1.0;
  npix = 8;

  /*
   * Get other info
   */

  printf("\nEnter r or R magnitude of the object: ");
  fgets(line,MAXC,stdin);
  while(sscanf(line,"%f",&rmag) != 1) {
    fprintf(stderr,"ERROR: Bad input.  Enter mag again:  ");
    fgets(line,MAXC,stdin);
  }

  get_info(&r,"Enter resolution expressed as lambda/dlambda");
  get_info(&sn,"Enter desired signal-to-noise ratio");
  get_info(&slitwidth,"Enter slitwidth in arcsec");

  /*
   * Calculate flux density from source
   */

  fsource = pow(10.0,-0.4*(rmag+48.8));

  /*
   * Calculate flux density from sky.  Note that the sky
   *  magnitudes are entered in mag/arcsec^2, so you have to
   *  multiply by the number of pixels in the resolution element
   *  and the area of each pixel in arcsec^2.  Assume that the
   *  spectral dimension is 1/2 of the slitwidth (i.e. that the
   *  slit projects to two pixels. (??)
   */

  fsky = pow(10.0,-0.4*(sky+48.8)) * pixscale * npix * slitwidth / 2.0;

  /*
   * Set the temporary variable c.  This is the constant, 
   *  non-flux-dependent term.  Also set x to be a repeated
   *  quantity in the equation.
   */

  c = atel * effic / (H * r);
  x = sn * (fsource + fsky);

  /*
   * Calculate the exposure time
   */

  texp = sn * (x + sqrt(x*x - 4*fsource*fsource*npix*rdnoise*rdnoise))
    / (2 * c * fsource * fsource);

  /*
   * Print out results
   */

  printf("\n--------------------------------------------------\n");
  printf("Parameter summary:\n");
  printf("  A_tel        =  %7.1f\n",atel);
  printf("  Readnoise    =  %7.3f\n",rdnoise);
  printf("  n_pix        =  %d\n",npix);
  printf("  Source mag   =  %6.3f\n",rmag);
  printf("  Sky /sq asec =  %6.3f\n",sky);
  printf("Calculated results:\n");
  printf("  S_source     =  %e ergs/cm^2/sec/Hz\n",fsource);
  printf("  S_sky        =  %e ergs/cm^2/sec/Hz\n",fsky);
  printf("  S_sky/sq asec=  %e ergs/cm^2/sec/Hz/arcsec^2\n",
	 fsky*2.0/(pixscale*npix*slitwidth));
  printf("  N_source     =  %e photons\n",fsource*c*texp);
  printf("  N_sky        =  %e photons\n",fsky*c*texp);
  printf("  N_rdnoise    =  %e photons\n",rdnoise*rdnoise*npix);
  printf("--------------------------------------------------\n");

  printf("\nEstimated exposure time is %.0f sec\n",texp);

  return 0;
}

/*.......................................................................
 *
 * Function get_info
 *
 * Handles I/O for float variables, given a pointer to the variable and
 *  a prompt string.
 *
 * Inputs: float *fvar         float variable to be filled.  Should have
 *                              a default value set already.
 *         char *prompt        string used to prompt the user
 *
 * Output: none                function won't be allowed to exit until
 *                              proper input has been received
 */

void get_info(float *fvar, char *prompt)
{
  char line[MAXC];          /* String for getting input */

  printf("\n%s: [%f] ",prompt,*fvar);
  fgets(line,MAXC,stdin);
  if(line[0] == '\n')
    return;
  else {
    while(sscanf(line,"%f",fvar) != 1) {
      fprintf(stderr,"ERROR: get_info.  Bad input.  Enter value again: ");
      fgets(line,MAXC,stdin);
    }
  }
}
