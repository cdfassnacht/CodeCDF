/* hst_sn.c
 *
 * Calculates a S/N ratio for STIS given a WF/PC-2 flux
 *
 */

#include <stdio.h>
#include <math.h>

#define SKY 1.0e-17
#define DARK 0.007
#define RN 4.0
#define PTSENS 9.6e14
#define DIFFAC 0.05

int main(int argc, char *argv[])
{
  int npix_lambda;
  int npix_spat;
  int npix;
  int n_read;
  int n_bin;
  float t;
  float wfpcflux;
  float sens;
  float slit_width;
  float cs;
  float bsky;
  float sn;
  char line[1000];

  printf("Enter WF/PC-2 flux in ergs/sec/cm^2/Ang/arcsec^2:  ");
  gets(line);
  sscanf(line,"%f",&wfpcflux);

  printf("Enter STIS slit width in arcsec:  ");
  gets(line);
  sscanf(line,"%f",&slit_width);

  printf("Enter number of CCD reads:  ");
  gets(line);
  sscanf(line,"%d",&n_read);

  printf("Enter number of spectral pixels in a resolution element:  ");
  gets(line);
  sscanf(line,"%d",&npix_lambda);

  printf("Enter number of spatial pixels in a resolution element:  ");
  gets(line);
  sscanf(line,"%d",&npix_spat);

  printf("Enter binning:  ");
  gets(line);
  sscanf(line,"%d",&n_bin);

  printf("Enter exposure time in sec:  ");
  gets(line);
  sscanf(line,"%f",&t);

  npix = npix_lambda * npix_spat;
  sens = PTSENS * DIFFAC * slit_width;
  cs = wfpcflux * sens * npix;
  bsky = SKY * sens;
  sn = cs * t / 
    sqrt(cs*t + (bsky+DARK)*npix*t + npix*n_read*RN*RN/n_bin);

  printf("\nFor this observation the parameters are:\n");
  printf("  WF/PC-2 Flux      = %g erg/sec/cm^2/Ang/arcsec^2\n",wfpcflux);
  printf("  Sensitivity       = %g\n\n",sens);
  printf("So in one resolution element of %d x %d pixels:\n",npix_lambda,
	 npix_spat);
  printf("  Source counts     = %g cts\n",cs*t);
  printf("  Sky counts        = %g cts\n",bsky*npix*t);
  printf("  Dark counts       = %g cts\n",DARK*npix*t);
  printf("  Read noise counts = %g cts\n\n",npix*n_read*RN*RN/n_bin);
  printf("Which gives a S/N of %6.2f:1 in %6.0f sec\n",sn,t);

  return 0;
}
