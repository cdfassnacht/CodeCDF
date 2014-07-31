/* nicmos_sn.c
 *
 * Calculates a S/N ratio for NICMOS given a flux
 *
 */

#include <stdio.h>
#include <math.h>

#define SKY 1.0e-17
#define DARK 0.007
#define RN 30.0
#define PTSENS 9.6e14
#define DIFFAC 0.05
#define ETA 1.86e6
#define B 0.13
#define ID 0.2

int main(int argc, char *argv[])
{
  int npix;
  int n_read;
  float t;
  float flux;
  float cs;
  float bsky;
  float sn;
  char line[1000];

  printf("Enter expected flux density in Jy/pix:  ");
  gets(line);
  sscanf(line,"%f",&flux);

  printf("Enter number of CCD reads:  ");
  gets(line);
  sscanf(line,"%d",&n_read);

  printf("Enter exposure time in sec:  ");
  gets(line);
  sscanf(line,"%f",&t);

  cs = flux * ETA;
  sn = cs * t / sqrt((cs+B+ID)*t + RN*RN/n_read);

  printf("\nFor this observation the parameters are:\n");
  printf("  Flux              = %g Jy/pix\n",flux);
  printf("  Source counts     = %g cts\n",cs*t);
  printf("  Sky counts        = %g cts\n",B*t);
  printf("  Dark counts       = %g cts\n",ID*t);
  printf("  Read noise counts = %g cts\n\n",RN*RN/n_read);
  printf("Which gives a S/N of %6.2f:1 in %6.0f sec\n",sn,t);

  return 0;
}
