/* nicmos_texp.c
 *
 * Calculates an exposure time for NICMOS given a flux and desired
 *  S/N.
 *
 */

#include <stdio.h>
#include <math.h>

#define RN 30.0
#define ETA 1.86e6
#define B 0.13
#define ID 0.2

int main(int argc, char *argv[])
{
  int n_read;
  float t;
  float flux;
  float cs;
  float sn;
  float sn2;
  float tmp;
  char line[1000];

  printf("Enter expected flux density in Jy/pix:  ");
  gets(line);
  sscanf(line,"%f",&flux);

  printf("Enter number of CCD reads:  ");
  gets(line);
  sscanf(line,"%d",&n_read);

  printf("Enter desired S/N:  ");
  gets(line);
  sscanf(line,"%f",&sn);

  cs = flux * ETA;
  tmp = cs+B+ID;
  sn2 = sn*sn;
  t = (sn2*tmp + sqrt(sn2*sn2*tmp*tmp + 4*sn2*cs*cs*RN*RN/n_read)) / (2*cs*cs);

  printf("\nFor this observation the parameters are:\n");
  printf("  Source flux       = %g Jy/pix\n\n",flux);
  printf("Giving:\n");
  printf("  Source counts     = %g cts\n",cs*t);
  printf("  Sky counts        = %g cts\n",B*t);
  printf("  Dark counts       = %g cts\n",ID*t);
  printf("  Read noise counts = %g cts\n\n",RN*RN/n_read);
  printf("Which means you need an exposure time of %6.0f sec to get\n",t);
  printf(" a S/N of %6.2f:1\n",sn);

  return 0;
}
