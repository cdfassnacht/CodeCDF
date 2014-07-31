#include<stdio.h>
#include<math.h>

/* vlanoise.c
 *
 * Usage: vlanoise
 *
 * Description: Gives rms noise in a VLA map given the band, bandwidth
 *               and integration time.
 *
 */

#define N 27
#define SMALLN 2
#define Fw 1
#define DELSL 0.38
#define DELSC 0.35
#define DELSX 0.29
#define DELSU 1.14
#define DELSK 2.0
#define MAX 1000

main(int argc, char *argv[])
{
  int band,junk;
  float tint,delnu,dels,deli;
  char line[MAX];

  /* Get input */
  printf("Enter the band of the observation (L, C, X, etc.):  ");
  band = getchar();
  while((junk = getchar()) != '\n') 
    ;

  switch(band) {
  case 'l': case 'L':
    dels = DELSL;
    break;
  case 'c': case 'C':
    dels = DELSC;
    break;
  case 'x': case 'X':
    dels = DELSX;
    break;
  case 'u': case 'U':
    dels = DELSU;
    break;
  case 'k': case 'K':
    dels = DELSK;
    break;
  default:
    printf("\nThat band is not supported.  Exiting program\n\n");
    return 0;
  }

  printf("\nEnter the bandwidth in MHz:  ");
  gets(line);
  sscanf(line,"%f",&delnu);
  printf("The bandwidth is %6.1f MHz.\n\n",delnu);
  printf("Enter the integration time in seconds:  ");
  gets(line);
  sscanf(line,"%f",&tint);
  printf("The integration time is %5.0f sec.\n\n",tint);

  /* Make the calculations.
   *
   * The RMS noise is found using equation 19-3 from Lecture 19 of
   *  the 1995 VLA Summer School.  That is,
   *
   *                   Fw * del_S
   *   del_I = ----------------------------
   *          sqrt(n*N*(N-1)*t_int*del_nu/2)
   *
   *   where Fw  =   1.0 for natural weighting and is > 1 for other weightings
   *         del_S = Single interferometer sensitivity per sec per MHz of
   *                     IF bandwidth (value from Lecture 19)
   *         n   =   2 for images of Stokes I from two orthogonal polarization
   *                     states at one sky freq.
   *         N   =   number of antennas
   *         t_int = integration time in sec
   *         del_nu = bandwidth in MHz
   *
   */

  deli = 1000*Fw * dels / sqrt(SMALLN*N*(N-1)*tint*delnu/2);
  /*  Factor of 1000 converts del_I to mJy */

  printf("The theoretical noise with this setup is %f mJy/beam.\n\n",deli);

  return 0;
}

