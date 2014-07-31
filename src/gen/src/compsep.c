#include<stdio.h>
#include<string.h>
#include<math.h>

/* compsep.c
 *
 * Usage:  compsep [filename]
 *
 * Description:  Takes as input the *.gmod files produced by the CLASS
 *                automapping script and gives the separation between
 *                the two brightest components as well as the flux ratio
 *                between them.
 *
 * 07Nov94, CDF   A modification of pos.c
 * v26Mar95, CDF  Cleaned things up
 * v03Jul95, CDF  Ignore sources with only one component
 */


#define MAXLINE 1000
#define PI 3.141592653589793

main(int argc, char *argv[])
{
  int i,count;
  float flux; 
  float highflux=-1.0;
  float flux2=-1.0;
  float fluxrat;         /* Flux ratio between brightest and next */
                         /*  brightest components */
  float r, theta;
  float maj,ratio,phi;
  float x, y, highx, highy;
  float x2=0.;
  float y2=0.;
  float sep;             /* Separation (in arcsec) between brightest and */
                         /*  second brightest component */
  char source[MAXLINE];  /* Source name */
  char line[MAXLINE];
  char *cptr;
  FILE *ifp;              /* Input file pointer */


  /* Get input */
  ifp = fopen(argv[1], "r" );

  /* Get source name from input file name */
  cptr = argv[1];
  for(i=0; i<=MAXLINE && *cptr != '.' && *cptr; i++,cptr++)
    source[i] = *cptr;
  source[i] = '\0';

  count=0;
  while (fgets(line, MAXLINE, ifp) != NULL) {
    if (line[0] != '!') {
      /* Initialize */
      maj = ratio = phi = 0.0;
      count++;

      /* change all "v"s to spaces */
      cptr = line;
      for(i=0; i<=MAXLINE && *cptr != '\n' && *cptr; i++,cptr++) {
	if (*cptr == 'v')
	  *cptr = ' ';
      }

      sscanf(line,"%f %f %f %f %f %f",&flux,&r,&theta,&maj,&ratio,&phi);
      
      /* Convert flux to mJy */
      flux = flux * 1000;

      theta = theta * PI / 180;
      x = r * sin(theta) / 1000;
      y = r * cos(theta) / 1000;

      if (flux > highflux) {
	flux2 = highflux;
	x2 = highx;
	y2 = highy;
	highflux = flux;
	highx = x;
	highy = y;
      }
      else if (flux > flux2) {
	flux2 = flux;
	x2 = x;
	y2 = y;
      }
    }
  }

  if(count>1) {
    sep = sqrt(pow((highx-x2),2.) + pow((highy-y2),2.));
    fluxrat = highflux/flux2;
    printf("%s %f %f\n",source, sep, fluxrat);
  }
  fclose(ifp);  

  return 0;
}


