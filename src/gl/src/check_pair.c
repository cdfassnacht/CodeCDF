/*
 * ray_trace.c
 *
 * Usage: check_pair [model_params_file] [image_pos_file]
 *
 * Description:  Uses a known potential and observed image positions
 *                to see if the images could be produced by a background
 *                source with redshift factor zf.  Finds the redshift
 *                factor for which rays through the images come closest
 *                to each other when passed back through the potential.
 *
 * 21Mar96, CDF Modification of ray_trace.c
 * v29Apr96, CDF Improve file I/O to match that in ray_trace and model_fit
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "structdef.h"
#include "file_io.h"
#include "list_tools.h"
#include "lensprogs.h"
#include "plotfuncs.h"

#define PI 3.141592653589793
#define D2R PI/180
#define MAXLINE 1000
#define NSTEPS 100

int main(int argc, char *argv[])
{
  int i,j;
  int nobs=0;             /* Number of observed images */
  int nmod=0;             /* Number of model types in lens potential */
  float zfmin,zfmax;      /* Range for redshift factors to be tried */
  float zf;               /* Redshift scaling factor */
  float delx,dely;        /* Distance offsets in x and y */
  float tmpdel,del;       /* Total distance offsets */
  char line[MAXLINE];     /* General character array for reading in info */
  char name[MAXLINE];     /* Source name */
  Lensparams *params;     /* Model potential parameters */
  Pos *sp,*xptr;          /* Source positions and a general pointer */
  Posinfo *piobs,*piptr;  /* Image positions and a general pointer */

  /*
   * Check input on command line.
   */

  if (argc < 3) {
    fprintf(stderr,"Usage:  check_pair [model_params_file] [image_pos_file]\n");
    return 1;
  }

  /*
   * Read image positions, fluxes and redshift factors
   *   from the first input file.
   */

  if(!(piobs = posinfo_io(argv[1],&nobs)))
    return 1;

  /*
   * Read source name and lens parameters from second input file.
   */

  if(!(params = lp_io(argv[2],name,&nmod)))
    return 1;

  if(nobs != 2) {
    fprintf(stderr,"Currently this program only works with image pairs.\n");
    return 1;
  }
  
  /*
   * Calculate source position based on lens parameters, potentials
   *  and image positions.
   */

  printf("Enter min value for redshift factor range:  ");
  gets(line);
  sscanf(line,"%f",&zfmin);
  printf("Enter max value for redshift factor range:  ");
  gets(line);
  sscanf(line,"%f",&zfmax);

  del = 1.0e7;
  for(i=0;i<NSTEPS;i++) {
    zf = zfmin + i*(zfmax-zfmin)/NSTEPS;
    if(!(sp=new_pos(nobs,1)))
      return 1;
    xptr = sp;
    piptr = piobs;
    for(j=0;j<nobs;j++,xptr++,piptr++) {
      piptr->zf = zf;
      if(rough_source(xptr,piptr,1,params,nmod,0)) {
	fprintf(stderr,"Exiting program\n");
	return 1;
      }
    }
    delx = sp->x - sp[1].x;
    dely = sp->y - sp[1].y;
    tmpdel = sqrt(delx*delx + dely*dely);
    if (tmpdel < del) {
      del = tmpdel;
      printf("%d. For zf = %6.4f source positions are %g pixels apart.\n",
	     i,zf,del);
    }
    else
      printf("%d. zf = %6.4f\n",i,zf);
    sp = del_pos(sp);
  }

  params = del_lp(params);

  return 0;
}

int MAIN_(void)
{
  return 0;
}
