#include<stdio.h>
#include<math.h>
#include "cpgplot.h"
#include "structdef.h"

#define MAXLINE 1000
#define c 3.0e10/3.1e24  /* To convert to Mpc/sec */
#define H 3.23e-18       /* Because H=100 h km/sec/Mpc = 3.2e-18 h 1/sec */

void set_omega(float *omega, int npoints);
float *calc_f(float zl, float zs, float *omega, int npoints);
int plot_graph(int npoints,float *omega,float *f);

int main(int argc, char *argv[])
{
  int npoints;
  float zl,zs;
  float *omega;
  float *f;
  char line[MAXLINE];

  npoints = 1000;
  if(!(omega = new_array(npoints,1)))
    return 1;

  set_omega(omega,npoints);

  printf("Enter the redshift of the lens:  ");
  gets(line);
  sscanf(line,"%f",&zl);
  printf("Enter the redshift of the source:  ");
  gets(line);
  sscanf(line,"%f",&zs);

  if(!(f = calc_f(zl,zs,omega,npoints)))
    return 1;

  if(plot_graph(npoints,omega,f) == 1)
    return 1;

  omega = del_array(omega);
  f = del_array(f);

  return 0;
}

void set_omega(float *omega, int npoints)
{
  int i;
  float *optr;

  optr = omega;
  for(i=0;i<npoints;i++,optr++)
    *optr = (0.1+i*0.9/(npoints-1));
}

float *calc_f(float zl, float zs, float *omega, int npoints)
{
  int i;
  float gl,gs;
  float dl,ds,dls;
  float *f_array,*fptr,*optr;

  if(!(f_array = new_array(npoints,1)))
    return NULL;

  dl = 2*c*((1+zl) - sqrt(1+zl))/(H*(1+zl)*(1+zl));
  ds = 2*c*((1+zs) - sqrt(1+zs))/(H*(1+zs)*(1+zs));
  dls = 2*c*sqrt(1+zl)*sqrt(1+zs)*(sqrt(1+zs)-sqrt(1+zl))/(H*(1+zl)*(1+zs)*(1+zs));
  printf("\nFor this lens D_l=%4.0f, D_s=%4.0f and D_ls=%4.0f h^{-1} Mpc\n",
	 dl,ds,dls);
  printf(" assuming Omega_0 = 1.0\n\n");
  optr = omega;
  fptr = f_array;
  for(i=0;i<npoints;i++,optr++,fptr++) {
    gl = sqrt(1 + *optr * zl);
    gs = sqrt(1 + *optr * zs);
    *fptr = (1 - *optr - gl)*(1 - *optr - gs)*(1-gl)*(1-gs) /
      (*optr * *optr * (1 - *optr - gl*gs) * (gl - gs));
  }

  return f_array;
}

int plot_graph(int npoints,float *omega,float *f)
{
  int i;
  float yhi,ylo,f1;
  float *fnorm,*yptr,*fptr;

  if(!(fnorm = new_array(npoints,1)))
    return 1;

  fptr = f;
  f1 = *(fptr+npoints-1);
  yptr = fnorm;
  yhi = ylo = *yptr;
  for(i = 0;i<npoints;i++,yptr++,fptr++) {
    *yptr = *fptr / f1;
    if(yhi < *yptr)
      yhi = *yptr;
    if(ylo > *yptr)
      ylo = *yptr;
  }

  if(cpgbeg(0, "?", 1, 1) != 1) {
    fprintf(stderr,"Error opening device.\n");
    return 1;
  }
  cpgenv(0.0,1.0,ylo,yhi,0,0);
  cpgline(npoints,omega,fnorm);
  cpglab("\\gW","\\gDt(\\gW)/\\gDt(1.0)",
	 "Effect of \\gW on Time Delay");
  cpgend();
  
  return 0;
}

int MAIN_(void)
{
  return 0;
}
