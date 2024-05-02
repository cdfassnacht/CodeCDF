#include<stdio.h>
#include<math.h>

#define MAXLINE 1000
#define STEP 10000
#define PI 3.141592653589793
#define c 3.0e10/3.1e24      /* To convert to Mpc/sec */
#define H 3.23e-18           /* Because H=100 h km/sec/Mpc = 3.2e-18 h 1/sec */

int main(int argc, char *argv[])
{
  int i;
  float omega;
  float n0;
  float N;
  float dz;
  float z;
  float zmax;
  float f;
  float sum;
  float tau;
  char line[MAXLINE];

  printf("Enter z_max  ");
  gets(line);
  sscanf(line,"%f",&zmax);
  printf("You have entered z_max = %f\n",zmax);
  printf("Enter Omega_0:  ");
  gets(line);
  sscanf(line,"%f",&omega);
  printf("Enter n_0 in units of 1/Mpc^3:  ");
  gets(line);
  sscanf(line,"%f",&n0);

  dz = zmax / STEP;
  printf("dz = %f\n",dz);
  z = dz/2;
  sum = 0.0;
  for(i=0;i<STEP;i++) {
    z += dz;
    f = (omega*z+(omega-2)*(sqrt(1+omega*z)-1))
      * (omega*z+(omega-2)*(sqrt(1+omega*z)-1))
	/ ((1+z)*(1+z)*(1+z)*sqrt(1+omega*z));
    sum += f*dz;
  }

  N = 16*PI*n0*pow((c/H),3.0)*sum/pow(omega,4.0);

  /* Calculate the optical depth assuming all lenses split by 1 arcsec */
  tau = 10e-10*PI*n0*pow((c/H),3.0)*sum/pow(omega,4.0);

  printf("N = %g, tau=%g\n",N,tau);

  return 0;
}
