/*
 * specindx.c
 *
 * Usage: specindx
 *
 * Calculates a two-point spectral index.
 *
 * Original version: CDF, ancient
 * v29Sep99, CDF  Added looping possibility.
 *                Better error checking.
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAXC 1000

int main(int argc, char *argv[])
{
  int go_again=1;       /* Flag set to 0 to break loop */
  float nu1,nu2;        /* Frequencies */
  float f1,f2;          /* Flux densities */
  float alph;           /* Spectral index */
  char line[MAXC];      /* General string for reading input */

  /*
   * Get the frequencies of the two endpoints
   */

  printf("Enter nu1 and nu2:  ");
  fgets(line,MAXC,stdin);
  while(sscanf(line,"%f %f",&nu1,&nu2) != 2) {
    fprintf(stderr,"ERROR: Bad input.  Enter nu1 and nu2 again:  ");
    fgets(line,MAXC,stdin);
  }
  printf("\n****************************************\n");
  printf("nu1 = %f, flux1 = %f\n",nu1,f1);
  printf("****************************************\n\n");

  /*
   * Now loop through, getting flux densities and calculating
   *  spectral indices for as long as the user wants.
   */

  while(go_again) {
    printf("Enter flux1 and flux2:  ");
    fgets(line,MAXC,stdin);
    while(sscanf(line,"%f %f",&f1,&f2) != 2) {
      fprintf(stderr,"ERROR: Bad input.  Enter flux1 and flux2 again:  ");
      fgets(line,MAXC,stdin);
    }

    alph = (log10(f2) - log10(f1))/(log10(nu2) -  log10(nu1));

    printf("\nnu2 = %f, flux2 = %f\n",nu2,f2);
    printf("Spectral index = %6.3f\n",alph);

    /*
     * Now see if further calculations are desired.
     */

    printf("\nCalculate another spectral index for the same pair of ");
    printf("frequencies? [y] ");
    fgets(line,MAXC,stdin);
    if(line[0] == 'n' || line[0] == 'N')
      go_again = 0;

    printf("\n****************************************\n\n");
  }

  return 0;
}
