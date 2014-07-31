/* compinfo.c                                                   
 *                                                              
 * Program to take as input the output file from pos.c and      
 *   report whether there is only one component in the source   
 *   or, if there are more than one, the flux ratio between the 
 *   brightest and next brightest components.                   
 *                                                              
 * Usage: compinfo [input file name]                            
 *                                                              
 * 01Nov94, CDF  A modification of ptfind.c                      
 * v22Mar95, CDF Fixed bug that didn't update next brightest    
 *               component flux                                 
 * v27Mar95, CDF Cleaned up 
 * v26Sep98, CDF Slight stylistic changes.                                    
 */


#include <stdio.h>
#include <string.h>

#define MAXLINE 1000

int main(int argc, char *argv[])
{
  int i;
  int compno,compno2;    /* Component number */
  float flux[25];
  float flux2;
  float highflux,secflux;
  float fratio;
  float maj,ratio,phi;
  float maj2,ratio2;
  float minor;
  char source[MAXLINE];  /* Source name */
  char source2[MAXLINE];
  char hr[MAXLINE],min[MAXLINE],sec[MAXLINE];
  char deg[MAXLINE],amin[MAXLINE],asec[MAXLINE];
  char hr2[MAXLINE],min2[MAXLINE],sec2[MAXLINE];
  char deg2[MAXLINE],amin2[MAXLINE],asec2[MAXLINE];
  char line[MAXLINE];
  char outname[MAXLINE];
  FILE *ifp;              /* Input file pointer */
  FILE *ofp;

  /* Get input */
  compno = 0;
  ifp = fopen(argv[1], "r" );
  printf("\n Enter the name of the output file:  ");
  gets(outname);
  ofp = fopen(outname, "w" );
  fputs("Source   Flux (mJy)    RA (2000)     Dec (2000)    Size(mas)\n",ofp);
  fputs("-------  ---------- -------------  -------------  ---------\n",ofp);

  while (fgets(line, MAXLINE, ifp) != NULL) {
      maj2 = ratio = 0;
      sscanf(line,"%s %d %f %s %s %s %s %s %s %f %f",source2,&compno2,&flux2,
	     hr2,min2,sec2,deg2,amin2,asec2,&maj2,&ratio2);
      if(compno2 == 1) {
	  /* Previous source has only one component */
	  if(compno == 1) {
	      sprintf(line,
	             "%s %6.1f     %s %s %s  %s %s %s %5.0f x %-5.0f (one component)\n",
		     source,highflux,hr,min,sec,deg,amin,asec,maj,minor);
	      fputs(line,ofp);
	  }

	  else if(compno > 1) {
	      secflux = 0;
	      for(i = 0; i < compno; i++)
		if(flux[i] != highflux && flux[i] > secflux) {
		  fratio = highflux / flux[i];
		  secflux = flux[i];
		}
	      fprintf(ofp,
		      "%s %6.1f     %s %s %s  %s %s %s %5.0f x %-5.0f (%3.0f:1 double)\n",
		      source,highflux,hr,min,sec,deg,amin,asec,maj,minor,fratio);
	    }

	  highflux = flux2;
	  strcpy(hr,hr2);
	  strcpy(min,min2);
	  strcpy(sec,sec2);
	  strcpy(deg,deg2);
	  strcpy(amin,amin2);
	  strcpy(asec,asec2);
	  maj = maj2;
	  minor = maj2 * ratio2;
	  strcpy(source,source2);
      }

      /* If this isn't the first component, check the component flux */
      /*   against the previous component fluxes for this source. */
      else if(flux2 > highflux) {
	  highflux = flux2;
	  strcpy(hr,hr2);
	  strcpy(min,min2);
	  strcpy(sec,sec2);
	  strcpy(deg,deg2);
	  strcpy(amin,amin2);
	  strcpy(asec,asec2);
	  maj = maj2;
	  minor = maj2 * ratio2;
      }

      flux[compno2-1] = flux2;
      compno = compno2;
  }

  if(compno == 1) {
      fprintf(ofp,
	     "%s %6.1f     %s %s %s  %s %s %s %5.0f x %-5.0f (one component)\n",
	     source,highflux,hr,min,sec,deg,amin,asec,maj,minor);
  }
  else if(compno > 1) {
      secflux = 0;
      for(i = 0; i < compno; i++)
	  if(flux[i] != highflux && flux[i] > secflux) {
	    fratio = highflux / flux[i];
	    secflux = flux[i];
	  }
      fprintf(ofp,
	      "%s %6.1f     %s %s %s  %s %s %s %5.0f x %-5.0f (%3.0f:1 double)\n",
	      source,highflux,hr,min,sec,deg,amin,asec,maj,minor,fratio);
  }
      

  fclose(ifp);
  fclose(ofp);

  return 0;
}

