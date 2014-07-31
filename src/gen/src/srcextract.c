#include <stdio.h>
#include <string.h>

#define MAX 1000

int main(int argc, char *argv[])
{
  int i;
  char inname1[MAX],inname2[MAX],outname[MAX];
  char line[MAX],line2[MAX];
  char source[MAX],source2[MAX];
  FILE *ifp,*ifp2,*ofp;

  /* Open files */
  printf("Enter the name of the file containing the source list:  ");
  gets(inname1);
  printf("Enter the name of the file containing the fluxes/indices\n");
  printf("  to be extracted:  ");
  gets(inname2);
  printf("Enter the name of the output file:  ");
  gets(outname);
  
  ifp = fopen(inname1,"r");
  ofp = fopen(outname,"w");

  i=0;
  while (fgets(line,MAX,ifp) != NULL) {
    sscanf(line,"%s",source);
    i++;
    printf("%d %s\n",i,source);
    ifp2 = fopen(inname2,"r");
    while(fgets(line2,MAX,ifp2) != NULL) {
      sscanf(line2,"%s",source2);
      if (strcmp(source,source2) == 0) {
	fputs(line2,ofp);
	break;
      }
    }
    fclose(ifp2);
  }
  
  fclose(ofp);
  fclose(ifp);

  return 0;
}
