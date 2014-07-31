#include<stdio.h>

/* detrate.c
 *
 * Usage: detrate
 *
 * Description: Determines the number of CLASS detections in different
 *                regions of parameter space, separated by the following
 *                criteria:
 *                 1) Texas-87GB spectral index = -0.6
 *                 2) W&B-87GB spectral index = -0.5
 *                 3) Declination = furthest northern range of Texas survey
 *
 * 05May95, CDF
 * v13May95, CDF  Modified Texas Declination limit and included detection
 *                 percentage.
 */

#define TEXLIM -0.6
#define WBLIM -0.5
#define DECLIM 71.6
#define MAXCHAR 1000

main(int argc, char *argv[])
{
  int count1,count2;
  int rahr,ramin;
  int decdeg;
  int decmin;
  int texyes,texno;
  int wbyes,wbno;
  int texwbyes,texwbno;
  int northyes,northno;
  int southyes,southno;
  float rasec,decsec;
  float atex,awb,fluxrat;
  char source[MAXCHAR];
  char sign1[MAXCHAR],sign2[MAXCHAR];
  char line[MAXCHAR];
  FILE *ifp1,*ifp2,*ofp;

  /* Open files */
  ifp1 = fopen("indices.list","r");
  ofp = fopen("outrate.list","w");

  /* initialize */
  texyes = 0;
  texno = 0;
  wbyes = 0;
  wbno = 0;
  texwbyes = 0;
  texwbno = 0;
  northyes = 0;
  northno = 0;
  southyes = 0;
  southno = 0;
  count1=0;

  while(fgets(line,MAXCHAR,ifp1) != NULL) {
    if (line[0] != '#') {
      count1++;
      sscanf(line,"%s %s %f %s %f %f",source,sign1,&atex,sign2,&awb,&fluxrat);
      if (fluxrat == 0.0) {
	if (atex > TEXLIM)
	  texno++;
	if (awb > WBLIM)
	  wbno++;
	if (atex > TEXLIM && awb > WBLIM)
	  texwbno++;
      }
      else {
	if (atex > TEXLIM)
	  texyes++;
	if (awb > WBLIM)
	  wbyes++;
	if (atex > TEXLIM && awb > WBLIM)
	  texwbyes++;
      }
    }
  }
  fclose(ifp1);

  fprintf(ofp,"Total Number of Sources:  %d\n\n",count1);
  fprintf(ofp,"tex-87gb > -0.6 | W&B-87gb > -0.5 | Detections | Nondetections | Detection rate\n");
  fprintf(ofp,"-------------------------------------------------------------------------------\n");
  fprintf(ofp,"       X        |        X        |  %4d      |   %4d        |    %4.1f%%\n",
	 texwbyes,texwbno,100.0*texwbyes/(texwbyes+texwbno));
  fprintf(ofp,"       X        |                 |  %4d      |   %4d        |    %4.1f%%\n",
	 texyes,texno,100.0*texyes/(texyes+texno));
  fprintf(ofp,"                |        X        |  %4d      |   %4d        |    %4.1f%%\n",
	 wbyes,wbno,100.0*wbyes/(wbyes+wbno));

  ifp1 = fopen("detections.list","r");
  while(fgets(line,MAXCHAR,ifp1) != NULL) {
    if(line[0] != '!') {
      sscanf(line,"%s %d:%d:%f %d:%d:%f",source,&rahr,&ramin,&rasec,
	     &decdeg,&decmin,&decsec);
      if ((1.0*decdeg+decmin/60.0) >= DECLIM)
	northyes++;
      else
	southyes++;
    }
  }
  fclose(ifp1);

  ifp1 = fopen("nondetections.list","r");
  while(fgets(line,MAXCHAR,ifp1) != NULL) {
    if(line[0] != '!') {
      sscanf(line,"%s %d:%d:%f %d:%d:%f",source,&rahr,&ramin,&rasec,
	     &decdeg,&decmin,&decsec);
      if ((1.0*decdeg+decmin/60.0) >= DECLIM)
	northno++;
      else
	southno++;
    }
  }
  fclose(ifp1);
  
  fprintf(ofp,"\nNorth of Texas range (dec > %4.1f):\n",DECLIM);
  fprintf(ofp,"  Detections:     %4d\n",northyes);
  fprintf(ofp,"  Nondetections:  %4d\n",northno);
  fprintf(ofp,"  Detection Rate: %4.2f%%\n\n",100.0*northyes/(northyes+northno));

  fprintf(ofp,"Inside of Texas range (dec < %4.1f):\n",DECLIM);
  fprintf(ofp,"  Detections:     %4d\n",southyes);
  fprintf(ofp,"  Nondetections:  %4d\n",southno);
  fprintf(ofp,"  Detection Rate: %4.2f%%\n\n",100.0*southyes/(southyes+southno));
  
  fclose(ofp);

  return 0;
}
