#include<stdio.h>
#include<math.h>

/* hoggles.c
 *
 * Usage: hoggles
 *
 * Description: Program to find lenses in HST fields
 *
 * 13Sep95, CDF
 */

#define DR 0.5
#define MAXSIZE 0.5
#define MAXLINE 1000
#define MAXARR 2000
#define ELLIPMIN 0.5
#define PI 3.141592653589793
#define DATA(L) ((L)->datapointer)
#define LISTNEXT(L) ((L)->listnext)
#define NEXT(L) ((L)->next)
#define PREV(L) ((L)->prev)
#define WFPIXSC 0.1

typedef enum {OK, ERROR } status;
typedef enum { FALSE=0, TRUE=1} bool;

struct fileline {
  int sourcenum;
  float x;
  float y;
  float iap;
  float iiso;
  float sigap;
  float sigiso;
  float a;
  float b;
  float theta;
  int groupnum;
  int groupflag;
  int doneflag;
};
 
typedef struct node node, *list;
struct node {
  struct fileline datapointer;
  list listnext;
  list next;
  list prev;
};

status init_list(list *p_L);
status allocate_node(list *p_L, struct fileline data);
void free_node(list *p_L);
bool empty_list(list L);
status append(list *p_L, struct fileline data);
list find_first(list *p_L);
list find_last(list *p_L);
void print_list(list *p_L);

main(int argc, char *argv[])
{
  int i,j;
  int count,compact,group;
  int totnum;
  int gpmatx[100,100];
  float r2,delx,dely;
  char filename[MAXLINE],line[MAXLINE];
  char ellipoutname[MAXLINE];
  struct fileline info,*p_start,*p_i,*p_j,*p_end;
  list filelist, tmpptr,inlist,holdpt;
  FILE *ifp,*eop,*ofp;

  init_list(&filelist);

  printf("Enter the name of the input file:  ");
  gets(filename);
  while((ifp = fopen(filename, "r")) == NULL) {
    printf("  %s is not a valid file, please enter a new filename.\n",filename);
    printf("       ");
    gets(filename);
  }
  printf("Enter the name of the output file for elliptical sources:  ");
  gets(ellipoutname);
  while((eop = fopen(ellipoutname, "w")) == NULL) {
    printf("  %s is not opening.  Please try again");
    printf("       ");
    gets(ellipoutname);
  }

  /* Read in info and put it into pointer array */
  totnum = 0;
  while(fgets(line,MAXLINE,ifp) != NULL)
    if(line[0] != '!')
      totnum++;
  fclose(ifp);

  printf("%d sources in input file.\n\n",totnum);

  ifp = fopen(filename,"r");
  i = 0;
  while(fgets(line,MAXLINE,ifp) != NULL) {
    if(line[0] != '!') {
      i++;
      sscanf(line,"%f %f %f %f %f %f %f %f %f",
	     &info.x,&info.y,&info.iap,&info.iiso,&info.sigap,&info.sigiso,
	     &info.a,&info.b,&info.theta);
      info.sourcenum = i;
      info.groupnum = 0;
      info.doneflag = 0;
      append(&filelist,info);
    }
  }
  fclose(ifp);


  /*
   * if find_first(list1) != find_first(list2) &&
   *  r*r < DR*DR
   *  then find_last(list1) -> find_first(list2)
   *
   */
  

  for(tmpptr=filelist; tmpptr; tmpptr=LISTNEXT(tmpptr)) {
    if(DATA(tmpptr).a > (4 * DATA(tmpptr).b))
      fprintf(eop,"%3d. %7.2f %7.2f %8.4f %8.4f\n",DATA(tmpptr).sourcenum,
	      DATA(tmpptr).x,DATA(tmpptr).y,
	      WFPIXSC*DATA(tmpptr).a,WFPIXSC*DATA(tmpptr).b);
    for(inlist=LISTNEXT(tmpptr); inlist; inlist=LISTNEXT(inlist))
      if(find_first(&tmpptr) != find_first(&inlist)) {
	delx = WFPIXSC * (DATA(tmpptr).x - DATA(inlist).x);
	dely = WFPIXSC * (DATA(tmpptr).y - DATA(inlist).y);
	if((delx*delx + dely*dely) < DR*DR) {
	  holdpt = find_last(&tmpptr);
	  NEXT(holdpt) = find_first(&inlist);
	  PREV(find_first(&inlist)) = holdpt;
	  printf("%2d paired with %2d",DATA(tmpptr).sourcenum,
		 DATA(inlist).sourcenum);
	  printf("  separated by %5.2f arcsec\n",sqrt(delx*delx + dely*dely));
	}
      }
  }
  printf("\n");
  fclose(eop);

  count=0;
  for(tmpptr=filelist; tmpptr; tmpptr=LISTNEXT(tmpptr)) {
    if(NEXT(find_first(&tmpptr)) && !DATA(find_first(&tmpptr)).doneflag) {
      count++;
      compact=0;
      printf("Group #%d\n",count);
      print_list(&tmpptr);
      DATA(find_first(&tmpptr)).doneflag = 1;
      printf("\n");
    }
  }

  return 0;

}

status init_list(list *p_L)
{
  *p_L = NULL;
  return OK;
}

status allocate_node(list *p_L, struct fileline data)
{
  list L = (list) malloc(sizeof(node));

  if (L == NULL)
    return ERROR;

  *p_L = L;
  DATA(L) = data;
  LISTNEXT(L) = NULL;
  NEXT(L) = NULL;
  PREV(L) = NULL;
  return OK;
}

void free_node(list *p_L)
{
  free(*p_L);
  *p_L = NULL;
}

bool empty_list(list L)
{
  return (L == NULL) ? TRUE : FALSE;
}

status append(list *p_L, struct fileline data)
{
  list L, tmplist;

  if(allocate_node(&L, data) == ERROR)
    return ERROR;

  if (empty_list(*p_L) == TRUE)
    *p_L = L;
  else {
    for (tmplist = *p_L; LISTNEXT(tmplist)!= NULL; tmplist=LISTNEXT(tmplist));
    LISTNEXT(tmplist) = L;
  }
  return OK;

}

list find_first(list *p_L)
{
  list tmplist;

  for (tmplist = *p_L; PREV(tmplist)!=NULL; tmplist=PREV(tmplist));
  return tmplist;
}

list find_last(list *p_L)
{
  list tmplist;

  for (tmplist = *p_L; NEXT(tmplist)!=NULL; tmplist=NEXT(tmplist));
  return tmplist;
}

void print_list(list *p_L)
{
  list tmplist;
  
  for (tmplist = find_first(p_L); tmplist; tmplist=NEXT(tmplist)) {
    printf(" %3d.  ",DATA(tmplist).sourcenum);
    printf("%6.2f %6.2f\n",DATA(tmplist).x,DATA(tmplist).y);
  }
}
    
