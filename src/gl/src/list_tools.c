#include<stdio.h>
#include<stdlib.h>
#include "structdef.h"
#include "list_tools.h"

int init_list(List *list)
{
  *list = NULL;
  return 0;
}

int allocate_node(List *list, Posinfo x, double del)
{
  List L;

  L = (List) malloc(sizeof(Node));
  if (!L)
    return 1;

  *list = L;
  L->impos = x;
  L->delb = del;
  L->next = NULL;
  return 0;
}

void free_node(List *list)
{
  free(*list);
  *list=NULL;
}

void free_list(List *list)
{
  struct Node *node;  /* The node being freed */
  struct Node *next;  /* The next node to be freed */
/*
 * Free each node of the list.
 */
  for(node=*list; node; node=next) {
    next=node->next;
    free_node(&node);
  }
  *list=NULL;
}

int empty_list(List L)
{
  return (L == NULL) ? 1 : 0;
}

int append_list(List *list, Posinfo posinfo, double del)
{
  List L,tmplist;

  if(allocate_node(&L,posinfo,del))
    return 1;
  if(empty_list(*list))
    *list = L;
  else {
    for(tmplist=*list;tmplist->next; tmplist=tmplist->next);
    tmplist->next = L;
  }

  return 0;
}
  
