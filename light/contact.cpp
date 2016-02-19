#include "contact.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include <omp.h>

//using namespace ispc;

/* allocate new master contact point that can be written to */
master_conpnt * newcon (master_conpnt * master, int *k)
{
  master_conpnt * con = master; 
  while (con->size == CONBUF && con->next != NULL) con = con->next; // find available item or rewind to end
   
  if (con->size < CONBUF)
  {
    *k = con->size++;
  }
  else
  {
    master_conpnt * ptr = new master_conpnt;
    ptr->size = 0;
    ptr->next = NULL;
    con->next = ptr; // append new item at the end 
    con = ptr; // return new item 
    *k = 0;
  }
  return con;
}

/* allocate new slave contact point that can be written to */
slave_conpnt * newcon (slave_conpnt * slave, int *k)
{
  printf("Slave alloc");
  slave_conpnt * con = slave;
 
  printf("Slave alloc");
  while (con->size == CONBUF) con = con->next; /* rewind to the end */
  printf("Slave alloc");
  *k = con->size++; 
  if (con->size == CONBUF)
  {
    slave_conpnt * ptr = new slave_conpnt;
    ptr->size = 0;
    ptr->next = NULL;
    con->next = ptr; /* append new item at the end */
  } 
  
  con->size = 0;
  return con;
}

/* free global array of master contact points */
void master_free (master_conpnt * con)
{
  master_conpnt * ptr = con->next;
  while (ptr)
  {
    master_conpnt * next = ptr->next;
    delete ptr;
    ptr = next;
  }
  con->size = 0;
  con->next = NULL;
}

/* free global array of slave contact points */
void slave_free (slave_conpnt * con)
{
  slave_conpnt * ptr = con->next;
  while (ptr)
  {
    slave_conpnt * next = ptr->next;
    delete ptr;
    ptr = next;
  }
  con->next = NULL;
}

//s1 and e1 mean start of section 1 and end of section 1, same for s2,e2 and nt size nts1, nts2
void contact_detection (int s1, int e1, int s2, int e2, 
                        iREAL *t[6][3], int *tid, int *pid, iREAL *v[3], 
                        iREAL *p[3], iREAL *q[3], master_conpnt *con)
{
  iREAL a[3], b[3], c[3];

  //Set triangle 1 points A,B,C
  for(int i=s1;i<e1;i++)
  { 
    a[0] = t[0][0][i];
    a[1] = t[0][1][i];
    a[2] = t[0][2][i];
    
    b[0] = t[1][0][i];
    b[1] = t[1][1][i];
    b[2] = t[1][2][i];
    
    c[0] = t[2][0][i];
    c[1] = t[2][1][i];
    c[2] = t[2][2][i];
    
    bf (s2, e2, a, b, c, t[0], t[1], t[2], p, q);
     
    iREAL margin = 10E-2;
    
    for(int j=s2;j<e2;j++) //careful; range can overflow due to ghosts particles
    {
      iREAL dist = sqrt(pow((q[0][j]-p[0][j]),2)+pow((q[1][j]-p[1][j]),2)+pow((q[2][j]-p[2][j]),2));
      
      if(dist < margin  && (pid[i] != pid[j]))
      { //contact found 
        iREAL midpt[3], normal[3];
        
        midpt[0] = (p[0][j]+q[0][j])/2; //x
        midpt[1] = (p[1][j]+q[1][j])/2; //y
        midpt[2] = (p[2][j]+q[2][j])/2; //z
    
        iREAL depth = margin-dist;
        
        //iREAL mul = 1/sqrt(pow(midpt[0],2)+pow(midpt[1],2)+pow(midpt[2],2));
        
        iREAL mul = sqrt(pow(q[0][j]-p[0][j], 2)+pow(q[1][j] - p[1][j], 2)+pow(q[2][j]-p[2][j], 2));
        normal[0] = ((q[0][j] - p[0][j])/mul);//*depth for inclusion to normal
        normal[1] = ((q[1][j] - p[1][j])/mul);
        normal[2] = ((q[2][j] - p[2][j])/mul);
        
        master_conpnt * conpiv = &con[pid[i]];   

        int idx;
        conpiv = newcon (conpiv, &idx);

        conpiv->master[idx] = tid[i];
        conpiv->slave[0][idx] = pid[j];
        conpiv->slave[1][idx] = tid[j];

        //store contact point;
        conpiv->point[0][idx] = midpt[0];
        conpiv->point[1][idx] = midpt[1];
        conpiv->point[2][idx] = midpt[2];

        conpiv->normal[0][idx] = normal[0];
        conpiv->normal[1][idx] = normal[1];
        conpiv->normal[2][idx] = normal[2];

        conpiv->color[0][idx] = 0;
        conpiv->color[1][idx] = 0;

        conpiv->depth[idx] = depth;
      }
    }
  }
}

