#ifndef HEADER_FILE
#define HEADER_FILE

#include "algo.h"
#include "bf.h"

#define CONBUF 8

/* master contact points; global array is used because
 * constitutive data at contacts can be persistent */
struct master_conpnt
{
  int master[CONBUF];  //element
  int slave[2][CONBUF]; // particle, element
  int color[2][CONBUF];
  iREAL point[3][CONBUF];
  iREAL normal[3][CONBUF];
  iREAL depth[CONBUF];
  iREAL force[3][CONBUF];
  int size;

  struct master_conpnt * next; /* local list */
};

/* slave contact points; they are created by
 * symmetrically coppying master contact points */
struct slave_conpnt
{
  int master[2][CONBUF]; /* particle, triangle */
  iREAL point[3][CONBUF];
  iREAL force[3][CONBUF];
  int size;

  struct slave_conpnt * next; /* local list */
};

/* calculate distances */
void contact_detection (int s1, int e1, int s2, int e2, 
                        iREAL *t[6][3], int *tid, int *pid, iREAL *v[3], 
                        iREAL *p[3], iREAL *q[3], master_conpnt *con);

master_conpnt * newcon (master_conpnt * master, int *k);
slave_conpnt * newcon (slave_conpnt * slave, int *k);

void master_free (master_conpnt * con);
void slave_free (slave_conpnt * con);

#endif
