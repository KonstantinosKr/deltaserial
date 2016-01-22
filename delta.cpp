#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include "error.h"
#include "tmr.h"
#include "input.h"
#include "output.h"
#include "contact.h"
#include "dynamics.h"
#include "forces.h"

int main (int argc, char **argv)
{
  iREAL *t[3][3]; /* triangles */
  iREAL *v[3]; /* velocities */
  iREAL *mass; /* scalar mass */
  iREAL *force[3]; /* total spatial force */
  iREAL *torque[3]; /* total spatial torque */
 
  int *parmat; /* particle material */
  iREAL *mparam[NMAT]; /* material parameters */ 
  for(int i = 0; i<NMAT; i++)//number of material parameters
  {
    mparam[i] = (iREAL *) malloc(1*sizeof(iREAL));//n number of materials
  }

  int pairnum = 1;
  int *pairs; /* color pairs */
  pairs = (int *) malloc(pairnum*sizeof(pairnum));

  int *ikind; /* interaction kind */
  ikind = (int *) malloc(1*sizeof(int)); //number of interaction kinds/types

  iREAL *iparam[NINT]; // interaction parameters // parameters per interaction type
  for(int i=0;i<NINT;i++)
  {
    iparam[i] = (iREAL *) malloc(1*sizeof(iREAL));
  }


  //set first kind
  ikind[0] = GRANULAR;
  
  //GRANULAR interaction type parameters 
  iparam[SPRING][GRANULAR] = 0;
  iparam[DAMPER][GRANULAR] = 0;
  iparam[FRISTAT][GRANULAR] = 0;
  iparam[FRIDYN][GRANULAR] = 0;
  iparam[FRIROL][GRANULAR] = 0;
  iparam[FRIDRIL][GRANULAR] = 0;
  iparam[KSKN][GRANULAR] = 0;
  iparam[LAMBDA][GRANULAR] = 0;
  iparam[YOUNG2][GRANULAR] = 0;
  iparam[KSKN2][GRANULAR] = 0;
  iparam[SIGC][GRANULAR] = 0;
  iparam[TAUC][GRANULAR] = 0;
  iparam[ALPHA][GRANULAR] = 0;

  iREAL *angular[6]; /* angular velocities (referential, spatial) */
  //iREAL *linear[3]; /* linear velocities */
  iREAL *rotation[9]; /* rotation operators */
  iREAL *position[6]; /* mass center current and reference positions */
  iREAL *inertia[9]; /* inertia tensors */
  iREAL *inverse[9]; /* inverse inertia tensors */
  iREAL *invm; /* inverse scalar mass */
 
  iREAL gravity[3];
  gravity[0] = 0;
  gravity[1] = 0;
  gravity[2] = 0;
  iREAL *distance; /*distance */
  iREAL *p[3],*q[3];//p and q points
  
  unsigned int nt = 0; /* number of triangles */
  unsigned int *pid; /*particle identifier */
  unsigned int *tid; /* triangle identifiers */
  master_conpnt *con = 0; slave_conpnt *slave = 0;
  iREAL lo[3] = {-500, -500, -500}; /* lower corner */
  iREAL hi[3] = {500, 500, 500}; /* upper corner */
  
  unsigned int size = 27000000; /* memory buffer size */
   
	for (int i = 0; i < 3; i ++)
	{ 
		t[0][i] = (iREAL *) malloc (size*sizeof(iREAL));
		t[1][i] = (iREAL *) malloc (size*sizeof(iREAL));
		t[2][i] = (iREAL *) malloc (size*sizeof(iREAL));
		v[i] = (iREAL *) malloc (size*sizeof(iREAL));
    
    position[i] = (iREAL *) malloc (size*sizeof(iREAL)); 
    torque[i] = (iREAL *) malloc (size*sizeof(iREAL));
    force[i] = (iREAL *) malloc (size*sizeof(iREAL));

		p[i] = (iREAL *) malloc (size*sizeof(iREAL));
		q[i] = (iREAL *) malloc (size*sizeof(iREAL));
  }
	
	for (int i = 0; i < 6; i++)
	{
		angular[i] = (iREAL *) malloc (size*sizeof(iREAL));
	}
 
  for (int i = 0; i<9; i++)
  {
    inverse[i] = (iREAL *) malloc (size*sizeof(iREAL));
    inertia[i] = (iREAL *) malloc (size*sizeof(iREAL));
    rotation[i] = (iREAL *) malloc (size*sizeof(iREAL));
  }

  con = (master_conpnt *) malloc (size*sizeof(master_conpnt));
  slave = (slave_conpnt *) malloc (size*sizeof(slave_conpnt));
	
  parmat = (int *) malloc (size*sizeof(int));
	
	tid = (unsigned int *) malloc (size*sizeof(unsigned int));
	pid = (unsigned int *) malloc (size*sizeof(unsigned int));

  invm = (iREAL *) malloc(size*sizeof(iREAL));
	mass = (iREAL *) malloc(size*sizeof(iREAL));

	for(unsigned int i=0;i<size;i++) tid[i] = UINT_MAX; 
	
	unsigned int nb;
	init_enviroment(&nt, &nb, t, v, angular, inertia, inverse, mass, invm, parmat, tid, pid, position, lo, hi);  
	printf("NT:%i NB:%i\n", nt, nb);
  
  unsigned long long int ncontacts = 0;
  
  /* perform time stepping */
  iREAL step = 1E-3, time; unsigned int timesteps=0; 
  
  for(time = 0; time < 0.1; time+=step)
  {
    printf("TIMESTEP: %i\n", timesteps); 
   
    contact_detection (0, nt, 0, nt, t, tid, pid, v, step, p, q, con, &ncontacts);
		 
    forces(con, slave, nt, position, angular, v, mass, invm, parmat, mparam, pairnum, pairs, ikind, iparam);
  
    dynamics(con, slave, nt, angular, v, rotation, position, inertia, inverse, mass, invm, force, torque, gravity, step);
    
    output_state(nt, t, v, timesteps);
    
    timesteps++;
  }
	printf("\nComputation Finished.\n");

  for (int i = 0; i < 3; i ++)
  {
    free (t[0][i]);
    free (t[1][i]);
    free (t[2][i]);
    free (v[i]);
    free (p[i]);
    free (q[i]);
  }

  return 0;
}
