#include <stdlib.h>
#include "contact.h"
#include "algo.h"
#include "math.h"

void gen_velocities (iREAL lo[3], iREAL hi[3], unsigned int nt, iREAL * v[3]);
void integrate (iREAL step, iREAL lo[3], iREAL hi[3], unsigned int nt, iREAL * t[3][3], iREAL * v[3]);

void dynamics (master_conpnt master[], slave_conpnt slave[],
  int nt, iREAL * angular[6], iREAL * v[3],
  iREAL * rotation[9], iREAL * position[6],
  iREAL * inertia[9], iREAL * inverse[9],
  iREAL mass[], iREAL invm[], iREAL * force[3],
  iREAL * torque[3], iREAL gravity[3], iREAL step);


static void expmap (iREAL Omega1, iREAL Omega2, iREAL Omega3,
                iREAL &Lambda1, iREAL &Lambda2, iREAL &Lambda3,
			          iREAL &Lambda4, iREAL &Lambda5, iREAL &Lambda6,
			          iREAL &Lambda7, iREAL &Lambda8, iREAL &Lambda9);
