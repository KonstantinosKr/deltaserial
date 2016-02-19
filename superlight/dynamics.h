
#include <stdlib.h>
#include "algo.h"
#include "math.h"
#include "contact.h"
#include "material.h"


void dynamics (master_conpnt master[], slave_conpnt slave[],
  int nt, int nb, iREAL *t[6][3], int *pid, iREAL * linear[3], iREAL * position[6],
  iREAL mass[], iREAL * force[3], iREAL gravity[3], iREAL step);

