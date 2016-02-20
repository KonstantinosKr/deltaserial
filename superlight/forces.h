#include "math.h"
#include "material.h"
#include "stdio.h"
#include "contact.h"

void granular(iREAL n[3], iREAL vij[3], iREAL depth, int i, int j, iREAL mass[], iREAL *iparam[NINT], int ij, iREAL f[3]);

void forces (master_conpnt master[], slave_conpnt slave[], std::vector<contact> conpnt[],
  int nb, iREAL * position[3], iREAL * linear[3],
  iREAL mass[], int parmat[], iREAL * mparam[NMAT],
  int pairnum, int pairs[], int ikind[], iREAL * iparam[NINT]);

