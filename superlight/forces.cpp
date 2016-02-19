#include "forces.h"

void granular(iREAL n[3], iREAL vij[3], iREAL depth, int i, int j, iREAL mass[], iREAL *iparam[NINT], int ij, iREAL f[3])
{
  iREAL ma = 1.0 / (1/mass[i] + 1/mass[j]);
  
  iREAL kn = iparam[SPRING][ij];
  iREAL en = iparam[DAMPER][ij] * 2.0 * sqrt(kn*ma);
  iREAL vn = DOT(vij,n);
  iREAL fn = kn*depth + en*vn;
  //printf("kn:%f, en:%f, vn:%f, fn:%f depth:%f, vij[0]:%f vij[1]:%f vij[2]:%f\n", kn, en, vn, fn, depth, vij[0], vij[1], vij[2]); 
  f[0] = fn*n[0];
  f[1] = fn*n[1];
  f[2] = fn*n[2];
  printf("CONTACT F[0]: %f, F[1]: %f, F[2]: %f\n", f[0], f[1], f[2]); 
}

/*
 void granular(iREAL normal[3], iREAL relativeVelocity[3], iREAL depth, int i, int j, iREAL mass[], iREAL force[3])
 {//contact normal has to be unit vector
 
    iREAL ma = 1.0 / ((1/mass[i]) + (1/mass[j]);
    iREAL SPRING = 1E6; //dashpot-spring parameter
    iREAL DAMPER = 1E2 * sqrt(kn*ma); //damping parameter
    iREAL vn = DOT(relativeVelocity,normal); //maybe not needed for linear contact
    iREAL fn = kn*depth + en*vn; //force magnitude
 
    force[0] = fn*normal[0];
    force[1] = fn*normal[1];
    force[2] = fn*normal[2];
 }
 */

/* return pairing index based on (i,j) pairing of colors */
int pairing (int nummat, int pairs[], int i, int j)
{
  return 0; /* default material */
}

void forces (master_conpnt master[], slave_conpnt slave[],
    int nb, iREAL * position[3], iREAL * linear[3],
    iREAL mass[], int parmat[], iREAL * mparam[NMAT],
    int pairnum, int pairs[], int ikind[], iREAL * iparam[NINT])
{
  for (int i = 0; i < nb; i++)
  {
    iREAL v[3], x[3];

    v[0] = linear[0][i];
    v[1] = linear[1][i];
    v[2] = linear[2][i];

    x[0] = position[0][i];
    x[1] = position[1][i];
    x[2] = position[2][i];

    /* update contact forces */
    for (master_conpnt * con = &master[i]; con; con = con->next)
    {
      for(int k = 0; k<con->size; k++)
      {
        iREAL p[3], n[3], vi[3], vj[3], vij[3];

        p[0] = con->point[0][k];
        p[1] = con->point[1][k];
        p[2] = con->point[2][k];

        n[0] = con->normal[0][k];
        n[1] = con->normal[1][k];
        n[2] = con->normal[2][k];

        vi[0] = v[0];
        vi[1] = v[1];
        vi[2] = v[2];

        int j = con->slave[0][k]; //get index from slave body contact

        vj[0] = linear[0][j];
        vj[1] = linear[1][j];
        vj[2] = linear[2][j];

        SUB (vj, vi, vij); // relative linear velocity

        int ij = pairing (pairnum, pairs, con->color[0][k], con->color[1][k]);//get material from colours
      
        iREAL f[3];

        switch (ikind[ij])
        {
          case GRANULAR:
            granular (n, vij, con->depth[k], i, j, mass, iparam, ij, f);
            break;
          default:
            printf ("ERROR: invalid pairing kind");
            break;
        }
       
        con->force[0][k] = f[0];
        con->force[1][k] = f[1];
        con->force[2][k] = f[2];
      }
    }
    
    /* symmetrical copy into slave contact points */
    for (master_conpnt * con = &master[i]; con; con = con->next)
    {
      for (int j = 0; j < con->size; j++)
      {
        slave_conpnt *ptr;
        int k;

        ptr = newcon (&slave[con->slave[0][j]], &k);

        ptr->master[0][k] = i;
        ptr->master[1][k] = con->master[j];
        
        ptr->point[0][k] = con->point[0][j];
        ptr->point[1][k] = con->point[1][j];
        ptr->point[2][k] = con->point[2][j];
        
        ptr->force[0][k] = -con->force[0][j];
        ptr->force[1][k] = -con->force[1][j];
        ptr->force[2][k] = -con->force[2][j];
      }
    }
  }
}


