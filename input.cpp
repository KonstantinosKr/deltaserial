#include "input.h" 
#include <float.h>
#include <algorithm>
#include <vector>
#include <set>
#include <map>
#include <iostream>
#include <fstream>
#include <iomanip>

// invert inertia properties 
void invert (int nb, iREAL * inertia[9], iREAL * inverse[9], iREAL mass[], iREAL invm[])
{
  for(int i = 0; i<nb;i++)
  {
    iREAL a[9], x[9], det;

    for (int j = 0; j < 9; j ++)
    {
      a[j] = inertia[j][i];
    }
    INVERT (a, x, det);

    for (int j = 0; j < 9; j ++)
    {
      inverse[j][i] = x[j];
    }

    invm[i] = 1.0 / mass[i];
  }
}

void getCentroid(int pid, unsigned int range1, unsigned int range2, iREAL *t[6][3], iREAL *centroid[6])
{
  iREAL cx=0;
  iREAL cy=0;
  iREAL cz=0;

  for(unsigned int i=range1;i<range2;i++)
  {
    cx = cx+t[0][0][i];
    cy = cy+t[0][1][i];
    cz = cz+t[0][2][i];
    
    cx = cx+t[1][0][i];
    cy = cy+t[1][1][i];
    cz = cz+t[1][2][i];
    
    cx = cx+t[2][0][i];
    cy = cy+t[2][1][i];
    cz = cz+t[2][2][i];
  }

  cx = cx/((range2-range1)*3);
  cy = cy/((range2-range1)*3);
  cz = cz/((range2-range1)*3);
    
  centroid[0][pid] = cx;
  centroid[1][pid] = cy;
  centroid[2][pid] = cz;
  
  centroid[3][pid] = cx;
  centroid[4][pid] = cy;
  centroid[5][pid] = cz;
}

void translate_enviroment(unsigned int i, unsigned int pid, iREAL *t[6][3], iREAL *p[6])
{
  t[0][0][i] = t[0][0][i] + p[0][pid];
  t[0][1][i] = t[0][1][i] + p[1][pid];
  t[0][2][i] = t[0][2][i] + p[2][pid];

  t[1][0][i] = t[1][0][i] + p[0][pid];
  t[1][1][i] = t[1][1][i] + p[1][pid];
  t[1][2][i] = t[1][2][i] + p[2][pid];

  t[2][0][i] = t[2][0][i] + p[0][pid];
  t[2][1][i] = t[2][1][i] + p[1][pid];
  t[2][2][i] = t[2][2][i] + p[2][pid];
  
  t[3][0][i] = t[3][0][i] + p[0][pid];
  t[3][1][i] = t[3][1][i] + p[1][pid];
  t[3][2][i] = t[3][2][i] + p[2][pid];

  t[4][0][i] = t[4][0][i] + p[0][pid];
  t[4][1][i] = t[4][1][i] + p[1][pid];
  t[4][2][i] = t[4][2][i] + p[2][pid];

  t[5][0][i] = t[5][0][i] + p[0][pid];
  t[5][1][i] = t[5][1][i] + p[1][pid];
  t[5][2][i] = t[5][2][i] + p[2][pid];
}

void condition_enviroment(unsigned int nb, iREAL *v[3], iREAL *angular[6], iREAL *rotation[9],
                        iREAL *mass, iREAL *inertia[9], iREAL *inverse[9], iREAL *invm, int *parmat)
{
  for(unsigned int i = 0; i < nb; i++)
  {
    iREAL rand = drand48();//random pull velocity
    
    //v[0][i] = 500 * rand;
    //v[1][i] = 500 * rand;
    //v[2][i] = 500 * rand;
    if(i == 0)
    {
      v[0][i] = 0;
      v[1][i] = 0;
      v[2][i] = 100;
    }else
    {
      v[0][i] = 0;
      v[1][i] = 0;
      v[2][i] = 0;
    }

    iREAL randomrot = 5 *rand;
    angular[0][i] = 0;
    angular[1][i] = 0;
    angular[2][i] = 0;
    angular[3][i] = 0; 
    angular[4][i] = 0;
    angular[5][i] = 0;
   
    mass[i] = 0.5;
    parmat[i] = 0;
    
    //double volume = (4./3.)*M_PI*rad*rad*rad;

    //mass[i] = volume*mparam[DENSITY][material];

    inertia[0][i] = 0.4*mass[i];
    inertia[1][i] = 0;
    inertia[2][i] = 0;
    inertia[3][i] = 0;
    inertia[4][i] = 0.4*mass[i];
    inertia[5][i] = 0;
    inertia[6][i] = 0;
    inertia[7][i] = 0;
    inertia[8][i] = 0.4*mass[i];

    rotation[0][i] = 1;
    rotation[1][i] = 0;
    rotation[2][i] = 0;
    rotation[3][i] = 0;
    rotation[4][i] = 1;
    rotation[5][i] = 0;
    rotation[6][i] = 0;
    rotation[7][i] = 0;
    rotation[8][i] = 1;
  }
  invert (nb, inertia, inverse, mass, invm);
}

void init_enviroment(unsigned int *nt, unsigned int *nb, 
                    iREAL *t[6][3], iREAL *v[3], iREAL *angular[6], iREAL *inertia[9], iREAL *inverse[9], iREAL *rotation[9], 
                    iREAL *mass, iREAL *invm, int *parmat, unsigned int tid[], unsigned int pid[], iREAL *position[6], iREAL lo[3], iREAL hi[3])
{
  //Input Type
  //0: Triangulated Mesh
  //1: Particle
  //2: Triangle
  //3: Sphere
  //4: Square
  //5: Hexahedron
  
    
  //load non-spherical meshes from vtk files
  /* *nParticles = 1000;
  int ptype[*nParticles];
  for(int i = 0; i < *nParticles; i++)
  {
    ptype[i] = 1;
  }*/
   
   
  //non-spherical particles generation and loading
  *nb = 2;
  int ptype[*nb];
  for(unsigned int i = 0; i < *nb; i++){ptype[i] = 6;}
  
  iREAL mint, maxt;
  load_enviroment(ptype, nt, *nb, t, tid, pid, position, &mint, &maxt);
  
  lo[0] = -250; // lower corner
  lo[1] = -250; // lower corner
  lo[2] = -250; // lower corner
  
  hi[0] = 250; // upper corner
  hi[1] = 250; // upper corner
  hi[2] = 250; // upper corner
  
  int radius = 10;

  unsigned int idx = 0; lo[0] =0; lo[1] = 0; lo[2] = 0;
  for(int ii = lo[0]; ii < hi[0]; ii=ii+radius)
  {
    for(int jj = lo[1]; jj < hi[1]; jj=jj+radius)
    {
      for(int kk = lo[2]; kk < hi[2]; kk=kk+radius)
      {
        position[0][idx] = ii+(radius/2);
        position[1][idx] = jj+(radius/2);
        position[2][idx] = kk+(radius/2);
        
        position[3][idx] = ii+(radius/2);
        position[4][idx] = jj+(radius/2);
        position[5][idx] = kk+(radius/2);
        
        //compute position to translate
        for(unsigned int j = 0; j < *nt; j++)
        {
          if(pid[j] == idx)
          {
            translate_enviroment(j, idx, t, position);
          }
        }
        idx++;
        if(idx > *nb) break;
      }
      if(idx > *nb) break;
    }
    if(idx > *nb) break;
  }

  condition_enviroment(*nb, v, angular, rotation, mass, inertia, inverse, invm, parmat);
}

void load_enviroment(int ptype[], unsigned int *nt, unsigned int nParticles, 
                    iREAL *t[6][3], unsigned int tid[], unsigned int pid[], iREAL *position[6], iREAL *mint, iREAL *maxt)
{
  unsigned int n = 0;
  *nt = 0;
  srand48(time(NULL));
  for(unsigned int i = 0; i < nParticles; i++)
  {
    if(ptype[i] != 6)
    {//load particle from file
        load_points(ptype[i], &n, i, *nt, t, tid, pid, position, mint, maxt);
    } else
    {//create point cloud and do delaunay hull triangulation
        //0.25 eps is the roundness degree, 5 is the radius, 50 are the point of the point cloud
        gen_nonsphericalparticle(0.25, 2.5, 50, &n, i, *nt, t, tid, pid, position, mint, maxt);
    }
    *nt = n + *nt;
    n = 0;
  }
}

void gen_nonsphericalparticle(iREAL eps, iREAL radius, int pointsize, unsigned int *nt, unsigned int bodyidx, unsigned int initidx, 
                              iREAL *t[6][3], unsigned int tid[], unsigned int pid[], iREAL *position[6], iREAL *mint, iREAL *maxt)
{
  iREAL v[100][3];
  for(int i = 0; i<pointsize; i++)
  {
    iREAL rng1 = 0;
    iREAL rng2 = 360;
    iREAL phi = (rng2-rng1) * (drand48()) + rng1;
    iREAL theta = (rng2-rng1) * (drand48()) + rng1;
    iREAL myradius = ((drand48()*eps)+1.0-eps) * radius;
            
    v[i][0] = myradius*sin(phi) * cos(theta);
    v[i][1] = myradius*sin(phi) * sin(theta);
    v[i][2] = myradius*cos(phi);
  }
    
  TRI* tr = NULL;
  free(tr);int pointlength = 0;
  tr = hull((iREAL *)v, pointsize, &pointlength);
  int counter = 0;
  
  for(TRI *tri = tr, *e = tri + pointlength; tri < e; tri ++){counter++;}
  
  unsigned int n = counter*3;
  *nt = pointlength;
  iREAL *point[3];
  
  point[0] = (iREAL *)malloc (n*sizeof(iREAL));
  point[1] = (iREAL *)malloc (n*sizeof(iREAL));
  point[2] = (iREAL *)malloc (n*sizeof(iREAL));
  iREAL min = DBL_MAX;
  iREAL max = DBL_MIN;
  
  counter = 0;
  for(TRI *tri = tr, *e = tri + pointlength; tri < e; tri ++)
  {
    point[0][counter] = tri->ver [0][0];
    point[1][counter] = tri->ver [0][1];
    point[2][counter] = tri->ver [0][2];
   
    if(point[0][counter] < min)
    {
        min = point[0][counter];
    }
    
    if(point[1][counter] < min)
    {
        min = point[1][counter];
    }
    
    if(point[2][counter] < min)
    {
        min = point[2][counter];
    }
    
    /////////////////////
    
    if(point[0][counter] > max)
    {
        max = point[0][counter];
    }
    
    if(point[1][counter] > max)
    {
        max = point[1][counter];
    }
    
    if(point[2][counter] > max)
    {
        max = point[2][counter];
    }
    
    counter++;
    point[0][counter] = tri->ver [1][0];
    point[1][counter] = tri->ver [1][1];
    point[2][counter] = tri->ver [1][2];
    
    
    if(point[0][counter] < min)
    {
        min = point[0][counter];
    }
    
    if(point[1][counter] < min)
    {
        min = point[1][counter];
    }
    
    if(point[2][counter] < min)
    {
        min = point[2][counter];
    }
    
    /////////////////////
    
    if(point[0][counter] > max)
    {
        max = point[0][counter];
    }
    
    if(point[1][counter] > max)
    {
        max = point[1][counter];
    }
    
    if(point[2][counter] > max)
    {
        max = point[2][counter];
    }
    
    counter++;
    point[0][counter] = tri->ver [2][0];
    point[1][counter] = tri->ver [2][1];
    point[2][counter] = tri->ver [2][2];
    
    if(point[0][counter] < min)
    {
        min = point[0][counter];
    }
    
    if(point[1][counter] < min)
    {
        min = point[1][counter];
    }
    
    if(point[2][counter] < min)
    {
        min = point[2][counter];
    }
    
    /////////////////////
    
    if(point[0][counter] > max)
    {
        max = point[0][counter];
    }
    
    if(point[1][counter] > max)
    {
        max = point[1][counter];
    }
    
    if(point[2][counter] > max)
    {
        max = point[2][counter];
    }
    
    counter++;
  }
  
  counter = 0;
  for(unsigned int i=initidx;i<initidx+n;i++)
  {
    //SPATIAL POINT A
    t[0][0][i] = point[0][counter];
    t[0][1][i] = point[1][counter];
    t[0][2][i] = point[2][counter];
  
    //REFERENTIAL POINT A
    t[3][0][i] = point[0][counter];
    t[3][1][i] = point[1][counter];
    t[3][2][i] = point[2][counter];

    counter++;
    //SPATIAL POINT B
    t[1][0][i] = point[0][counter];
    t[1][1][i] = point[1][counter];
    t[1][2][i] = point[2][counter];
    
    //REFENTIAL POINT B
    t[4][0][i] = point[0][counter];
    t[4][1][i] = point[1][counter];
    t[4][2][i] = point[2][counter];

    counter++;
    //SPATIAL POINT C
    t[2][0][i] = point[0][counter];
    t[2][1][i] = point[1][counter];
    t[2][2][i] = point[2][counter];
    
    //SPATIAL POINT C
    t[5][0][i] = point[0][counter];
    t[5][1][i] = point[1][counter];
    t[5][2][i] = point[2][counter];
    
    counter++;
    
    tid[i] = i;
    pid[i] = bodyidx;
  }
  
  getCentroid(bodyidx, initidx, initidx+n, t, position);
  
  *mint = min;
  *maxt = max;
}


void load_points(int ptype, unsigned int *nt, unsigned int bodyID, unsigned int startIDX, 
                iREAL *t[6][3], unsigned int tid[], unsigned int pid[], iREAL *position[6], iREAL *mint, iREAL *maxt)
{
  //////////VTK format////////////

  //Input Type
  //0: Triangulated Mesh
  //1: Particle
  //2: Triangle
  //3: Sphere
  //4: Square
  //5: Hexahedron
  //6: point cloud 
  iREAL min = DBL_MAX;
  iREAL max = DBL_MIN;
  FILE *fp1;
  if(ptype == 0)
  {
    char filename[100] = "input/mesh/mesh";
    char strtmp[100];
    sprintf(strtmp, "%i.vtk", bodyID);
    strcat(filename, strtmp);
    fp1 = fopen(filename, "r+");
    printf("%s\n", filename);
  } else if(ptype == 1)
  {
    char filename[100] = "input/particles/par_";
    char strtmp[100];
    sprintf(strtmp, "%i.vtk", bodyID);
    strcat(filename, strtmp);
    fp1 = fopen(filename, "r+");
    printf("%s\n", filename);
  } else if(ptype == 2)
  {
    fp1 = fopen("input/shapes/triangle.vtk", "r");
  } else if(ptype == 3)
  {
    fp1 = fopen("input/shapes/sphere.vtk", "r");
  } else if(ptype == 4)
  {
    fp1 = fopen("input/shapes/square.vtk", "r");
  } else if(ptype == 5)
  {
    fp1 = fopen("input/shapes/hexahedron.vtk", "r");
  }

  if( fp1 == NULL )
  {
      perror("Error while opening the file.\n");
      exit(EXIT_FAILURE);
  }
  
  char ch, word[500];
  iREAL *point[3];
  
  do {
      ch = fscanf(fp1,"%s",word);
      if(strcmp(word, "POINTS")==0)
      {
          printf("found!\n");
          ch = fscanf(fp1,"%s",word);
          unsigned int n = atol(word);
          //get points
          ch = fscanf(fp1,"%s",word);
          //printf("will read: %llu\n",n); 
          point[0] = (iREAL *)malloc (n*sizeof(iREAL));
          point[1] = (iREAL *)malloc (n*sizeof(iREAL));
          point[2] = (iREAL *)malloc (n*sizeof(iREAL));
          
          for(unsigned int i=0;i<n;i++)
          {
              fscanf(fp1, "%lf", &point[0][i]);
              fscanf(fp1, "%lf", &point[1][i]);
              fscanf(fp1, "%lf", &point[2][i]);
              //printf("POINT[0] = %f | POINT[1] = %f | POINT[2] = %f\n", point[0][i], point[1][i], point[2][i]);
              
              if(point[0][i] < min) 
              {
                min = point[0][i];
              }

              if(point[1][i] < min) 
              {
                min = point[1][i];
              }

              if(point[2][i] < min)
              {
                min = point[2][i];
              }

              /////////////////////
              
              if(point[0][i] > max) 
              {
                max = point[0][i];
              }

              if(point[1][i] > max) 
              {
                max = point[1][i];
              }

              if(point[2][i] > max)
              {
                max = point[2][i];
              }
          }
      }
      if(strcmp(word, "CELLS")==0 || strcmp(word, "POLYGONS") == 0)
      { 
          ch = fscanf(fp1,"%s",word);
          unsigned int n = atol(word);
          *nt = n;
          ch = fscanf(fp1,"%s",word);
          printf(":::%u::\n",n);
          for(unsigned int i=startIDX;i<startIDX+n;i++)
          {
              ch = fscanf(fp1,"%s",word);
              ch = fscanf(fp1,"%s",word);
              
              unsigned int index = atol(word);
              t[0][0][i] = point[0][index];
              t[0][1][i] = point[1][index];
              t[0][2][i] = point[2][index];
              
              //printf("idx:%s T[0][0] = %f | T[0][1] = %f | T[0][2] = %f\n", word, t[0][0][i], t[0][1][i], t[0][2][i]);
              
              ch = fscanf(fp1,"%s",word);
              index = atol(word);
              t[1][0][i] = point[0][index];
              t[1][1][i] = point[1][index];
              t[1][2][i] = point[2][index];
              
              //printf("idx:%s T[1][0] = %f | T[1][1] = %f | T[1][2] = %f\n", word, t[1][0][i], t[1][1][i], t[1][2][i]);
              
              ch = fscanf(fp1,"%s",word);
              index = atol(word);
              t[2][0][i] = point[0][index];
              t[2][1][i] = point[1][index];
              t[2][2][i] = point[2][index];
              
              //printf("idx:%s T[2][0] = %f | T[2][1] = %f | T[2][2] = %f\n", word, t[2][0][i], t[2][1][i], t[2][2][i]);
              
              tid[i] = i;
              pid[i] = bodyID;
          }
          getCentroid(bodyID, startIDX, startIDX+n, t, position);
      }
  } while (ch != EOF);
  *mint = min;
  *maxt = max;
  fclose(fp1);
}

void normalize(unsigned int nt, iREAL *t[6][3], iREAL mint, iREAL maxt) 
{  
  //range -255 to 255
  iREAL inv_range = 510.0/(maxt-mint);

  for(unsigned int i=0; i < nt; i++)
  {
    for(int j=0; j < 3; j++)
    {
      for(int k=0; k < 3; k++)
      {
        t[j][k][i] = -255.0 + (t[j][k][i] - mint) * inv_range;
      }
    }
  }
}

/*
// temporary surface pairing 
struct pair
{
  int color1;
  int color2;
  REAL iparam[NINT];
};

// material comparison by color pair 
struct cmp
{
  bool operator() (const pair& a, const pair& b)
  {
    if (a.color1 == b.color1) return a.color2 < b.color2;
    else return a.color1 < b.color1;
  }
};

// sort materials according to surface color pairs 
static void sort_materials ()
{
  std::vector<pair> v;

  v.reserve (pairnum);

  for (int i = 0; i < pairnum; i++)
  {
    pair x;

    x.color1 = pairs[2*i];
    x.color2 = pairs[2*i+1];
    for (int j = 0; j < NINT; j ++) x.iparam[j] = iparam[j][i];

    v.push_back (x);
  }

  std::sort (v.begin(), v.end(), cmp());

  int i = 0;

  for (std::vector<pair>::const_iterator p = v.begin(); p != v.end(); ++p, ++i)
  {
    pairs[2*i] = p->color1;
    pairs[2*i+1] = p->color2;
    for (int j = 0; j < NINT; j ++) iparam[j][i] = p->iparam[j];
  }
}

// create material 
static PyObject* MATERIAL (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("density", "young", "poisson");
  double density, young, poisson;

  PARSEKEYS ("ddd", &density, &young, &poisson);

  TYPETEST (is_positive (density, kwl[0]) && is_positive (young, kwl[1]) && is_positive (poisson, kwl[2]));

  if (matnum >= material_buffer_size) material_buffer_grow ();

  int i = matnum ++;

  mparam[DENSITY][i] = density;
  mparam[YOUNG][i] = young;
  mparam[POISSON][i] = poisson;

  return PyInt_FromLong (i);
}

// define surface pairing 
static PyObject* PAIRING (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("color1", "color2", "model", "spring", "damper", "friction", "rolling", "drilling",
            "kskn", "uforce", "lambda", "young2", "kskn2", "sigc", "tauc", "alpha");
  int color1, color2;
  double spring, damper, rolling, drilling, kskn, lambda, young2, kskn2, sigc, tauc, alpha;
  PyObject *model, *friction, *ufo;
  double fri[2];
  int iki;

  spring = 1E6;
  damper = 1.0;
  rolling = 0.0;
  drilling = 0.0;
  kskn = 0.5;
  lambda = 1.0;
  young2 = 0.0;
  kskn2 = 0.5;
  sigc = 0.0;
  tauc = 0.0;
  alpha = 0.0;
  friction = NULL;
  fri[0] = 0.0;
  fri[1] = 0.0;
  ufo = NULL;

  PARSEKEYS ("iiO|ddOdddOdddddd", &color1, &color2, &model, &spring, &damper, &friction,
    &rolling, &drilling, &kskn, &ufo, &lambda, &young2, &kskn2, &sigc, &tauc, &alpha);

  TYPETEST (is_string (model, kwl[2]) && is_positive (spring, kwl[3]) && is_non_negative (damper, kwl[4]) &&
    is_non_negative (rolling, kwl[6]) && is_non_negative (drilling, kwl[7]) && is_positive (kskn, kwl[8]) &&
    is_callable (ufo, kwl[9]) && is_positive (lambda, kwl[10]) && is_non_negative (young2, kwl[11]) &&
    is_positive (kskn2, kwl[12]) && is_non_negative (sigc, kwl[13]) && is_non_negative (tauc, kwl[14]) &&
    is_non_negative (alpha, kwl[15]));

  IFIS (model, "granular")
  {
    iki = GRANULAR;
  }
  ELIF (model, "bonded")
  {
    iki = BONDED;
  }
  ELIF (model, "user")
  {
    iki = UFORCE;
  }
  ELSE
  {
    PyErr_SetString (PyExc_ValueError, "Invalid interaction model");
    return NULL;
  }

  if (PyTuple_Check (friction))
  {
    TYPETEST (is_tuple (friction, kwl[5], 2));

    fri[0] = PyFloat_AsDouble (PyTuple_GetItem (friction,0));
    fri[1] = PyFloat_AsDouble (PyTuple_GetItem (friction,1));

    TYPETEST (is_non_negative (fri[0], kwl[5]) && is_non_negative (fri[1], kwl[5]));
  }
  else
  {
    fri[0] = PyFloat_AsDouble (friction);

    TYPETEST (is_non_negative (fri[0], kwl[5]));

    fri[1] = fri[0];
  }

  int i = -1;

  if (color1 == 0 && color2 == 0) // default
  {
    i = 0;
  }
  else if (color1 > 0 && color2 > 0)
  {
    if (pairnum >= pair_buffer_size) pair_buffer_grow ();

    int i = pairnum ++;

    pairs[2*i] = std::min(color1, color2);
    pairs[2*i+1] = std::max(color1, color2);
  }
  else
  {
    PyErr_SetString (PyExc_ValueError, "invalid color values");
    return NULL;
  }

  ikind[i] = iki;
  iparam[SPRING][i] = spring;
  iparam[DAMPER][i] = damper;
  iparam[FRISTAT][i] = fri[0];
  iparam[FRIDYN][i] = fri[1];
  iparam[FRIROL][i] = rolling;
  iparam[FRIDRIL][i] = drilling;
  iparam[KSKN][i] = kskn;
  iparam[LAMBDA][i] = lambda;
  iparam[YOUNG2][i] = young2;
  iparam[KSKN2][i] = kskn2;
  iparam[SIGC][i] = sigc;
  iparam[TAUC][i] = tauc;
  iparam[ALPHA][i] = alpha;
  uforce[i] = (callback_t) ufo;

  Py_RETURN_NONE;
}

// set particle velocity 
static PyObject* VELOCITY (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("parnum", "linear", "angular");
  PyObject *lin, *ang;
  int i;

  ang = NULL;

  PARSEKEYS ("iO|O", &i, &lin, &ang);

  TYPETEST (is_in_range (i, 0, parnum, kwl[0]) && is_tuple (lin, kwl[1], 3) && is_tuple (ang, kwl[2], 3));

  linear[0][i] = PyFloat_AsDouble (PyTuple_GetItem (lin, 0)); 
  linear[1][i] = PyFloat_AsDouble (PyTuple_GetItem (lin, 1)); 
  linear[2][i] = PyFloat_AsDouble (PyTuple_GetItem (lin, 2)); 

  if (ang)
  {
    angular[0][i] = PyFloat_AsDouble (PyTuple_GetItem (ang, 0)); 
    angular[1][i] = PyFloat_AsDouble (PyTuple_GetItem (ang, 1)); 
    angular[2][i] = PyFloat_AsDouble (PyTuple_GetItem (ang, 2)); 
  }

  Py_RETURN_NONE;
}

// set gravity
static PyObject* GRAVITY (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("gx", "gy", "gz");
  double gx, gy, gz;

  PARSEKEYS ("ddd", &gx, &gy, &gz);

  gravity[0] = gx;
  gravity[1] = gy;
  gravity[2] = gz;

  Py_RETURN_NONE;
}*/
