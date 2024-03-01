#include "mex.h" //add for mex MMW May 2019
#include <cstdint> //add for mex MMW May 2019

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define M_PI 3.14159265358979323846

#define KARRAY_C(i,j) karray_c[(i) + (j)*karray_dims[0]] //4D

double haltonnumber(int index,int base)
{
  double result = 0;
  double f = 1.0;
  int i = index;
  while (i>0)
  {
    f = f/base;
    result += f*fmod(i,base);
    i = i/base;
  }
  return result;
}

void haltonSeq(double *m_adAzimuthalAngle,double *m_adPolarAngle,long num_frames,long num_projPerFrame)
{
  int p1 = 2;
  int p2 = 3;
  double z;
  double phi;
  int linter;
  int lk;                  // index for projections in a frame
  int lFrame;              // index for number of frames

  for (lFrame = 0; lFrame < num_frames; lFrame ++ )
  {
    for (lk = 0; lk < num_projPerFrame; lk ++ )
    {
      linter  = lk + lFrame * num_projPerFrame;
      z = haltonnumber(lk+1,p1) * 2 - 1;
      phi = 2 * M_PI * haltonnumber(lk+1,p2);
      //calculate polar and azimuthal angle on a unit sphere
      m_adPolarAngle[linter] = acos(z);
      m_adAzimuthalAngle[linter] = phi;
    }
  }
}

void spiralSeq(double *m_azi,double *m_polar,long num_frames,long num_projPerFrame)
{
  long   lk;
  long   llin;                // linearly increasing index
  long   linter;              // interleaved increasing index
  long   lFrame;              // index for number of frames
  double dPreviousAngle;      // previous angle value
  double dH;                  // parameter for the calculation

  long num_totalProjections = num_frames * num_projPerFrame;

  dPreviousAngle = 0;

  for (lk = 0; lk < num_projPerFrame; lk ++ )
  {
    for (lFrame = 0; lFrame < num_frames; lFrame ++ )
    {
      llin    = lFrame + lk * num_frames;
      linter  = lk + lFrame * num_projPerFrame;

      dH = -1.0 + 2.0 * llin / (double)num_totalProjections;
      m_polar[linter] = acos (dH);
      if (llin == 0)
      m_azi[linter] = 0;
      else
      m_azi[linter] = fmod ( dPreviousAngle + 3.6/ ( sqrt( num_totalProjections * (1. - dH*dH) ) ) , 2.0 * M_PI );

      dPreviousAngle = m_azi[linter];

    }//endfor lk
  }//endfor lInt
}

void archimedianSeq(double *m_azi,double *m_polar,long num_frames,long num_projPerFrame)
{
  long   lk;
  long   linter;              // interleaved increasing index
  long   lFrame;              // index for number of frames
  double dZ;
  double dAngle;

  dAngle = (3.0-sqrt(5.0))*M_PI;
  dZ = 2.0/(num_projPerFrame - 1.0);

  for (lFrame = 0; lFrame < num_frames; lFrame ++ )
  {
    for (lk = 0; lk < num_projPerFrame; lk ++ )
    {
      linter  = lk + lFrame * num_projPerFrame;
      m_polar[linter] = acos(1.0 - dZ*lk);
      m_azi[linter] = lk*dAngle;
    }
  }
}

void dgoldenMSeq(double *m_azi,double *m_polar,long num_frames,long num_projPerFrame)
{
  long   lk;
  long   linter;              // interleaved increasing index
  long   lFrame;              // index for number of frames
  double goldmean1;
  double goldmean2;
  goldmean1 = 0.465571231876768;
  goldmean2 = 0.682327803828019;

  for (lFrame = 0; lFrame < num_frames; lFrame ++ )
  {
    for (lk = 0; lk < num_projPerFrame; lk ++ )
    {
      linter  = lk + lFrame * num_projPerFrame;
      m_polar[linter] = acos(2.0 * fmod(lk * goldmean1, 1) - 1);
      m_azi[linter] = 2 * M_PI * fmod(lk * goldmean2, 1);
    }
  }
}

// Functions for Random Spiral
void swap (double *a, double *b)
{
    double t = *a;
    *a = *b;
    *b = t;
}

long partition(double *ht_polar, double *sp_polar, double *sp_azi, long low, long high)
{
    double pivot = ht_polar[high];
    long i = low - 1;

    for(long j = low; j<=high-1;j++)
    {
        if(ht_polar[j] <= pivot)
        {
            i++;
            swap(&ht_polar[i],&ht_polar[j]);
            swap(&sp_polar[i],&sp_polar[j]);
            swap(&sp_azi[i],&sp_azi[j]);
        }
    }
    swap(&ht_polar[i+1],&ht_polar[high]);
    swap(&sp_polar[i+1],&sp_polar[high]);
    swap(&sp_azi[i+1],&sp_azi[high]);

    return (i+1);
}

void quickSort(double *ht_polar, double *sp_polar, double *sp_azi, long low, long high)
{
    if(low < high)
    {
        long pi = partition(ht_polar,sp_polar,sp_azi,low,high);

        quickSort(ht_polar,sp_polar,sp_azi,low,pi-1);
        quickSort(ht_polar,sp_polar,sp_azi,pi+1,high);
    }
}

// Ziyi function randomizes spiral traj with Halton sequence
void randomSpiral(double *m_adAzimuthalAngle,double *m_adPolarAngle,long num_projPerFrame)
{
  double* ht_adAzimu = (double*) calloc(num_projPerFrame, sizeof(double));
  double* ht_adPolar = (double*) calloc(num_projPerFrame, sizeof(double));

  haltonSeq(ht_adAzimu,ht_adPolar,1,num_projPerFrame);

  quickSort(ht_adPolar,m_adPolarAngle,m_adAzimuthalAngle,0,num_projPerFrame-1);

  free(ht_adPolar);
  free(ht_adAzimu);
}

// Function generates trajectory according to the type
void gen_traj(long m_lProjectionsPerFrame, long m_lTrajectoryType, mxArray** karray)
// input:
// m_lProjectionsPerFrame: total number of radial views
// m_lTrajectoryType: flag for differnet trajectory types: 1: Spiral, 2: Halton, 3: Randomized Spiral using Halton, 4:Archimedian Spiral, 5: double golden mean.
// We use #3 for our clinical acquistion   
{
  const mwSize *_karray_dims;
  _karray_dims = mxGetDimensions(*karray);
  double *karray_c = mxGetPr(*karray);

  int karray_dims[3];
  karray_dims[0] = _karray_dims[0];
  karray_dims[1] = _karray_dims[1];
  karray_dims[2] = _karray_dims[2];

  //double m_adAzimuthalAngle = NULL;
  //double m_adPolarAngle = NULL;
  //double coordinates = NULL;
 
  double* m_adAzimuthalAngle =(double*) calloc(m_lProjectionsPerFrame, sizeof(double));
  double* m_adPolarAngle =(double*) calloc(m_lProjectionsPerFrame, sizeof(double));
  double* coordinates =(double*) calloc(m_lProjectionsPerFrame*3, sizeof(double));
  
  //m_adAzimuthalAngle = (double*) calloc(m_lProjectionsPerFrame, sizeof(double)); 
  //m_adPolarAngle = (double*) calloc(m_lProjectionsPerFrame, sizeof(double));
  //coordinates = (double*) calloc(m_lProjectionsPerFrame*3, sizeof(double));

  long m_lNumberOfFrames = 1;

  switch(m_lTrajectoryType)
  {
    //Ziyi: siemens spiral
    case 1: {spiralSeq(m_adAzimuthalAngle,m_adPolarAngle,m_lNumberOfFrames,m_lProjectionsPerFrame);break;}

    //Ziyi, case 2: Halton sequence
    case 2: {haltonSeq(m_adAzimuthalAngle,m_adPolarAngle,m_lNumberOfFrames,m_lProjectionsPerFrame);break;}

    //Ziyi: case 3 haltoned spiral
    case 3:
    {
      spiralSeq(m_adAzimuthalAngle,m_adPolarAngle,m_lNumberOfFrames,m_lProjectionsPerFrame);
      randomSpiral(m_adAzimuthalAngle,m_adPolarAngle,m_lProjectionsPerFrame);
      break;
    }

    //Ziyi, case 4: archimedian Seq
    case 4: {archimedianSeq(m_adAzimuthalAngle,m_adPolarAngle,m_lNumberOfFrames,m_lProjectionsPerFrame);break;}

    //Ziyi, case 5: double golden means
    default: {dgoldenMSeq(m_adAzimuthalAngle,m_adPolarAngle,m_lNumberOfFrames,m_lProjectionsPerFrame);break;}
  }

  //convert angles to x, y ,z coordinate
  for(int k = 0; k<m_lProjectionsPerFrame; k++)
  {
    //Xï¼ŒYï¼ŒZ coordinates
      //PJN - just save the azimuthal and polar angles for ease of calculating trajectories after the fact.
    KARRAY_C(k,0) = m_adAzimuthalAngle[k];
    KARRAY_C(k,1) = m_adPolarAngle[k];
    coordinates[k] = sin(m_adPolarAngle[k])*cos(m_adAzimuthalAngle[k]);
    coordinates[k+m_lProjectionsPerFrame] = sin(m_adPolarAngle[k])*sin(m_adAzimuthalAngle[k]);
    coordinates[k+m_lProjectionsPerFrame*2] = cos(m_adPolarAngle[k]);
  }

  free(m_adAzimuthalAngle);
  free(m_adPolarAngle);

//  return(coordinates);
}


void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    long NPro  = mxGetScalar(prhs[0]);
    long Traj_Type = mxGetScalar(prhs[1]);
    //1: Spiral, 2: Halton, 3: Randomized Spiral using Halton, 4:Archimedian Spiral, 5: double golden mean
    
    //Create karrays via mex
    mwSize *mwk_dims;//create pointer
    mwk_dims = (mwSize *) mxMalloc (2 * sizeof (mwSize)); //allocate memory
    mwk_dims[0] = NPro; //Number of Projections
    mwk_dims[1] = 2; //Azimuthal and Polar Angles


    mxArray *karray = mxCreateNumericArray (2, mwk_dims, mxDOUBLE_CLASS, mxREAL);

  
    gen_traj(NPro, Traj_Type, &karray);
    
    plhs[0] = karray;

}
