/***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mex.h>
#include <math.h>
#include <time.h>

#include "cuda.h"
#include "device_functions.h"
#include "CoReg_TVL1_Newton_Ind_kernels.cu"

#define BLOCKSIZE 512
#define MAX_GRID_SIZE 65535
#define NUMTHREADS 512

dim3 GetGrid(int size){
    size = (size-1) / NUMTHREADS + 1;
    dim3 grid( size, 1, 1 );
    if( grid.x > MAX_GRID_SIZE ) grid.x = grid.y = (int) sqrt( (double)(size-1) ) + 1;
    else if( grid.y > MAX_GRID_SIZE ) grid.x = grid.y = grid.z = (int) pow( (double)(size-1), (double)1.0/3.0 ) + 1;
    return grid;
}

extern void mexFunction(int iNbOut, mxArray *pmxOut[],
        int iNbIn, const mxArray *pmxIn[]){
    
    /* iNbOut: number of outputs */
    /* pmxOut: array of pointers to output arguments */
    
    /* iNbIn: number of inputs
    /* pmxIn: array of pointers to input arguments */
    
    /*  host arrays and variables */
    float   *h_u1x, *h_u1y, *h_u1z, *h_cvg, *h_U1x, *h_U1y, *h_U1z;
    float   *h_u2x, *h_u2y, *h_u2z, *h_U2x, *h_U2y, *h_U2z;
    float   *h_VecParameters,*h_G1x, *h_G1y, *h_G1z, *h_G1f, *h_G1t;
    float   *h_G2x, *h_G2y, *h_G2z, *h_G2f, *h_G2t;
    float   *h_b1x1, *h_b1x2, *h_b1x3, *h_b1y1, *h_b1y2, *h_b1y3, *h_b1z1, *h_b1z2, *h_b1z3;
    float   *h_b2x1, *h_b2x2, *h_b2x3, *h_b2y1, *h_b2y2, *h_b2y3, *h_b2z1, *h_b2z2, *h_b2z3;
    float   *h_q1, *h_q2, *h_gkx, *h_gky, *h_gkz, *tt, *h_dv1x, *h_dv1y, *h_dv1z;
    float   *h_dv2x, *h_dv2y, *h_dv2z;
    float   fError, cc, steps, fPenalty1, fPenalty2, fPenalty3, fps;
    float   weight1, weight2;
    /*
    int     *punum, iNy, iNx, iNz, iNdim, iDim[3], iNI;
    int     maxIter, SZF, iDev;
    */
    int     *punum, iNy, iNx, iNz, iNdim, iDim[3], maxIter;
    int     flag = 1;
    
    
    cudaSetDevice(1);
    
    /* Timing */
    cudaEvent_t start, stop;
    float time;
    
    /*  device arrays */
    float   *d_b1x1, *d_b1y1, *d_b1z1, *d_b1x2, *d_b1y2, *d_b1z2, *d_b1x3, *d_b1y3, *d_b1z3, *d_dv1x, *d_dv1y, *d_dv1z;
    float   *d_b2x1, *d_b2y1, *d_b2z1, *d_b2x2, *d_b2y2, *d_b2z2, *d_b2x3, *d_b2y3, *d_b2z3, *d_dv2x, *d_dv2y, *d_dv2z;

    float   *d_q1, *d_q2;
    float   *d_u1x, *d_u1y, *d_u1z, *d_U1x, *d_U1y, *d_U1z, *d_G1x, *d_G1y, *d_G1z, *d_G1f, *d_G1t;
    float   *d_u2x, *d_u2y, *d_u2z, *d_U2x, *d_U2y, *d_U2z, *d_G2x, *d_G2y, *d_G2z, *d_G2f, *d_G2t;
    float   *d_gkx, *d_gky, *d_gkz;

    float   *h_FPS, *d_FPS;
    float   *h_qq1, *h_qq2,*h_qq3, *d_qq1, *d_qq2, *d_qq3;
    
    /* CUDA event-based timer start */
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord( start, 0 );
    
    /* input interface with matlab arrays */
    h_VecParameters = (float *)mxGetData(pmxIn[0]); /* Vector of parameters */
    h_U1x = (float *)mxGetData(pmxIn[1]);
    h_U1y = (float *)mxGetData(pmxIn[2]);
    h_U1z = (float *)mxGetData(pmxIn[3]);
    h_G1x = (float *)mxGetData(pmxIn[4]);
    h_G1y = (float *)mxGetData(pmxIn[5]);
    h_G1z = (float *)mxGetData(pmxIn[6]);
    h_G1t = (float *)mxGetData(pmxIn[7]);
    h_G1f = (float *)mxGetData(pmxIn[8]);

    h_U2x = (float *)mxGetData(pmxIn[9]);
    h_U2y = (float *)mxGetData(pmxIn[10]);
    h_U2z = (float *)mxGetData(pmxIn[11]);
    h_G2x = (float *)mxGetData(pmxIn[12]);
    h_G2y = (float *)mxGetData(pmxIn[13]);
    h_G2z = (float *)mxGetData(pmxIn[14]);
    h_G2t = (float *)mxGetData(pmxIn[15]);
    h_G2f = (float *)mxGetData(pmxIn[16]);
    
    
    /* dimensions */
    iNy = (int) h_VecParameters[0];
    iNx = (int) h_VecParameters[1];
    iNz = (int) h_VecParameters[2];
    
    unsigned int imageSize = iNx*iNy*iNz;
    
    /* parameters */
    maxIter = (int) h_VecParameters[3]; /* total number of iterations */
    fError = (float) h_VecParameters[4]; /* error criterion */
    cc = (float) h_VecParameters[5]; /* cc for ALM */
    steps = (float) h_VecParameters[6]; /* steps for each iteration */
    fPenalty1 = (float) h_VecParameters[7];
    fPenalty2 = (float) h_VecParameters[8];
    fPenalty3 = (float) h_VecParameters[9];
    weight1 = (float) h_VecParameters[10];
    weight2 = (float) h_VecParameters[11];
    
    /* output interface with matlab */
    /* u1x */
    iNdim = 3;
    iDim[0] = iNy;
    iDim[1] = iNx;
    iDim[2] = iNz;
    
    pmxOut[0] = mxCreateNumericArray(iNdim,(const int*)iDim,mxSINGLE_CLASS,mxREAL);
    h_u1x = (float*)mxGetData(pmxOut[0]);
    
    /* u1y */
    iNdim = 3;
    iDim[0] = iNy;
    iDim[1] = iNx;
    iDim[2] = iNz;
    
    
    pmxOut[1] = mxCreateNumericArray(iNdim,(const int*)iDim,mxSINGLE_CLASS,mxREAL);
    h_u1y = (float*)mxGetData(pmxOut[1]);
    
    /* u1z */
    iNdim = 3;
    iDim[0] = iNy;
    iDim[1] = iNx;
    iDim[2] = iNz;
    
    pmxOut[2] = mxCreateNumericArray(iNdim,(const int*)iDim,mxSINGLE_CLASS,mxREAL);
    h_u1z = (float*)mxGetData(pmxOut[2]);

    /* u2x */
    iNdim = 3;
    iDim[0] = iNy;
    iDim[1] = iNx;
    iDim[2] = iNz;
    
    pmxOut[3] = mxCreateNumericArray(iNdim,(const int*)iDim,mxSINGLE_CLASS,mxREAL);
    h_u2x = (float*)mxGetData(pmxOut[3]);
    
    /* u2y */
    iNdim = 3;
    iDim[0] = iNy;
    iDim[1] = iNx;
    iDim[2] = iNz;
    
    
    pmxOut[4] = mxCreateNumericArray(iNdim,(const int*)iDim,mxSINGLE_CLASS,mxREAL);
    h_u2y = (float*)mxGetData(pmxOut[4]);
    
    /* u2z */
    iNdim = 3;
    iDim[0] = iNy;
    iDim[1] = iNx;
    iDim[2] = iNz;
    
    pmxOut[5] = mxCreateNumericArray(iNdim,(const int*)iDim,mxSINGLE_CLASS,mxREAL);
    h_u2z = (float*)mxGetData(pmxOut[5]);
    
    /* convergence rate */
    iNdim = 2;
    iDim[0] = 1;
    iDim[1] = maxIter;
    pmxOut[6] = mxCreateNumericArray(iNdim,(const int*)iDim,mxSINGLE_CLASS,mxREAL);
    h_cvg = (float*)mxGetData(pmxOut[6]);
    
    /* number of iterations */
    iNdim = 2;
    iDim[0] = 1;
    iDim[1] = 1;
    pmxOut[7] = mxCreateNumericArray(iNdim,(const int*)iDim,mxUINT16_CLASS,mxREAL);
    punum = (int*)mxGetData(pmxOut[7]);
    
    /* computation time */
    iNdim = 2;
    iDim[0] = 1;
    iDim[1] = 1;
    pmxOut[8] = mxCreateNumericArray(iNdim,(const int*)iDim,mxSINGLE_CLASS,mxREAL);
    tt = (float*)mxGetData(pmxOut[8]);
    
    /* allocate host memory */
    /* bx1, bx2, bx3 */
    h_b1x1 = (float *) calloc( (unsigned)imageSize, sizeof(float) );
    h_b1x2 = (float *) calloc( (unsigned)imageSize, sizeof(float) );
    h_b1x3 = (float *) calloc( (unsigned)imageSize, sizeof(float) );
    if (!h_b1x1 || !h_b1x2 || !h_b1x3) mexPrintf("calloc: Memory allocation failure\n");
    
    /* by1, by2, by3 */
    h_b1y1 = (float *) calloc( (unsigned)imageSize, sizeof(float) );
    h_b1y2 = (float *) calloc( (unsigned)imageSize, sizeof(float) );
    h_b1y3 = (float *) calloc( (unsigned)imageSize, sizeof(float) );
    if (!h_b1y1 || !h_b1y2 || !h_b1y3) mexPrintf("calloc: Memory allocation failure\n");
    
    /* bz1, bz2, bz3 */
    h_b1z1 = (float *) calloc( (unsigned)imageSize, sizeof(float) );
    h_b1z2 = (float *) calloc( (unsigned)imageSize, sizeof(float) );
    h_b1z3 = (float *) calloc( (unsigned)imageSize, sizeof(float) );
    if (!h_b1z1 || !h_b1z2 || !h_b1z3) mexPrintf("calloc: Memory allocation failure\n");

    /* bx1, bx2, bx3 */
    h_b2x1 = (float *) calloc( (unsigned)imageSize, sizeof(float) );
    h_b2x2 = (float *) calloc( (unsigned)imageSize, sizeof(float) );
    h_b2x3 = (float *) calloc( (unsigned)imageSize, sizeof(float) );
    if (!h_b2x1 || !h_b2x2 || !h_b2x3) mexPrintf("calloc: Memory allocation failure\n");
    
    /* by1, by2, by3 */
    h_b2y1 = (float *) calloc( (unsigned)imageSize, sizeof(float) );
    h_b2y2 = (float *) calloc( (unsigned)imageSize, sizeof(float) );
    h_b2y3 = (float *) calloc( (unsigned)imageSize, sizeof(float) );
    if (!h_b2y1 || !h_b2y2 || !h_b2y3) mexPrintf("calloc: Memory allocation failure\n");
    
    /* bz1, bz2, bz3 */
    h_b2z1 = (float *) calloc( (unsigned)imageSize, sizeof(float) );
    h_b2z2 = (float *) calloc( (unsigned)imageSize, sizeof(float) );
    h_b2z3 = (float *) calloc( (unsigned)imageSize, sizeof(float) );
    if (!h_b2z1 || !h_b2z2 || !h_b2z3) mexPrintf("calloc: Memory allocation failure\n");
    
    
    /* q1 */
    h_q1 = (float *) calloc( (unsigned)imageSize, sizeof(float) );
    if (!h_q1) mexPrintf("calloc: Memory allocation failure\n");
    
    /* q2 */
    h_q2 = (float *) calloc( (unsigned)imageSize, sizeof(float) );
    if (!h_q2) mexPrintf("calloc: Memory allocation failure\n");
    
    /* gk */
    h_gkx = (float *) calloc( (unsigned)imageSize, sizeof(float) );
    h_gky = (float *) calloc( (unsigned)imageSize, sizeof(float) );
    h_gkz = (float *) calloc( (unsigned)imageSize, sizeof(float) );
    if (!h_gkx || !h_gky || !h_gkz) mexPrintf("calloc: Memory allocation failure\n");
    
    /* div */
    h_dv1x = (float *) calloc( (unsigned)imageSize, sizeof(float) );
    h_dv1y = (float *) calloc( (unsigned)imageSize, sizeof(float) );
    h_dv1z = (float *) calloc( (unsigned)imageSize, sizeof(float) );
    if (!h_dv1x || !h_dv1y || !h_dv1z ) mexPrintf("calloc: Memory allocation failure\n");

    /* div */
    h_dv2x = (float *) calloc( (unsigned)imageSize, sizeof(float) );
    h_dv2y = (float *) calloc( (unsigned)imageSize, sizeof(float) );
    h_dv2z = (float *) calloc( (unsigned)imageSize, sizeof(float) );
    if (!h_dv2x || !h_dv2y || !h_dv2z ) mexPrintf("calloc: Memory allocation failure\n");
    

    /*I am here now~~. continue tomorrow.*/
    /* h_FPS */
    h_FPS = (float *) calloc( (unsigned)imageSize, sizeof(float) );
    if (!h_FPS) mexPrintf("calloc: Memory allocation failure\n");

    /* h_qq1 */
    h_qq1 = (float *) calloc( (unsigned)imageSize, sizeof(float) );
    if (!h_qq1) mexPrintf("calloc: Memory allocation failure\n");
    
    /* h_qq2 */
    h_qq2 = (float *) calloc( (unsigned)imageSize, sizeof(float) );
    if (!h_qq2) mexPrintf("calloc: Memory allocation failure\n");

    /* h_qq3 */
    h_qq3 = (float *) calloc( (unsigned)imageSize, sizeof(float) );
    if (!h_qq3) mexPrintf("calloc: Memory allocation failure\n");
    
    
    
    /* device memory allocation */
    cudaMalloc( (void**) &d_b1x1, sizeof(float)*(unsigned)imageSize);
    cudaMalloc( (void**) &d_b1x2, sizeof(float)*(unsigned)imageSize);
    cudaMalloc( (void**) &d_b1x3, sizeof(float)*(unsigned)imageSize);
    cudaMalloc( (void**) &d_b1y1, sizeof(float)*(unsigned)imageSize);
    cudaMalloc( (void**) &d_b1y2, sizeof(float)*(unsigned)imageSize);
    cudaMalloc( (void**) &d_b1y3, sizeof(float)*(unsigned)imageSize);
    cudaMalloc( (void**) &d_b1z1, sizeof(float)*(unsigned)imageSize);
    cudaMalloc( (void**) &d_b1z2, sizeof(float)*(unsigned)imageSize);
    cudaMalloc( (void**) &d_b1z3, sizeof(float)*(unsigned)imageSize);
    cudaMalloc( (void**) &d_gkx, sizeof(float)*(unsigned)imageSize);
    cudaMalloc( (void**) &d_gky, sizeof(float)*(unsigned)imageSize);
    cudaMalloc( (void**) &d_gkz, sizeof(float)*(unsigned)imageSize);
    cudaMalloc( (void**) &d_dv1x, sizeof(float)*(unsigned)imageSize);
    cudaMalloc( (void**) &d_dv1y, sizeof(float)*(unsigned)imageSize);
    cudaMalloc( (void**) &d_dv1z, sizeof(float)*(unsigned)imageSize);
    cudaMalloc( (void**) &d_q1,  sizeof(float)*(unsigned)imageSize);
    cudaMalloc( (void**) &d_u1x, sizeof(float)*(unsigned)imageSize);
    cudaMalloc( (void**) &d_u1y, sizeof(float)*(unsigned)imageSize);
    cudaMalloc( (void**) &d_u1z, sizeof(float)*(unsigned)imageSize);
    cudaMalloc( (void**) &d_U1x, sizeof(float)*(unsigned)imageSize);
    cudaMalloc( (void**) &d_U1y, sizeof(float)*(unsigned)imageSize);
    cudaMalloc( (void**) &d_U1z, sizeof(float)*(unsigned)imageSize);
    cudaMalloc( (void**) &d_G1x, sizeof(float)*(unsigned)imageSize);
    cudaMalloc( (void**) &d_G1y, sizeof(float)*(unsigned)imageSize);
    cudaMalloc( (void**) &d_G1z, sizeof(float)*(unsigned)imageSize);
    cudaMalloc( (void**) &d_G1t, sizeof(float)*(unsigned)imageSize);
    cudaMalloc( (void**) &d_G1f, sizeof(float)*(unsigned)imageSize);

    cudaMalloc( (void**) &d_b2x1, sizeof(float)*(unsigned)imageSize);
    cudaMalloc( (void**) &d_b2x2, sizeof(float)*(unsigned)imageSize);
    cudaMalloc( (void**) &d_b2x3, sizeof(float)*(unsigned)imageSize);
    cudaMalloc( (void**) &d_b2y1, sizeof(float)*(unsigned)imageSize);
    cudaMalloc( (void**) &d_b2y2, sizeof(float)*(unsigned)imageSize);
    cudaMalloc( (void**) &d_b2y3, sizeof(float)*(unsigned)imageSize);
    cudaMalloc( (void**) &d_b2z1, sizeof(float)*(unsigned)imageSize);
    cudaMalloc( (void**) &d_b2z2, sizeof(float)*(unsigned)imageSize);
    cudaMalloc( (void**) &d_b2z3, sizeof(float)*(unsigned)imageSize);
    cudaMalloc( (void**) &d_dv2x, sizeof(float)*(unsigned)imageSize);
    cudaMalloc( (void**) &d_dv2y, sizeof(float)*(unsigned)imageSize);
    cudaMalloc( (void**) &d_dv2z, sizeof(float)*(unsigned)imageSize);
    cudaMalloc( (void**) &d_q2,  sizeof(float)*(unsigned)imageSize);
    cudaMalloc( (void**) &d_u2x, sizeof(float)*(unsigned)imageSize);
    cudaMalloc( (void**) &d_u2y, sizeof(float)*(unsigned)imageSize);
    cudaMalloc( (void**) &d_u2z, sizeof(float)*(unsigned)imageSize);
    cudaMalloc( (void**) &d_U2x, sizeof(float)*(unsigned)imageSize);
    cudaMalloc( (void**) &d_U2y, sizeof(float)*(unsigned)imageSize);
    cudaMalloc( (void**) &d_U2z, sizeof(float)*(unsigned)imageSize);
    cudaMalloc( (void**) &d_G2x, sizeof(float)*(unsigned)imageSize);
    cudaMalloc( (void**) &d_G2y, sizeof(float)*(unsigned)imageSize);
    cudaMalloc( (void**) &d_G2z, sizeof(float)*(unsigned)imageSize);
    cudaMalloc( (void**) &d_G2t, sizeof(float)*(unsigned)imageSize);
    cudaMalloc( (void**) &d_G2f, sizeof(float)*(unsigned)imageSize);
    
    cudaMalloc( (void**) &d_FPS, sizeof(float)*(unsigned)imageSize);
    cudaMalloc( (void**) &d_qq1, sizeof(float)*(unsigned)imageSize);
    cudaMalloc( (void**) &d_qq2, sizeof(float)*(unsigned)imageSize);
    cudaMalloc( (void**) &d_qq3, sizeof(float)*(unsigned)imageSize);
    
    /* copy arrays from host to device */
    cudaMemcpy( d_b1x1, h_b1x1, sizeof(float)*(unsigned)imageSize, cudaMemcpyHostToDevice);
    cudaMemcpy( d_b1x2, h_b1x2, sizeof(float)*(unsigned)imageSize, cudaMemcpyHostToDevice);
    cudaMemcpy( d_b1x3, h_b1x3, sizeof(float)*(unsigned)imageSize, cudaMemcpyHostToDevice);
    cudaMemcpy( d_b1y1, h_b1y1, sizeof(float)*(unsigned)imageSize, cudaMemcpyHostToDevice);
    cudaMemcpy( d_b1y2, h_b1y2, sizeof(float)*(unsigned)imageSize, cudaMemcpyHostToDevice);
    cudaMemcpy( d_b1y3, h_b1y3, sizeof(float)*(unsigned)imageSize, cudaMemcpyHostToDevice);
    cudaMemcpy( d_b1z1, h_b1z1, sizeof(float)*(unsigned)imageSize, cudaMemcpyHostToDevice);
    cudaMemcpy( d_b1z2, h_b1z2, sizeof(float)*(unsigned)imageSize, cudaMemcpyHostToDevice);
    cudaMemcpy( d_b1z3, h_b1z3, sizeof(float)*(unsigned)imageSize, cudaMemcpyHostToDevice);
    cudaMemcpy( d_gkx, h_gkx, sizeof(float)*(unsigned)imageSize, cudaMemcpyHostToDevice);
    cudaMemcpy( d_gky, h_gky, sizeof(float)*(unsigned)imageSize, cudaMemcpyHostToDevice);
    cudaMemcpy( d_gkz, h_gkz, sizeof(float)*(unsigned)imageSize, cudaMemcpyHostToDevice);
    cudaMemcpy( d_dv1x, h_dv1x, sizeof(float)*(unsigned)imageSize, cudaMemcpyHostToDevice);
    cudaMemcpy( d_dv1y, h_dv1y, sizeof(float)*(unsigned)imageSize, cudaMemcpyHostToDevice);
    cudaMemcpy( d_dv1z, h_dv1z, sizeof(float)*(unsigned)imageSize, cudaMemcpyHostToDevice);
    cudaMemcpy( d_q1,  h_q1,  sizeof(float)*(unsigned)imageSize, cudaMemcpyHostToDevice);
    cudaMemcpy( d_u1x, h_u1x, sizeof(float)*(unsigned)imageSize, cudaMemcpyHostToDevice);
    cudaMemcpy( d_u1y, h_u1y, sizeof(float)*(unsigned)imageSize, cudaMemcpyHostToDevice);
    cudaMemcpy( d_u1z ,h_u1z, sizeof(float)*(unsigned)imageSize, cudaMemcpyHostToDevice);
    cudaMemcpy( d_U1x, h_U1x, sizeof(float)*(unsigned)imageSize, cudaMemcpyHostToDevice);
    cudaMemcpy( d_U1y, h_U1y, sizeof(float)*(unsigned)imageSize, cudaMemcpyHostToDevice);
    cudaMemcpy( d_U1z, h_U1z, sizeof(float)*(unsigned)imageSize, cudaMemcpyHostToDevice);
    cudaMemcpy( d_G1x, h_G1x, sizeof(float)*(unsigned)imageSize, cudaMemcpyHostToDevice);
    cudaMemcpy( d_G1y, h_G1y, sizeof(float)*(unsigned)imageSize, cudaMemcpyHostToDevice);
    cudaMemcpy( d_G1z, h_G1z, sizeof(float)*(unsigned)imageSize, cudaMemcpyHostToDevice);
    cudaMemcpy( d_G1t, h_G1t, sizeof(float)*(unsigned)imageSize, cudaMemcpyHostToDevice);
    cudaMemcpy( d_G1f, h_G1f, sizeof(float)*(unsigned)imageSize, cudaMemcpyHostToDevice);

    cudaMemcpy( d_b2x1, h_b2x1, sizeof(float)*(unsigned)imageSize, cudaMemcpyHostToDevice);
    cudaMemcpy( d_b2x2, h_b2x2, sizeof(float)*(unsigned)imageSize, cudaMemcpyHostToDevice);
    cudaMemcpy( d_b2x3, h_b2x3, sizeof(float)*(unsigned)imageSize, cudaMemcpyHostToDevice);
    cudaMemcpy( d_b2y1, h_b2y1, sizeof(float)*(unsigned)imageSize, cudaMemcpyHostToDevice);
    cudaMemcpy( d_b2y2, h_b2y2, sizeof(float)*(unsigned)imageSize, cudaMemcpyHostToDevice);
    cudaMemcpy( d_b2y3, h_b2y3, sizeof(float)*(unsigned)imageSize, cudaMemcpyHostToDevice);
    cudaMemcpy( d_b2z1, h_b2z1, sizeof(float)*(unsigned)imageSize, cudaMemcpyHostToDevice);
    cudaMemcpy( d_b2z2, h_b2z2, sizeof(float)*(unsigned)imageSize, cudaMemcpyHostToDevice);
    cudaMemcpy( d_b2z3, h_b2z3, sizeof(float)*(unsigned)imageSize, cudaMemcpyHostToDevice);
    cudaMemcpy( d_dv2x, h_dv2x, sizeof(float)*(unsigned)imageSize, cudaMemcpyHostToDevice);
    cudaMemcpy( d_dv2y, h_dv2y, sizeof(float)*(unsigned)imageSize, cudaMemcpyHostToDevice);
    cudaMemcpy( d_dv2z, h_dv2z, sizeof(float)*(unsigned)imageSize, cudaMemcpyHostToDevice);
    cudaMemcpy( d_q2,  h_q2,  sizeof(float)*(unsigned)imageSize, cudaMemcpyHostToDevice);
    cudaMemcpy( d_u2x, h_u2x, sizeof(float)*(unsigned)imageSize, cudaMemcpyHostToDevice);
    cudaMemcpy( d_u2y, h_u2y, sizeof(float)*(unsigned)imageSize, cudaMemcpyHostToDevice);
    cudaMemcpy( d_u2z ,h_u2z, sizeof(float)*(unsigned)imageSize, cudaMemcpyHostToDevice);
    cudaMemcpy( d_U2x, h_U2x, sizeof(float)*(unsigned)imageSize, cudaMemcpyHostToDevice);
    cudaMemcpy( d_U2y, h_U2y, sizeof(float)*(unsigned)imageSize, cudaMemcpyHostToDevice);
    cudaMemcpy( d_U2z, h_U2z, sizeof(float)*(unsigned)imageSize, cudaMemcpyHostToDevice);
    cudaMemcpy( d_G2x, h_G2x, sizeof(float)*(unsigned)imageSize, cudaMemcpyHostToDevice);
    cudaMemcpy( d_G2y, h_G2y, sizeof(float)*(unsigned)imageSize, cudaMemcpyHostToDevice);
    cudaMemcpy( d_G2z, h_G2z, sizeof(float)*(unsigned)imageSize, cudaMemcpyHostToDevice);
    cudaMemcpy( d_G2t, h_G2t, sizeof(float)*(unsigned)imageSize, cudaMemcpyHostToDevice);
    cudaMemcpy( d_G2f, h_G2f, sizeof(float)*(unsigned)imageSize, cudaMemcpyHostToDevice);
    
    cudaMemcpy( d_qq1, h_qq1, sizeof(float)*(unsigned)imageSize, cudaMemcpyHostToDevice);
    cudaMemcpy( d_qq2, h_qq2, sizeof(float)*(unsigned)imageSize, cudaMemcpyHostToDevice);
    cudaMemcpy( d_qq3, h_qq3, sizeof(float)*(unsigned)imageSize, cudaMemcpyHostToDevice);
    
    
    /* run optimization */
    
    /* iNI = 0; */
    dim3 threads(BLOCKSIZE,1,1);
    dim3 grid = GetGrid(imageSize);
    
    for( int i = 0; i < maxIter; i++)
    {
        
        /* update p1 */
        krnl_1<<<grid, threads>>>(d_dv1x, d_dv1y, d_dv1z, 
                d_u1x, d_u1y, d_u1z,
                d_gkx, d_gky, d_gkz, 
                d_G1x, d_G1y, d_G1z,
                d_G1t, d_G1f, d_q1, 
                d_U1x, d_U1y, d_U1z,
                cc, iNx, iNy, iNz, d_qq1, d_qq2, d_qq3, flag, weight1);

        /* update p1x p1y, p1z */
        krnl_23z<<<grid, threads>>>(d_b1x1, d_b1y1, d_b1z1, 
                d_b1x2, d_b1y2, d_b1z2,
                d_b1x3, d_b1y3, d_b1z3,
                d_gkx, d_gky, d_gkz,
                steps, iNx, iNy, iNz);
        /*
//         krnl_2<<<grid, threads>>>(d_bx1, d_by1, d_bz1, d_gkx, d_gky, d_gkz,
//                 steps, iNx, iNy, iNz);
// 
//         krnl_3<<<grid, threads>>>(d_bx2, d_by2, d_bz2, d_gkx, d_gky, d_gkz,
//                 steps, iNx, iNy, iNz);
//         
//         krnl_z<<<grid, threads>>>(d_bx3, d_by3, d_bz3, d_gkx, d_gky, d_gkz,
//                 steps, iNx, iNy, iNz);
        */
        
        /* projection step */
        krnl_4<<<grid, threads>>>(d_b1x1, d_b1x2, d_b1x3, 
                d_b1y1, d_b1y2, d_b1y3,
                d_b1z1, d_b1z2, d_b1z3, 
                d_gkx, d_gky, d_gkz,
                fPenalty1*weight1, iNx, iNy, iNz);
        
        krnl_56zp<<<grid, threads>>>(d_b1x1, d_b1y1, d_b1z1,
                d_b1x2, d_b1y2, d_b1z2,
                d_b1x3, d_b1y3, d_b1z3,
                d_gkx, d_gky, d_gkz,
                iNx, iNy, iNz);
        /*
//         krnl_5<<<grid, threads>>>(d_bx1, d_by1, d_bz1, d_gkx, d_gky, d_gkz,
//                 iNx, iNy, iNz);
//         
//         krnl_6<<<grid, threads>>>(d_bx2, d_by2, d_bz2, d_gkx, d_gky, d_gkz,
//                 iNx, iNy, iNz);
//         
//         krnl_zp<<<grid, threads>>>(d_bx3, d_by3, d_bz3, d_gkx, d_gky, d_gkz,
//                 iNx, iNy, iNz);
        */
        krnl_7<<<grid, threads>>>(d_b1x1, d_b1x2, d_b1x3, 
                d_b1y1, d_b1y2, d_b1y3,
                d_b1z1, d_b1z2, d_b1z3, 
                d_dv1x, d_dv1y, d_dv1z, 
                iNx, iNy, iNz);
        
        /* update p2 */
        krnl_1<<<grid, threads>>>(d_dv2x, d_dv2y, d_dv2z, 
                d_u2x, d_u2y, d_u2z,
                d_gkx, d_gky, d_gkz, 
                d_G2x, d_G2y, d_G2z,
                d_G2t, d_G2f, d_q2, 
                d_U2x, d_U2y, d_U2z,
                cc, iNx, iNy, iNz, d_qq1, d_qq2, d_qq3, -flag, weight2);

        /* update p2x, p2y, p2z */
        krnl_23z<<<grid, threads>>>(d_b2x1, d_b2y1, d_b2z1, 
                d_b2x2, d_b2y2, d_b2z2,
                d_b2x3, d_b2y3, d_b2z3,
                d_gkx, d_gky, d_gkz,
                steps, iNx, iNy, iNz);
        /*
//         krnl_2<<<grid, threads>>>(d_bx1, d_by1, d_bz1, d_gkx, d_gky, d_gkz,
//                 steps, iNx, iNy, iNz);
// 
//         krnl_3<<<grid, threads>>>(d_bx2, d_by2, d_bz2, d_gkx, d_gky, d_gkz,
//                 steps, iNx, iNy, iNz);
//         
//         krnl_z<<<grid, threads>>>(d_bx3, d_by3, d_bz3, d_gkx, d_gky, d_gkz,
//                 steps, iNx, iNy, iNz);
        */
        
        /* projection step */
        krnl_4<<<grid, threads>>>(d_b2x1, d_b2x2, d_b2x3, 
                d_b2y1, d_b2y2, d_b2y3,
                d_b2z1, d_b2z2, d_b2z3, 
                d_gkx, d_gky, d_gkz,
                fPenalty2*weight2, iNx, iNy, iNz);
        
        krnl_56zp<<<grid, threads>>>(d_b2x1, d_b2y1, d_b2z1,
                d_b2x2, d_b2y2, d_b2z2,
                d_b2x3, d_b2y3, d_b2z3,
                d_gkx, d_gky, d_gkz,
                iNx, iNy, iNz);
        /*
//         krnl_5<<<grid, threads>>>(d_bx1, d_by1, d_bz1, d_gkx, d_gky, d_gkz,
//                 iNx, iNy, iNz);
//         
//         krnl_6<<<grid, threads>>>(d_bx2, d_by2, d_bz2, d_gkx, d_gky, d_gkz,
//                 iNx, iNy, iNz);
//         
//         krnl_zp<<<grid, threads>>>(d_bx3, d_by3, d_bz3, d_gkx, d_gky, d_gkz,
//                 iNx, iNy, iNz);
        */
        krnl_7<<<grid, threads>>>(d_b2x1, d_b2x2, d_b2x3, 
                d_b2y1, d_b2y2, d_b2y3,
                d_b2z1, d_b2z2, d_b2z3, 
                d_dv2x, d_dv2y, d_dv2z, 
                iNx, iNy, iNz);
        
       /* compute qq */
       krnl_8<<<grid, threads>>>(d_U1x, d_U1y, d_U1z, d_U2x, d_U2y, d_U2z,
                                 d_q1, d_q2,
                                 d_G1x, d_G1y, d_G1z, d_G2x, d_G2y, d_G2z,
                                 d_dv1x, d_dv1y, d_dv1z, d_dv2x, d_dv2y, d_dv2z, 
                                 d_u1x, d_u1y, d_u1z, d_u2x, d_u2y, d_u2z,
                                 cc, fPenalty3, iNx, iNy, iNz, d_qq1, d_qq2, d_qq3);
       /*compute u1,2x, u1,2y, u1,2z*/
       krnl_9<<<grid, threads>>>(d_dv1x, d_dv1y, d_dv1z, d_dv2x, d_dv2y, d_dv2z, 
                                 d_G1x, d_G1y, d_G1z, d_G2x, d_G2y, d_G2z,
                                 d_q1, d_q2, d_u1x, d_u1y, d_u1z, d_u2x, d_u2y, d_u2z,
                                 d_FPS, cc, iNx, iNy, iNz, d_qq1, d_qq2, d_qq3);
       
        
        /* compute convergence */
        cudaMemcpy( h_FPS, d_FPS, sizeof(float)*unsigned(imageSize), cudaMemcpyDeviceToHost);
        
        fps = 0;
        for (int j=0; j< imageSize; j++){
            fps += abs(h_FPS[j]);
        }
        
        h_cvg[i] = fps / (float)imageSize /6;
        
        if (h_cvg[i] <= fError){
            break; 
        }
        
        /*mexPrintf("cvg: %f\n",h_cvg[i]); */
        
        punum[0] = i+1;
        
    }
    
    /* copy arrays from device to host */
    cudaMemcpy( h_u1x, d_u1x, sizeof(float)*(unsigned)(imageSize), cudaMemcpyDeviceToHost);
    cudaMemcpy( h_u1y, d_u1y, sizeof(float)*(unsigned)(imageSize), cudaMemcpyDeviceToHost);
    cudaMemcpy( h_u1z, d_u1z, sizeof(float)*(unsigned)(imageSize), cudaMemcpyDeviceToHost);
    cudaMemcpy( h_u2x, d_u2x, sizeof(float)*(unsigned)(imageSize), cudaMemcpyDeviceToHost);
    cudaMemcpy( h_u2y, d_u2y, sizeof(float)*(unsigned)(imageSize), cudaMemcpyDeviceToHost);
    cudaMemcpy( h_u2z, d_u2z, sizeof(float)*(unsigned)(imageSize), cudaMemcpyDeviceToHost);
    
    mexPrintf("number of iterations = %i\n",punum[0]);
    
    
    /* Free memory */
    free( (float *) h_b1x1 );
    free( (float *) h_b1x2 );
    free( (float *) h_b1x3 );
    free( (float *) h_b1y1 );
    free( (float *) h_b1y2 );
    free( (float *) h_b1y3 );
    free( (float *) h_b1z1 );
    free( (float *) h_b1z2 );
    free( (float *) h_b1z3 );
    free( (float *) h_gkx );
    free( (float *) h_gky );
    free( (float *) h_gkz );
    free( (float *) h_dv1x );
    free( (float *) h_dv1y );
    free( (float *) h_dv1z );
    free( (float *) h_q1 );

    free( (float *) h_b2x1 );
    free( (float *) h_b2x2 );
    free( (float *) h_b2x3 );
    free( (float *) h_b2y1 );
    free( (float *) h_b2y2 );
    free( (float *) h_b2y3 );
    free( (float *) h_b2z1 );
    free( (float *) h_b2z2 );
    free( (float *) h_b2z3 );
    free( (float *) h_dv2x );
    free( (float *) h_dv2y );
    free( (float *) h_dv2z );
    free( (float *) h_q2 );
    
    free( (float *) h_FPS );
    free( (float *) h_qq1 );
    free( (float *) h_qq2 );
    free( (float *) h_qq3 );
    
    /*    Free GPU Memory */
    cudaFree(d_b1x1);
    cudaFree(d_b1x2);
    cudaFree(d_b1x3);
    cudaFree(d_b1y1);
    cudaFree(d_b1y2);
    cudaFree(d_b1y3);
    cudaFree(d_b1z1);
    cudaFree(d_b1z2);
    cudaFree(d_b1z3);
    cudaFree(d_gkx);
    cudaFree(d_gky);
    cudaFree(d_gkz);
    cudaFree(d_dv1x);
    cudaFree(d_dv1y);
    cudaFree(d_dv1z);
    
    cudaFree(d_u1x);
    cudaFree(d_u1y);
    cudaFree(d_u1z);
    cudaFree(d_U1x);
    cudaFree(d_U1y);
    cudaFree(d_U1z);
    cudaFree(d_G1x);
    cudaFree(d_G1y);
    cudaFree(d_G1z);
    cudaFree(d_G1t);
    cudaFree(d_G1f);
    cudaFree(d_q1);

    cudaFree(d_b2x1);
    cudaFree(d_b2x2);
    cudaFree(d_b2x3);
    cudaFree(d_b2y1);
    cudaFree(d_b2y2);
    cudaFree(d_b2y3);
    cudaFree(d_b2z1);
    cudaFree(d_b2z2);
    cudaFree(d_b2z3);
    cudaFree(d_dv2x);
    cudaFree(d_dv2y);
    cudaFree(d_dv2z);
    
    cudaFree(d_u2x);
    cudaFree(d_u2y);
    cudaFree(d_u2z);
    cudaFree(d_U2x);
    cudaFree(d_U2y);
    cudaFree(d_U2z);
    cudaFree(d_G2x);
    cudaFree(d_G2y);
    cudaFree(d_G2z);
    cudaFree(d_G2t);
    cudaFree(d_G2f);
    cudaFree(d_q2);

    cudaFree(d_FPS);
    cudaFree(d_qq1);
    cudaFree(d_qq2);
    cudaFree(d_qq3);

    /* CUDA event-based timer */
    cudaEventRecord( stop, 0 );
    cudaEventSynchronize( stop );
    cudaEventElapsedTime( &time, start, stop );
    cudaEventDestroy( start );
    cudaEventDestroy( stop );
    
    
    tt[0] = time;
    
    mexPrintf("\nComputational Time for Dual Optimization = %.4f sec\n \n",tt[0]/1000000);
    
    
}
