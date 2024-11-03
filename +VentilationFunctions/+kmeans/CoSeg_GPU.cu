#include <stdio.h>
#include <stdlib.h>
#include <mex.h>
#include <math.h>
#include <time.h>
#include "cuda.h"
#include "cuda_runtime_api.h"
#include "CoSeg_kernels.cu"

#define YES 0
#define NO 1

#define PI 3.1415926

#define MAX(a,b) ( a > b ? a : b )
#define MIN(a,b) ( a <= b ? a : b )
#define SIGN(x) ( x >= 0.0 ? 1.0 : -1.0 )
#define ABS(x) ( (x) > 0.0 ? x : -(x) )

#ifndef HAVE_RINT 
#define rint(A) floor((A)+(((A) < 0)? -0.5 : 0.5)) 
#endif


/**********************************************/
/************** MAIN FUNCTION *****************/
/**********************************************/

/****************************************/
extern void mexFunction(int iNbOut, mxArray *pmxOut[],
int iNbIn, const mxArray *pmxIn[])
{
    
  /* iNbOut: number of outputs
     pmxOut: array of pointers to output arguments */
    
  /* iNbIn: number of inputs
     pmxIn: array of pointers to input arguments */
    
    
    float   *pfpenalty1, *pfpenalty2, *pfu1, *pfu2, *pfCs1, *pfCs2;
    float   *pfCt1, *pfCt2, *pfqq, *pfcvg, *pfVecParameters;
    float   *pfbx1, *pfbx2, *pfby1, *pfby2, *pfbz1, *pfbz2, *pfps1, *pfps2;
    float   *pfpt1, *pfpt2, *pfgk1, *pfgk2, *tt, *pfdv1, *pfdv2;

    float   fError, cc, steps, fps, beta;
    int     *punum, iNy, iNx, iNz, iNdim, iDim[3], ix, iy, iNI;
    int     SZF, idz, iz;

    int iDev;
    int     iNbIters, szImg, idx, index;
    time_t  start_time, end_time;

    //    GPU Variables
    float   *pfbx1_GPU, *pfby1_GPU, *pfpenalty1_GPU, *pfdv1_GPU;
    float   *pfbx2_GPU, *pfby2_GPU, *pfpenalty2_GPU, *pfdv2_GPU;
    float   *pfbz1_GPU, *pfbz2_GPU;
    float   *pfps1_GPU, *pfpt1_GPU, *pfgk1_GPU, *pfu1_GPU, *pfCs1_GPU, *pfCt1_GPU;
    float   *pfps2_GPU, *pfpt2_GPU, *pfgk2_GPU, *pfu2_GPU, *pfCs2_GPU, *pfCt2_GPU;
    float   *FPS, *FPS_GPU, *pfqq_GPU;
    
    
    cudaDeviceProp prop;

    cudaGetDeviceCount(&iDev);

    if ((unsigned int)iDev == 0){
        printf("There is no CUDA device found!");
        return;
    }
    else{
        printf("There are %d CUDA devices in your computer. \n", iDev);
        for(int ii = 0; ii < iDev; ii ++){
            cudaGetDeviceProperties(&prop, ii);
            printf("------ General Information for CUDA device %d ------ \n", ii);
            printf("Name:  %s \n", prop.name);
            printf("Multiprocessor count:  %d \n", prop.multiProcessorCount);
            printf("Total global memory: %ld \n", prop.totalGlobalMem);
            printf("---------------------------------------------------- \n\n");
         }
    }

    /* Inputs */
    pfpenalty1 = (float*)mxGetData(pmxIn[0]); /* Given penalty1,2 */
    pfpenalty2 = (float*)mxGetData(pmxIn[1]);
    pfCs1 = (float*)mxGetData(pmxIn[2]); /* bound of source flows ps1 and ps2 */
    pfCs2 = (float*)mxGetData(pmxIn[3]);
    pfCt1 = (float*)mxGetData(pmxIn[4]); /* bound of sink flows pt1 and pt2 */
    pfCt2 = (float*)mxGetData(pmxIn[5]);
    pfVecParameters = (float*)mxGetData(pmxIn[6]); /* Vector of parameters */

     /* 
     *pfVecParameters Setting
     * [0] : number of columns 
     * [1] : number of rows
     * [2] : the maximum iteration number
     * [3] : error criterion
     * [4] : cc for the step-size of ALM
     * [5] : steps for the step-size of projected-gradient of p
     */
   
    /* Size */
    iNy = (int) pfVecParameters[0];
    iNx = (int) pfVecParameters[1];
    iNz = (int) pfVecParameters[2];
    szImg = iNy*iNx*iNz;
    SZF = iNy*iNx;

    /* Choice of region segmentation model */
    iNbIters = (int) pfVecParameters[3]; /* the maximum iteration number */
    fError = (float) pfVecParameters[4]; /* error bound for convergence */
    cc = (float) pfVecParameters[5]; /* the step-size for ALM */
    steps = (float) pfVecParameters[6]; /* the step-size for each projected-gradient step */
    beta = (float) pfVecParameters[7];

    printf("Initializing ................................................ \n\n");

    /* Outputs */
    /* outputs the computed u1(x)  */
    iNdim = 3;
    iDim[0] = iNy;
    iDim[1] = iNx;
    iDim[2] = iNz;
    
    pmxOut[0] = mxCreateNumericArray(iNdim,(const int*)iDim,mxSINGLE_CLASS,mxREAL);
    pfu1 = (float*)mxGetData(pmxOut[0]);
    
    /* outputs the computed u2(x)  */
    iNdim = 3;
    iDim[0] = iNy;
    iDim[1] = iNx;
    iDim[2] = iNz;

    pmxOut[1] = mxCreateNumericArray(iNdim,(const int*)iDim,mxSINGLE_CLASS,mxREAL);
    pfu2 = (float*)mxGetData(pmxOut[1]);
   
    /* outputs the convergence rate  */
    iNdim = 2;
    iDim[0] = 1;
    iDim[1] = iNbIters;
    
    pmxOut[2] = mxCreateNumericArray(iNdim,(const int*)iDim,mxSINGLE_CLASS,mxREAL);
    pfcvg = (float*)mxGetData(pmxOut[2]);
    
    /* outputs the iteration number  */
    iNdim = 2;
    iDim[0] = 1;
    iDim[1] = 1;
    pmxOut[3] = mxCreateNumericArray(iNdim,(const int*)iDim,mxUINT16_CLASS,mxREAL);
    punum = (int*)mxGetData(pmxOut[3]);
    
    /* outputs the computation time  */
    iNdim = 2;
    iDim[0] = 1;
    iDim[1] = 1;
    pmxOut[4] = mxCreateNumericArray(iNdim,(const int*)iDim,mxSINGLE_CLASS,mxREAL);
    tt = (float*)mxGetData(pmxOut[4]);
    
    /* Memory allocation */
    
    /* allocate the memory for px1 and px2 */
    pfbx1 = (float *) calloc( (unsigned)(iNy*(iNx+1)*iNz), sizeof(float) );
    if (!pfbx1)
        mexPrintf("Memory allocation failure\n");
    
    pfbx2 = (float *) calloc( (unsigned)(iNy*(iNx+1)*iNz), sizeof(float) );
    if (!pfbx2)
        mexPrintf("Memory allocation failure\n");
    
    /* allocate the memory for py1 and py2 */
    pfby1 = (float *) calloc( (unsigned)((iNy+1)*iNx*iNz), sizeof(float) );
    if (!pfby1)
        mexPrintf("Memory allocation failure\n");
    
    pfby2 = (float *) calloc( (unsigned)((iNy+1)*iNx*iNz), sizeof(float) );
    if (!pfby2)
        mexPrintf("Memory allocation failure\n");
    
    /* allocate the memory for pz1 and pz2 */
    pfbz1 = (float *) calloc( (unsigned)(iNy*iNx*(iNz+1)), sizeof(float) );
    if (!pfbz1)
        mexPrintf("Memory allocation failure\n");

    pfbz2 = (float *) calloc( (unsigned)(iNy*iNx*(iNz+1)), sizeof(float) );
    if (!pfbz2)
        mexPrintf("Memory allocation failure\n");

    /* allocate the memory for ps1 and ps2 */
    pfps1 = (float *) calloc( (unsigned)(iNy*iNx*iNz), sizeof(float) );
    if (!pfps1)
        mexPrintf("Memory allocation failure\n");
    
    pfps2 = (float *) calloc( (unsigned)(iNy*iNx*iNz), sizeof(float) );
    if (!pfps2)
        mexPrintf("Memory allocation failure\n");
    
    /* allocate the memory for pt1 and pt2 */
    pfpt1 = (float *) calloc( (unsigned)(iNy*iNx*iNz), sizeof(float) );
    if (!pfpt1)
        mexPrintf("Memory allocation failure\n");
    
    pfpt2 = (float *) calloc( (unsigned)(iNy*iNx*iNz), sizeof(float) );
    if (!pfpt2)
        mexPrintf("Memory allocation failure\n");
    
    /* allocate the memory for the coupled flow q */
    pfqq = (float *) calloc( (unsigned)(iNy*iNx*iNz), sizeof(float) );
    if (!pfqq)
        mexPrintf("Memory allocation failure\n");
    
    /* allocate the memory for gk1 */
    pfgk1 = (float *) calloc( (unsigned)(iNy*iNx*iNz), sizeof(float) );
    if (!pfgk1)
        mexPrintf("Memory allocation failure\n");
    
    /* allocate the memory for gk2 */
    pfgk2 = (float *) calloc( (unsigned)(iNy*iNx*iNz), sizeof(float) );
    if (!pfgk2)
        mexPrintf("Memory allocation failure\n");
    
    /* allocate the memory for div1 and div2 */
    pfdv1 = (float *) calloc( (unsigned)(iNy*iNx*iNz), sizeof(float) );
    if (!pfdv1)
        mexPrintf("Memory allocation failure\n");
    
    pfdv2 = (float *) calloc( (unsigned)(iNy*iNx*iNz), sizeof(float) );
    if (!pfdv2)
        mexPrintf("Memory allocation failure\n");

    /* allocate the memory for FPS */
    FPS = (float *) calloc( (unsigned)(iNy*iNx*iNz), sizeof(float) );
    if (!FPS)
        mexPrintf("Memory allocation failure\n");


     //    GPU Memory Allocation
    
    cudaMalloc( (void**) &pfbx1_GPU, sizeof(float)*(unsigned)((iNx+1)*iNy*iNz));
    cudaMalloc( (void**) &pfby1_GPU, sizeof(float)*(unsigned)(iNx*(iNy+1)*iNz));
    cudaMalloc( (void**) &pfbz1_GPU, sizeof(float)*(unsigned)(iNx*iNy*(iNz+1)));
    cudaMalloc( (void**) &pfpenalty1_GPU, sizeof(float)*(unsigned)(iNy*iNx*iNz));
    cudaMalloc( (void**) &pfdv1_GPU, sizeof(float)*(unsigned)(iNy*iNx*iNz));
    cudaMalloc( (void**) &pfps1_GPU, sizeof(float)*(unsigned)(iNy*iNx*iNz));
    cudaMalloc( (void**) &pfpt1_GPU, sizeof(float)*(unsigned)(iNy*iNx*iNz));
    cudaMalloc( (void**) &pfgk1_GPU, sizeof(float)*(unsigned)(iNy*iNx*iNz));
    cudaMalloc( (void**) &pfu1_GPU, sizeof(float)*(unsigned)(iNy*iNx*iNz));
    cudaMalloc( (void**) &pfCs1_GPU, sizeof(float)*(unsigned)(iNy*iNx*iNz));
    cudaMalloc( (void**) &pfCt1_GPU, sizeof(float)*(unsigned)(iNy*iNx*iNz));

    cudaMalloc( (void**) &pfbx2_GPU, sizeof(float)*(unsigned)((iNx+1)*iNy*iNz));
    cudaMalloc( (void**) &pfby2_GPU, sizeof(float)*(unsigned)(iNx*(iNy+1)*iNz));
    cudaMalloc( (void**) &pfbz2_GPU, sizeof(float)*(unsigned)(iNx*iNy*(iNz+1)));
    cudaMalloc( (void**) &pfpenalty2_GPU, sizeof(float)*(unsigned)(iNy*iNx*iNz));
    cudaMalloc( (void**) &pfdv2_GPU, sizeof(float)*(unsigned)(iNy*iNx*iNz));
    cudaMalloc( (void**) &pfps2_GPU, sizeof(float)*(unsigned)(iNy*iNx*iNz));
    cudaMalloc( (void**) &pfpt2_GPU, sizeof(float)*(unsigned)(iNy*iNx*iNz));
    cudaMalloc( (void**) &pfgk2_GPU, sizeof(float)*(unsigned)(iNy*iNx*iNz));
    cudaMalloc( (void**) &pfu2_GPU, sizeof(float)*(unsigned)(iNy*iNx*iNz));
    cudaMalloc( (void**) &pfCs2_GPU, sizeof(float)*(unsigned)(iNy*iNx*iNz));
    cudaMalloc( (void**) &pfCt2_GPU, sizeof(float)*(unsigned)(iNy*iNx*iNz));

    cudaMalloc( (void**) &pfqq_GPU, sizeof(float)*(unsigned)(iNy*iNx*iNz));
    cudaMalloc( (void**) &FPS_GPU, sizeof(float)*(unsigned)(iNy*iNx*iNz));

    /* Preprocessing initial values */
    for (iz=0; iz < iNz; iz++){
        idz = iz*SZF;
        for (ix=0; ix< iNx; ix++){
            idx = idz + ix*iNy;
            for (iy=0; iy< iNy; iy++){
               index = idx + iy; 

                if (pfCs1[index] < pfCt1[index]){
                    pfps1[index] = pfCs1[index];
                    pfpt1[index] = pfCs1[index];
                    pfdv1[index] = pfbx1[index+iNy] - pfbx1[index] 
                                + pfby1[index+1] - pfby1[index]
                                + pfbz1[index+SZF] - pfbz1[index];
                }
                else{
                    pfu1[index] = 1;
                    pfps1[index] = pfCt1[index];
                    pfpt1[index] = pfCt1[index];
                    pfdv1[index] = pfbx1[index+iNy] - pfbx1[index] 
                                + pfby1[index+1] - pfby1[index]
                                + pfbz1[index+SZF] - pfbz1[index];
                }

               if (pfCs2[index] < pfCt2[index]){
                    pfps2[index] = pfCs2[index];
                    pfpt2[index] = pfCs2[index];
                    pfdv2[index] = pfbx2[index+iNy] - pfbx2[index] 
                                 + pfby2[index+1] - pfby2[index]
                                 + pfbz2[index+SZF] - pfbz2[index];
                }
                else{
                    pfu2[index] = 1;
                    pfps2[index] = pfCt2[index];
                    pfpt2[index] = pfCt2[index];
                    pfdv2[index] = pfbx2[index+iNy] - pfbx2[index] 
                                + pfby2[index+1] - pfby2[index]
                                + pfbz2[index+SZF] - pfbz2[index];
                }
            }
        }
    }
    
    //    Copy Parameters from Host to Device

    cudaMemcpy( pfbx1_GPU, pfbx1, sizeof(float)*(unsigned)(iNy*(iNx+1)*iNz), cudaMemcpyHostToDevice);
    cudaMemcpy( pfby1_GPU, pfby1, sizeof(float)*(unsigned)((iNy+1)*iNx*iNz), cudaMemcpyHostToDevice);
    cudaMemcpy( pfbz1_GPU, pfbz1, sizeof(float)*(unsigned)(iNy*iNx*(iNz+1)), cudaMemcpyHostToDevice);
    cudaMemcpy( pfpenalty1_GPU, pfpenalty1, sizeof(float)*(unsigned)(iNy*iNx*iNz), cudaMemcpyHostToDevice);
    cudaMemcpy( pfdv1_GPU, pfdv1, sizeof(float)*(unsigned)(iNy*iNx*iNz), cudaMemcpyHostToDevice);
    cudaMemcpy( pfps1_GPU, pfps1, sizeof(float)*(unsigned)(iNy*iNx*iNz), cudaMemcpyHostToDevice);
    cudaMemcpy( pfpt1_GPU, pfpt1, sizeof(float)*(unsigned)(iNy*iNx*iNz), cudaMemcpyHostToDevice);
    cudaMemcpy( pfgk1_GPU, pfgk1, sizeof(float)*(unsigned)(iNy*iNx*iNz), cudaMemcpyHostToDevice);
    cudaMemcpy( pfu1_GPU, pfu1, sizeof(float)*(unsigned)(iNy*iNx*iNz), cudaMemcpyHostToDevice);
    cudaMemcpy( pfCs1_GPU, pfCs1, sizeof(float)*(unsigned)(iNy*iNx*iNz), cudaMemcpyHostToDevice);
    cudaMemcpy( pfCt1_GPU, pfCt1, sizeof(float)*(unsigned)(iNy*iNx*iNz), cudaMemcpyHostToDevice);

    cudaMemcpy( pfbx2_GPU, pfbx2, sizeof(float)*(unsigned)(iNy*(iNx+1)*iNz), cudaMemcpyHostToDevice);
    cudaMemcpy( pfby2_GPU, pfby2, sizeof(float)*(unsigned)((iNy+1)*iNx*iNz), cudaMemcpyHostToDevice);
    cudaMemcpy( pfbz2_GPU, pfbz2, sizeof(float)*(unsigned)(iNy*iNx*(iNz+1)), cudaMemcpyHostToDevice);
    cudaMemcpy( pfpenalty2_GPU, pfpenalty2, sizeof(float)*(unsigned)(iNy*iNx*iNz), cudaMemcpyHostToDevice);
    cudaMemcpy( pfdv2_GPU, pfdv2, sizeof(float)*(unsigned)(iNy*iNx*iNz), cudaMemcpyHostToDevice);
    cudaMemcpy( pfps2_GPU, pfps2, sizeof(float)*(unsigned)(iNy*iNx*iNz), cudaMemcpyHostToDevice);
    cudaMemcpy( pfpt2_GPU, pfpt2, sizeof(float)*(unsigned)(iNy*iNx*iNz), cudaMemcpyHostToDevice);
    cudaMemcpy( pfgk2_GPU, pfgk2, sizeof(float)*(unsigned)(iNy*iNx*iNz), cudaMemcpyHostToDevice);
    cudaMemcpy( pfu2_GPU, pfu2, sizeof(float)*(unsigned)(iNy*iNx*iNz), cudaMemcpyHostToDevice);
    cudaMemcpy( pfCs2_GPU, pfCs2, sizeof(float)*(unsigned)(iNy*iNx*iNz), cudaMemcpyHostToDevice);
    cudaMemcpy( pfCt2_GPU, pfCt2, sizeof(float)*(unsigned)(iNy*iNx*iNz), cudaMemcpyHostToDevice);
    
    cudaMemcpy( pfqq_GPU, pfqq, sizeof(float)*(unsigned)(iNy*iNx*iNz), cudaMemcpyHostToDevice);
    cudaMemcpy( FPS_GPU, FPS, sizeof(float)*(unsigned)(iNy*iNx*iNz), cudaMemcpyHostToDevice);

    /*  Main iterations */
    
    iNI = 0;

    /* Compute the execution configuration */

    dim3 dimBlock(BLOCK_SIZE,BLOCK_SIZE,BLOCK_SIZE);
  
    int blocksInX = (iNy/dimBlock.x) + (!(iNy%dimBlock.x)?0:1);
    int blocksInY = (iNx/dimBlock.y) + (!(iNx%dimBlock.y)?0:1);
    int blocksInZ = (iNz/dimBlock.z) + (!(iNz%dimBlock.z)?0:1);

    dim3 dimGrid ( blocksInX, blocksInY*blocksInZ);

    blocksInX = ((iNy-1)/dimBlock.x) + (!((iNy-1)%dimBlock.x)?0:1);
    int blocksInY_x = (iNx/dimBlock.y) + (!(iNx%dimBlock.y)?0:1);
    blocksInZ = (iNz/dimBlock.z) + (!(iNz%dimBlock.z)?0:1);

    dim3 dimGrid_x (blocksInX, blocksInY_x*blocksInZ);

    blocksInX = (iNy/dimBlock.x) + (!(iNy%dimBlock.x)?0:1);
    int blocksInY_y = ((iNx-1)/dimBlock.y) + (!((iNx-1)%dimBlock.y)?0:1);
    blocksInZ = (iNz/dimBlock.z) + (!(iNz%dimBlock.z)?0:1);

    dim3 dimGrid_y ( blocksInX, blocksInY_y*blocksInZ);

    blocksInX = (iNy/dimBlock.x) + (!(iNy%dimBlock.x)?0:1);
    int blocksInY_z = (iNx/dimBlock.y) + (!(iNx%dimBlock.y)?0:1);
    blocksInZ = ((iNz-1)/dimBlock.z) + (!((iNz-1)%dimBlock.z)?0:1);
 
    dim3 dimGrid_z (blocksInX, blocksInY_z*blocksInZ);
    
    start_time = clock();

    printf("Start computing ......................................... \n\n");

    while( iNI<iNbIters ) 
    { 

        /* update px */
        krnl_1<<< dimGrid, dimBlock>>>(pfpt1_GPU, pfps1_GPU, pfu1_GPU, 
                    pfgk1_GPU, pfdv1_GPU, pfpt2_GPU, pfps2_GPU, pfu2_GPU, 
                    pfgk2_GPU, pfdv2_GPU, pfqq_GPU, cc, iNx, iNy, iNz, SZF, 
                    blocksInY, 1.0f/(float)blocksInY);
       
        krnl_2<<< dimGrid_y, dimBlock>>>(pfbx1_GPU, pfgk1_GPU, pfbx2_GPU, 
                    pfgk2_GPU, steps, iNx, iNy, iNz, SZF, blocksInY_y, 1.0f/(float)blocksInY_y);

        krnl_3<<< dimGrid_x, dimBlock>>>(pfby1_GPU, pfgk1_GPU, pfby2_GPU, 
                    pfgk2_GPU, steps, iNx, iNy, iNz, SZF, blocksInY_x, 1.0f/(float)blocksInY_x);
      
        krnl_z<<<dimGrid_z, dimBlock>>>(pfbz1_GPU, pfgk1_GPU, pfbz2_GPU, pfgk2_GPU, steps, 
                  iNx, iNy, iNz, SZF, blocksInY_z, 1.0f/(float)blocksInY_z);

        /* projection step */
        krnl_4<<< dimGrid, dimBlock>>>(pfbx1_GPU, pfby1_GPU, pfbz1_GPU, pfgk1_GPU, pfpenalty1_GPU, 
                    pfbx2_GPU, pfby2_GPU, pfbz2_GPU, pfgk2_GPU, pfpenalty2_GPU, iNx, iNy,
                    iNz, SZF, blocksInY, 1.0f/(float)blocksInY);
      
        krnl_5<<< dimGrid_y, dimBlock >>>(pfbx1_GPU, pfgk1_GPU, pfbx2_GPU, 
                    pfgk2_GPU, iNx, iNy, iNz, SZF, blocksInY_y, 1.0f/(float)blocksInY_y);
    
        krnl_6<<< dimGrid_x, dimBlock >>>(pfby1_GPU, pfgk1_GPU, pfby2_GPU, 
                    pfgk2_GPU, iNx, iNy, iNz, SZF, blocksInY_x, 1.0f/(float)blocksInY_x);

        krnl_zp<<<dimGrid_z, dimBlock>>>(pfbz1_GPU, pfgk1_GPU, pfbz2_GPU, pfgk2_GPU, iNx, 
                    iNy, iNz, SZF, blocksInY_z, 1.0f/(float)blocksInY_z);

        /* compute the divergence  */
        krnl_7<<< dimGrid, dimBlock>>>(pfbx1_GPU, pfby1_GPU, pfbz1_GPU, pfdv1_GPU, 
                    pfbx2_GPU, pfby2_GPU, pfbz2_GPU, pfdv2_GPU, iNx, iNy, iNz, SZF, 
                    blocksInY, 1.0f/(float)blocksInY);

        /* update ps  */
        krnl_8<<< dimGrid, dimBlock>>>(pfps1_GPU, pfpt1_GPU, pfu1_GPU, pfdv1_GPU, 
                    pfCs1_GPU, pfps2_GPU, pfpt2_GPU, pfu2_GPU, pfdv2_GPU, pfCs2_GPU, 
                    pfqq_GPU, cc, iNx, iNy, iNz, SZF, blocksInY, 1.0f/(float)blocksInY);
        
        /* update pt  */
        krnl_9<<< dimGrid, dimBlock>>>(pfps1_GPU, pfpt1_GPU, pfu1_GPU, pfdv1_GPU, 
                    pfCt1_GPU, pfps2_GPU, pfpt2_GPU, pfu2_GPU, pfdv2_GPU, 
                    pfCt2_GPU, pfqq_GPU, cc, iNx, iNy, iNz, SZF, blocksInY, 1.0f/(float)blocksInY);
   
        /* update qq  */
        krnl_qq<<< dimGrid, dimBlock>>>(pfps1_GPU, pfpt1_GPU, pfu1_GPU, pfdv1_GPU, 
                    pfps2_GPU, pfpt2_GPU, pfu2_GPU, pfdv2_GPU, 
                    pfqq_GPU, beta, cc, iNx, iNy, iNz, SZF, blocksInY, 1.0f/(float)blocksInY);

        /* update multipliers */
        krnl_10<<< dimGrid, dimBlock>>>(pfpt1_GPU, pfdv1_GPU, pfps1_GPU, pfu1_GPU, 
                pfpt2_GPU, pfdv2_GPU, pfps2_GPU, pfu2_GPU, pfqq_GPU,
                FPS_GPU, cc, iNx, iNy, iNz, SZF, blocksInY, 1.0f/(float)blocksInY);

        cudaMemcpy( FPS, FPS_GPU, sizeof(float)*(unsigned)(szImg), cudaMemcpyDeviceToHost);

        fps = 0;
        for (int ii=0; ii< szImg; ii++){
                fps += FPS[ii];
        }

        pfcvg[iNI] = fps / szImg / 2;
        
        if (pfcvg[iNI] <= fError){
            break;
        }

        iNI ++;
     }   

    cudaMemcpy( pfu1, pfu1_GPU, sizeof(float)*(unsigned)(szImg), cudaMemcpyDeviceToHost);
    cudaMemcpy( pfu2, pfu2_GPU, sizeof(float)*(unsigned)(szImg), cudaMemcpyDeviceToHost);

    mexPrintf("Total iteration number = %i\n",iNI);
    end_time = clock();

    /* Outputs (see above) */
    punum[0] = iNI;
    
    /* Free memory */
    free( (float *) pfbx1 );
    free( (float *) pfby1 );
    free( (float *) pfbz1 );
    free( (float *) pfps1 );
    free( (float *) pfpt1 );
    free( (float *) pfgk1 );
    free( (float *) pfdv1 );

    free( (float *) pfbx2 );
    free( (float *) pfby2 );
    free( (float *) pfbz2 );
    free( (float *) pfps2 );
    free( (float *) pfpt2 );
    free( (float *) pfgk2 );
    free( (float *) pfdv2 );

    free( (float *) pfqq );
    free( (float *) FPS );

    //    Free GPU Memory
    cudaFree(pfbx1_GPU);
    cudaFree(pfby1_GPU);
    cudaFree(pfbz1_GPU);
    cudaFree(pfpenalty1_GPU);
    cudaFree(pfps1_GPU);
    cudaFree(pfpt1_GPU);
    cudaFree(pfgk1_GPU);
    cudaFree(pfdv1_GPU);
    cudaFree(pfu1_GPU);
    cudaFree(pfCs1_GPU);
    cudaFree(pfCt1_GPU);

    cudaFree(pfbx2_GPU);
    cudaFree(pfby2_GPU);
    cudaFree(pfbz2_GPU);
    cudaFree(pfpenalty2_GPU);
    cudaFree(pfps2_GPU);
    cudaFree(pfpt2_GPU);
    cudaFree(pfgk2_GPU);
    cudaFree(pfdv2_GPU);
    cudaFree(pfu2_GPU);
    cudaFree(pfCs2_GPU);
    cudaFree(pfCt2_GPU);

    cudaFree(pfqq_GPU);
    cudaFree(FPS_GPU);

    tt[0] = difftime(end_time,start_time)/1000000;
    mexPrintf("\nComputing Time for max-flow = %.4f sec\n \n",tt[0]);
    
    
}
/****************************************/






/**********************************************/
/************** END MAIN FUNCTION *************/
/**********************************************/
