#define BLOCK_SIZE 8
#define MAX(a,b) ( a > b ? a : b )
#define MIN(a,b) ( a <= b ? a : b )
#define SIGN(x) ( x >= 0.0 ? 1.0 : -1.0 )
#define ABS(x) ( (x) > 0.0 ? x : -(x) )
#define SQR(x) (x)*(x)

static __global__ void krnl_1(float *pfpt1, float *pfps1, float *pfu1, 
        float *pfgk1, float *pfdv1, float *pfpt2, float *pfps2, float *pfu2, 
        float *pfgk2, float *pfdv2, float *pfqq, float cc, int iNx, int iNy,
        int iNz, int SZF, int blocksInY, float invBlocksInY)
{
   int blockIdxz = __float2uint_rd(blockIdx.y * invBlocksInY);
   int blockIdxy = blockIdx.y - __umul24(blockIdxz,blocksInY);
   int idx   = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
   int idy   = __mul24(blockIdxy,blockDim.y)+threadIdx.y;
   int idz   = __mul24(blockIdxz,blockDim.z)+threadIdx.z;

   if( idx<iNy && idy<iNx && idz<iNz)
   {
    int index = idx + __mul24(idy, iNy) + __mul24(idz, SZF);
    
    pfgk1[index] = pfdv1[index] - (pfps1[index] - pfpt1[index] 
                + pfu1[index]/cc + pfqq[index]);

    pfgk2[index] = pfdv2[index] - (pfps2[index] - pfpt2[index] 
                + pfu2[index]/cc - pfqq[index]);
   }

}

static __global__ void krnl_2(float *pfbx1, float *pfgk1, float *pfbx2, 
            float *pfgk2, float steps, int iNx, int iNy,
            int iNz, int SZF, int blocksInY, float invBlocksInY)
{
    int blockIdxz = __float2uint_rd(blockIdx.y * invBlocksInY);
    int blockIdxy = blockIdx.y - __umul24(blockIdxz,blocksInY);
    int idx   = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
    int idy   = __mul24(blockIdxy,blockDim.y)+threadIdx.y;
    int idz   = __mul24(blockIdxz,blockDim.z)+threadIdx.z;

    if( idx<iNy && idy<(iNx-1) && idz<iNz)
    {
      int index = idx + __mul24(idy+1, iNy) + __mul24(idz, SZF);

      pfbx1[index] = steps*(pfgk1[index] - pfgk1[index-iNy]) + pfbx1[index];
      pfbx2[index] = steps*(pfgk2[index] - pfgk2[index-iNy]) + pfbx2[index];  
    }
}

static __global__ void krnl_3(float *pfby1, float *pfgk1, float *pfby2, float *pfgk2,
            float steps, int iNx, int iNy, int iNz, int SZF, int blocksInY, float invBlocksInY)
{
    int blockIdxz = __float2uint_rd(blockIdx.y * invBlocksInY);
    int blockIdxy = blockIdx.y - __umul24(blockIdxz,blocksInY);
    int idx   = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
    int idy   = __mul24(blockIdxy,blockDim.y)+threadIdx.y;
    int idz   = __mul24(blockIdxz,blockDim.z)+threadIdx.z;
    
    if( idx<(iNy-1) && idy<iNx && idz<iNz)
    {
      int index =idx + __mul24(idy, iNy) + __mul24(idz, SZF) + 1;
      
      pfby1[index] = steps*(pfgk1[index] - pfgk1[index-1]) + pfby1[index];
      pfby2[index] = steps*(pfgk2[index] - pfgk2[index-1]) + pfby2[index];
    }
}

static __global__ void krnl_z(float *pfbz1, float *pfgk1, float *pfbz2, float *pfgk2, 
        float steps, int iNx, int iNy, int iNz, int SZF, int blocksInY, float invBlocksInY)
{
    int blockIdxz = __float2uint_rd(blockIdx.y * invBlocksInY);
    int blockIdxy = blockIdx.y - __umul24(blockIdxz,blocksInY);
    int idx   = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
    int idy   = __mul24(blockIdxy,blockDim.y)+threadIdx.y;
    int idz   = __mul24(blockIdxz,blockDim.z)+threadIdx.z;
    
    if( idx<iNy && idy<iNx && idz<(iNz-1))
    {
      int index = idx + __mul24(idy, iNy) + __mul24(idz+1, SZF);
    
      pfbz1[index] = __fadd_rz(__fmul_rz(steps, __fadd_rz(pfgk1[index], - pfgk1[index-SZF])), pfbz1[index]);
      pfbz2[index] = __fadd_rz(__fmul_rz(steps, __fadd_rz(pfgk2[index], - pfgk2[index-SZF])), pfbz2[index]);
    }
}

static __global__ void krnl_4(float *pfbx1, float *pfby1, float *pfbz1, float *pfgk1, 
        float *pfpenalty1, float *pfbx2, float *pfby2, float *pfbz2, float *pfgk2, 
        float *pfpenalty2, int iNx, int iNy, int iNz, int SZF, int blocksInY, 
        float invBlocksInY)
{
    int blockIdxz = __float2uint_rd(blockIdx.y * invBlocksInY);
    int blockIdxy = blockIdx.y - __umul24(blockIdxz,blocksInY);
    int idx   = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
    int idy   = __mul24(blockIdxy,blockDim.y)+threadIdx.y;
    int idz   = __mul24(blockIdxz,blockDim.z)+threadIdx.z;
    float fpt;
   
    if( idx<iNy && idy<iNx && idz<iNz)
    {
      int index = idx + __mul24(idy, iNy) + __mul24(idz, SZF);

      fpt = sqrt((SQR(pfbx1[index]) + SQR(pfbx1[index+iNy]) 
            + SQR(pfby1[index]) + SQR(pfby1[index+1]) +
            + SQR(pfbz1[index])+ SQR(pfbz1[index+SZF]))*0.5);
                
      if (fpt > pfpenalty1[index])
          pfgk1[index] = pfpenalty1[index]/fpt;
      else
          pfgk1[index] = 1;

      fpt = sqrt((SQR(pfbx2[index]) + SQR(pfbx2[index+iNy]) 
            + SQR(pfby2[index]) + SQR(pfby2[index+1])
            + SQR(pfbz2[index])+ SQR(pfbz2[index+SZF]))*0.5);
                
      if (fpt > pfpenalty2[index])
          pfgk2[index] = pfpenalty2[index]/fpt;
      else
          pfgk2[index] = 1;

    }
}

static __global__ void krnl_5(float *pfbx1, float *pfgk1, float *pfbx2, 
            float *pfgk2, int iNx, int iNy, int iNz, int SZF, int blocksInY, float invBlocksInY)
{
    int blockIdxz = __float2uint_rd(blockIdx.y * invBlocksInY);
    int blockIdxy = blockIdx.y - __umul24(blockIdxz,blocksInY);
    int idx   = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
    int idy   = __mul24(blockIdxy,blockDim.y)+threadIdx.y;
    int idz   = __mul24(blockIdxz,blockDim.z)+threadIdx.z;

    if( idx<iNy && idy<(iNx-1) && idz<iNz)
    {
      int index = idx + __mul24(idy+1, iNy) + __mul24(idz, SZF);  
       
      pfbx1[index] = (pfgk1[index] + pfgk1[index-iNy])*0.5*pfbx1[index];
      pfbx2[index] = (pfgk2[index] + pfgk2[index-iNy])*0.5*pfbx2[index];
    }
}


static __global__ void krnl_6(float *pfby1, float *pfgk1, float *pfby2, 
            float *pfgk2, int iNx, int iNy, int iNz, int SZF, int blocksInY, float invBlocksInY)
{
    int blockIdxz = __float2uint_rd(blockIdx.y * invBlocksInY);
    int blockIdxy = blockIdx.y - __umul24(blockIdxz,blocksInY);
    int idx   = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
    int idy   = __mul24(blockIdxy,blockDim.y)+threadIdx.y;
    int idz   = __mul24(blockIdxz,blockDim.z)+threadIdx.z;
    
    if( idx<(iNy-1) && idy<iNx && idz<iNz)
    {
      int index = idx + __mul24(idy, iNy) + __mul24(idz, SZF)+1;
      
      pfby1[index] = 0.5*(pfgk1[index] + pfgk1[index-1])*pfby1[index];
      pfby2[index] = 0.5*(pfgk2[index-1] + pfgk2[index])*pfby2[index];
    }
}

static __global__ void krnl_zp(float *pfbz1, float *pfgk1, float *pfbz2, float *pfgk2, int iNx, 
        int iNy, int iNz, int SZF, int blocksInY, float invBlocksInY)
{
    int blockIdxz = __float2uint_rd(blockIdx.y * invBlocksInY);
    int blockIdxy = blockIdx.y - __umul24(blockIdxz,blocksInY);
    int idx   = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
    int idy   = __mul24(blockIdxy,blockDim.y)+threadIdx.y;
    int idz   = __mul24(blockIdxz,blockDim.z)+threadIdx.z;
    
    if( idx<iNy && idy<iNx && idz<(iNz-1))
    {
      int index = idx + __mul24(idy, iNy) + __mul24(idz+1, SZF);

      pfbz1[index] = __fmul_rz(__fmul_rz(__fadd_rz(pfgk1[index], pfgk1[index-SZF]), 0.5), pfbz1[index]);
      pfbz2[index] = __fmul_rz(__fmul_rz(__fadd_rz(pfgk2[index], pfgk2[index-SZF]), 0.5), pfbz2[index]);
    }
}

static __global__ void krnl_7(float *pfbx1, float *pfby1, float *pfbz1, float *pfdv1, 
                  float *pfbx2, float *pfby2, float *pfbz2, float *pfdv2, int iNx, int iNy,
                  int iNz, int SZF, int blocksInY, float invBlocksInY)
{
    int blockIdxz = __float2uint_rd(blockIdx.y * invBlocksInY);
    int blockIdxy = blockIdx.y - __umul24(blockIdxz,blocksInY);
    int idx   = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
    int idy   = __mul24(blockIdxy,blockDim.y)+threadIdx.y;
    int idz   = __mul24(blockIdxz,blockDim.z)+threadIdx.z;

   if( idx<iNy && idy<iNx && idz<iNz)
   {
     int index = idx + __mul24(idy, iNy) + __mul24(idz, SZF);

     pfdv1[index] = pfbx1[index+iNy] - pfbx1[index] 
                + pfby1[index+1] - pfby1[index]
                + pfbz1[index+SZF] - pfbz1[index];

     pfdv2[index] = pfbx2[index+iNy] - pfbx2[index] 
                + pfby2[index+1] - pfby2[index]
                + pfbz2[index+SZF] - pfbz2[index];
    }
}
 

static __global__ void krnl_8(float *pfps1, float *pfpt1, float *pfu1, float *pfdv1, 
                float *pfCs1, float *pfps2, float *pfpt2, float *pfu2, float *pfdv2, 
                float *pfCs2, float *pfqq, float cc, int iNx, int iNy,
                int iNz, int SZF, int blocksInY, float invBlocksInY)
{
    int blockIdxz = __float2uint_rd(blockIdx.y * invBlocksInY);
    int blockIdxy = blockIdx.y - __umul24(blockIdxz,blocksInY);
    int idx   = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
    int idy   = __mul24(blockIdxy,blockDim.y)+threadIdx.y;
    int idz   = __mul24(blockIdxz,blockDim.z)+threadIdx.z;
    float fpt;
    
    if( idx<iNy && idy<iNx && idz<iNz)
    {
      int index = idx + __mul24(idy, iNy) + __mul24(idz, SZF);

      fpt = pfpt1[index] - pfu1[index]/cc + pfdv1[index] - pfqq[index] + 1/cc;
      pfps1[index] = MIN(fpt, pfCs1[index]);

      fpt = pfpt2[index] - pfu2[index]/cc + pfdv2[index] + pfqq[index] + 1/cc;
      pfps2[index] = MIN(fpt , pfCs2[index]);
    }
}


static __global__ void krnl_9(float *pfps1, float *pfpt1, float *pfu1, float *pfdv1, 
                    float *pfCt1, float *pfps2, float *pfpt2, float *pfu2, float *pfdv2, 
                    float *pfCt2, float *pfqq, float cc, int iNx, int iNy,
                    int iNz, int SZF, int blocksInY, float invBlocksInY)
{
    int blockIdxz = __float2uint_rd(blockIdx.y * invBlocksInY);
    int blockIdxy = blockIdx.y - __umul24(blockIdxz,blocksInY);
    int idx   = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
    int idy   = __mul24(blockIdxy,blockDim.y)+threadIdx.y;
    int idz   = __mul24(blockIdxz,blockDim.z)+threadIdx.z;
    float fpt;
    
    if( idx<iNy && idy<iNx && idz<iNz)
    {
      int index = idx + __mul24(idy, iNy) + __mul24(idz, SZF);

      fpt = pfps1[index] + pfu1[index]/cc - pfdv1[index] + pfqq[index];
      pfpt1[index] = MIN(fpt, pfCt1[index]);
        
      fpt = pfps2[index] + pfu2[index]/cc - pfdv2[index] - pfqq[index];
      pfpt2[index] = MIN(fpt , pfCt2[index]);
    }
}

static __global__ void krnl_qq(float *pfps1, float *pfpt1, float *pfu1, float *pfdv1, 
                    float *pfps2, float *pfpt2, float *pfu2, float *pfdv2, 
                    float *pfqq, float beta, float cc, int iNx, int iNy,
                    int iNz, int SZF, int blocksInY, float invBlocksInY)
{
    int blockIdxz = __float2uint_rd(blockIdx.y * invBlocksInY);
    int blockIdxy = blockIdx.y - __umul24(blockIdxz,blocksInY);
    int idx   = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
    int idy   = __mul24(blockIdxy,blockDim.y)+threadIdx.y;
    int idz   = __mul24(blockIdxz,blockDim.z)+threadIdx.z;
    float fpt;
    
    if( idx<iNy && idy<iNx && idz<iNz)
    {
      int index = idx + __mul24(idy, iNy) + __mul24(idz, SZF);

      fpt = ((pfdv1[index] - pfps1[index] - pfu1[index]/cc + pfpt1[index]) +
           (pfps2[index] + pfu2[index]/cc - pfdv2[index] - pfpt2[index]))/2;
      pfqq[index] = MAX(MIN(fpt , beta), -beta);
    }
}

static __global__ void krnl_10(float *pfpt1, float *pfdv1, float *pfps1, float *pfu1,
                    float *pfpt2, float *pfdv2, float *pfps2, float *pfu2, float *pfqq,
                    float *FPS, float cc, int iNx, int iNy, 
                    int iNz, int SZF, int blocksInY, float invBlocksInY)
{
    int blockIdxz = __float2uint_rd(blockIdx.y * invBlocksInY);
    int blockIdxy = blockIdx.y - __umul24(blockIdxz,blocksInY);
    int idx   = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
    int idy   = __mul24(blockIdxy,blockDim.y)+threadIdx.y;
    int idz   = __mul24(blockIdxz,blockDim.z)+threadIdx.z;
    float fpt;

    if( idx<iNy && idy<iNx && idz<iNz)
    {
      int index = idx + __mul24(idy, iNy) + __mul24(idz, SZF);

      fpt = cc*(pfpt1[index] + pfdv1[index] - pfps1[index] - pfqq[index]);
      FPS[index] = ABS(fpt);

      pfu1[index] -= fpt;

      fpt = cc*(pfpt2[index] + pfdv2[index] - pfps2[index] + pfqq[index]);
      FPS[index] += ABS(fpt);

      pfu2[index] -= fpt;
    }
}