#define MAX(a,b) ( a > b ? a : b )
#define MIN(a,b) ( a <= b ? a : b )

static __global__ void krnl_1(float *dvx, float *dvy, float *dvz,
        float *ux, float *uy, float *uz,
        float *gkx, float *gky, float *gkz,
        float *Gx, float *Gy, float *Gz,
        float *Gt, float *Gf, float *q,
        float *Ux, float *Uy, float *Uz,
        float cc, int iNx, int iNy, int iNz, float *d_qq1, float *d_qq2, float *d_qq3, int flag, float weight){
    
    int idx = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * (blockIdx.y + gridDim.y * blockIdx.z));
    if (idx >= iNx*iNy*iNz) return;

    if( ( (idx%iNy) != (iNy-1) ) &&
            ( (idx/(iNx*iNy)) < (iNz-1) ) &&
            ( ((idx/iNy)%iNx) != (iNx-1))
            ){
        float tx = dvx[idx] + flag*d_qq1[idx] - ux[idx]/cc;
        float ty = dvy[idx] + flag*d_qq2[idx] - uy[idx]/cc;
        float tz = dvz[idx] + flag*d_qq3[idx] - uz[idx]/cc;
        
        q[idx] = (Gt[idx] - cc*(Gx[idx]*tx + Gy[idx]*ty + Gz[idx]*tz))/(1 + cc*Gf[idx]);
        q[idx] = MAX(MIN(q[idx],weight),-weight);     /*this was added by Fumin Guo, 2015/08/28*/  
        
        gkx[idx] = tx + q[idx]*Gx[idx] - Ux[idx]/cc*2;
        gky[idx] = ty + q[idx]*Gy[idx] - Uy[idx]/cc*2;
        gkz[idx] = tz + q[idx]*Gz[idx] - Uz[idx]/cc*2;
        
    }
}

static __global__ void krnl_2(float *bx1, float *by1, float *bz1,
        float *gkx, float *gky, float *gkz,
        float steps, int iNx, int iNy, int iNz){
    
    int idx = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * (blockIdx.y + gridDim.y * blockIdx.z));
    if (idx >= iNx*iNy*iNz) return;
    
    if( ( (idx%iNy) != (iNy-1) ) &&
            ( (idx/(iNx*iNy)) < (iNz-1) ) &&
            ( ((idx/iNy)%iNx) != (iNx-1))
            ){
        
        bx1[idx+iNy] = steps*(gkx[idx+iNy] - gkx[idx]) + bx1[idx+iNy];
        by1[idx+iNy] = steps*(gky[idx+iNy] - gky[idx]) + by1[idx+iNy];
        bz1[idx+iNy] = steps*(gkz[idx+iNy] - gkz[idx]) + bz1[idx+iNy];
    }
}

static __global__ void krnl_3(float *bx2, float *by2, float *bz2,
        float *gkx, float *gky, float *gkz,
        float steps, int iNx, int iNy, int iNz){
    
    int idx = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * (blockIdx.y + gridDim.y * blockIdx.z));
    if (idx >= iNx*iNy*iNz) return;
    
    if( ( (idx%iNy) != (iNy-1) ) &&
            ( (idx/(iNx*iNy)) < (iNz-1) ) &&
            ( ((idx/iNy)%iNx) != (iNx-1))
            ){
        
        bx2[idx+1] = steps*(gkx[idx+1] - gkx[idx]) + bx2[idx+1];
        by2[idx+1] = steps*(gky[idx+1] - gky[idx]) + by2[idx+1];
        bz2[idx+1] = steps*(gkz[idx+1] - gkz[idx]) + bz2[idx+1];
        
    }
}

static __global__ void krnl_z(float *bx3, float *by3, float *bz3,
        float *gkx, float *gky, float *gkz,
        float steps, int iNx, int iNy, int iNz){
    
    int idx = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * (blockIdx.y + gridDim.y * blockIdx.z));
    if (idx >= iNx*iNy*iNz) return;
    
    if( ( (idx%iNy) != (iNy-1) ) &&
            ( (idx/(iNx*iNy)) < (iNz-1) ) &&
            ( ((idx/iNy)%iNx) != (iNx-1))
            ){
        
        bx3[idx+(iNx*iNy)] = steps*(gkx[idx+(iNx*iNy)] - gkx[idx]) + bx3[idx+(iNx*iNy)];
        by3[idx+(iNx*iNy)] = steps*(gky[idx+(iNx*iNy)] - gky[idx]) + by3[idx+(iNx*iNy)];
        bz3[idx+(iNx*iNy)] = steps*(gkz[idx+(iNx*iNy)] - gkz[idx]) + bz3[idx+(iNx*iNy)];
    }
}

static __global__ void krnl_23z(float *bx1, float *by1, float *bz1,
        float *bx2, float *by2, float *bz2,
        float *bx3, float *by3, float *bz3,
        float *gkx, float *gky, float *gkz,
        float steps, int iNx, int iNy, int iNz){
    
    int idx = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * (blockIdx.y + gridDim.y * blockIdx.z));
    if (idx >= iNx*iNy*iNz) return;
    
    if( ( (idx%iNy) != (iNy-1) ) &&
            ( (idx/(iNx*iNy)) < (iNz-1) ) &&
            ( ((idx/iNy)%iNx) != (iNx-1))
            ){
        bx1[idx+iNy] = steps*(gkx[idx+iNy] - gkx[idx]) + bx1[idx+iNy];
        by1[idx+iNy] = steps*(gky[idx+iNy] - gky[idx]) + by1[idx+iNy];
        bz1[idx+iNy] = steps*(gkz[idx+iNy] - gkz[idx]) + bz1[idx+iNy];
        
        bx2[idx+1] = steps*(gkx[idx+1] - gkx[idx]) + bx2[idx+1];
        by2[idx+1] = steps*(gky[idx+1] - gky[idx]) + by2[idx+1];
        bz2[idx+1] = steps*(gkz[idx+1] - gkz[idx]) + bz2[idx+1];
        
        bx3[idx+(iNx*iNy)] = steps*(gkx[idx+(iNx*iNy)] - gkx[idx]) + bx3[idx+(iNx*iNy)];
        by3[idx+(iNx*iNy)] = steps*(gky[idx+(iNx*iNy)] - gky[idx]) + by3[idx+(iNx*iNy)];
        bz3[idx+(iNx*iNy)] = steps*(gkz[idx+(iNx*iNy)] - gkz[idx]) + bz3[idx+(iNx*iNy)];
        
    }
}

static __global__ void krnl_4(float *bx1, float *bx2, float *bx3,
        float *by1, float *by2, float *by3,
        float *bz1, float *bz2, float *bz3,
        float *gkx, float *gky, float *gkz,
        float fPenalty, int iNx, int iNy, int iNz){
    
    int idx = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * (blockIdx.y + gridDim.y * blockIdx.z));
    if (idx >= iNx*iNy*iNz) return;
    
    if( ( (idx%iNy) != (iNy-1) ) &&
            ( (idx/(iNx*iNy)) < (iNz-1) ) &&
            ( ((idx/iNy)%iNx) != (iNx-1))
            ){
        
        
        
        float fpt = sqrtf((powf(bx1[idx],2) + powf(bx1[idx+iNy],2)
        + powf(bx2[idx],2) + powf(bx2[idx+1],2)
        + powf(bx3[idx],2) + powf(bx3[idx+(iNx*iNy)],2))*0.5);
        
        if (fpt > fPenalty)
            gkx[idx] = fPenalty/fpt;
        else
            gkx[idx] = 1;
        
        fpt = sqrtf((powf(by1[idx],2) + powf(by1[idx+iNy],2)
        + powf(by2[idx],2) + powf(by2[idx+1],2)
        + powf(by3[idx],2) + powf(by3[idx+(iNx*iNy)],2))*0.5);
        
        if (fpt > fPenalty)
            gky[idx] = fPenalty/fpt;
        else
            gky[idx] = 1;
        
        fpt = sqrtf((powf(bz1[idx],2) + powf(bz1[idx+iNy],2)
        + powf(bz2[idx],2) + powf(bz2[idx+1],2)
        + powf(bz3[idx],2) + powf(bz3[idx+(iNx*iNy)],2))*0.5);
        
        if (fpt > fPenalty)
            gkz[idx] = fPenalty/fpt;
        else
            gkz[idx] = 1;
        
    }
}

static __global__ void krnl_5(float *bx1, float *by1, float *bz1,
        float *gkx, float *gky, float *gkz,
        int iNx, int iNy, int iNz){
    
    int idx = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * (blockIdx.y + gridDim.y * blockIdx.z));
    if (idx >= iNx*iNy*iNz) return;
    
    if( ( (idx%iNy) != (iNy-1) ) &&
            ( (idx/(iNx*iNy)) < (iNz-1) ) &&
            ( ((idx/iNy)%iNx) != (iNx-1))
            ){
        
        
        bx1[idx+iNy] = (gkx[idx+iNy] + gkx[idx])*0.5*bx1[idx+iNy];
        by1[idx+iNy] = (gky[idx+iNy] + gky[idx])*0.5*by1[idx+iNy];
        bz1[idx+iNy] = (gkz[idx+iNy] + gkz[idx])*0.5*bz1[idx+iNy];
    }
}

static __global__ void krnl_6(float *bx2, float *by2, float *bz2,
        float *gkx, float *gky, float *gkz,
        int iNx, int iNy, int iNz){
    
    int idx = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * (blockIdx.y + gridDim.y * blockIdx.z));
    if (idx >= iNx*iNy*iNz) return;
    
    if( ( (idx%iNy) != (iNy-1) ) &&
            ( (idx/(iNx*iNy)) < (iNz-1) ) &&
            ( ((idx/iNy)%iNx) != (iNx-1))
            ){
        
        bx2[idx+1] = 0.5*(gkx[idx+1] + gkx[idx])*bx2[idx+1];
        by2[idx+1] = 0.5*(gky[idx+1] + gky[idx])*by2[idx+1];
        bz2[idx+1] = 0.5*(gkz[idx+1] + gkz[idx])*bz2[idx+1];
    }
}

static __global__ void krnl_zp(float *bx3, float *by3, float *bz3,
        float *gkx, float *gky, float *gkz,
        int iNx, int iNy, int iNz){
    
    int idx = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * (blockIdx.y + gridDim.y * blockIdx.z));
    if (idx >= iNx*iNy*iNz) return;

    if( ( (idx%iNy) != (iNy-1) ) &&
            ( (idx/(iNx*iNy)) < (iNz-1) ) &&
            ( ((idx/iNy)%iNx) != (iNx-1))
            ){
        
        bx3[idx+(iNx*iNy)] = 0.5*(gkx[idx+(iNx*iNy)] + gkx[idx])*bx3[idx+(iNx*iNy)];
        by3[idx+(iNx*iNy)] = 0.5*(gky[idx+(iNx*iNy)] + gky[idx])*by3[idx+(iNx*iNy)];
        bz3[idx+(iNx*iNy)] = 0.5*(gkz[idx+(iNx*iNy)] + gkz[idx])*bz3[idx+(iNx*iNy)];
    }
}

static __global__ void krnl_56zp(float *bx1, float *by1, float *bz1,
        float *bx2, float *by2, float *bz2,
        float *bx3, float *by3, float *bz3,
        float *gkx, float *gky, float *gkz,
        int iNx, int iNy, int iNz){
    
    int idx = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * (blockIdx.y + gridDim.y * blockIdx.z));
    if (idx >= iNx*iNy*iNz) return;
    
    if( ( (idx%iNy) != (iNy-1) ) &&
            ( (idx/(iNx*iNy)) < (iNz-1) ) &&
            ( ((idx/iNy)%iNx) != (iNx-1))
            ){
        
        bx1[idx+iNy] = (gkx[idx+iNy] + gkx[idx])*0.5*bx1[idx+iNy];
        by1[idx+iNy] = (gky[idx+iNy] + gky[idx])*0.5*by1[idx+iNy];
        bz1[idx+iNy] = (gkz[idx+iNy] + gkz[idx])*0.5*bz1[idx+iNy];
        
        bx2[idx+1] = 0.5*(gkx[idx+1] + gkx[idx])*bx2[idx+1];
        by2[idx+1] = 0.5*(gky[idx+1] + gky[idx])*by2[idx+1];
        bz2[idx+1] = 0.5*(gkz[idx+1] + gkz[idx])*bz2[idx+1];
        
        bx3[idx+(iNx*iNy)] = 0.5*(gkx[idx+(iNx*iNy)] + gkx[idx])*bx3[idx+(iNx*iNy)];
        by3[idx+(iNx*iNy)] = 0.5*(gky[idx+(iNx*iNy)] + gky[idx])*by3[idx+(iNx*iNy)];
        bz3[idx+(iNx*iNy)] = 0.5*(gkz[idx+(iNx*iNy)] + gkz[idx])*bz3[idx+(iNx*iNy)];
        
        
    }
}

static __global__ void krnl_7(float *bx1, float *bx2, float *bx3,
        float *by1, float *by2, float *by3,
        float *bz1, float *bz2, float *bz3,
        float *dvx, float *dvy, float *dvz,
        int iNx, int iNy, int iNz){
    
    int idx = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * (blockIdx.y + gridDim.y * blockIdx.z));
    if (idx >= iNx*iNy*iNz) return;
    
    if( ( (idx%iNy) != (iNy-1) ) &&
            ( (idx/(iNx*iNy)) < (iNz-1) ) &&
            ( ((idx/iNy)%iNx) != (iNx-1))
            ){
        
        
        dvx[idx] = bx1[idx+iNy] - bx1[idx]
                + bx2[idx+1] - bx2[idx]
                + bx3[idx+(iNx*iNy)] - bx3[idx];
        
        dvy[idx] = by1[idx+iNy] - by1[idx]
                + by2[idx+1] - by2[idx]
                + by3[idx+(iNx*iNy)] - by3[idx];
        
        dvz[idx] = bz1[idx+iNy] - bz1[idx]
                + bz2[idx+1] - bz2[idx]
                + bz3[idx+(iNx*iNy)] - bz3[idx];        
    }
}

static __global__ void krnl_8(float *U1x, float *U1y, float *U1z, float *U2x, float *U2y, float *U2z,
        float *q1, float *q2,
        float *G1x, float *G1y, float *G1z, float *G2x, float *G2y, float *G2z,
        float *dv1x, float *dv1y, float *dv1z, float *dv2x, float *dv2y, float *dv2z,
        float *u1x, float *u1y, float *u1z, float *u2x, float *u2y, float *u2z,    
        float cc, float fPenalty, int iNx, int iNy, int iNz, float *d_qq1, float *d_qq2, float *d_qq3){
    
    int idx = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * (blockIdx.y + gridDim.y * blockIdx.z));
    if (idx >= iNx*iNy*iNz) return;
    
    if( ( (idx%iNy) != (iNy-1) ) &&
            ( (idx/(iNx*iNy)) < (iNz-1) ) &&
            ( ((idx/iNy)%iNx) != (iNx-1))
            ){
        

        d_qq1[idx] = U1x[idx] - U2x[idx]
                    + cc*(q2[idx]*G2x[idx] + dv2x[idx] - u2x[idx] - 2*U2x[idx])
                    - cc*(q1[idx]*G1x[idx] + dv1x[idx] - u1x[idx] - 2*U1x[idx]);
        d_qq1[idx] /=  (2*cc);
        d_qq1[idx] = MAX(MIN(d_qq1[idx],fPenalty),-fPenalty);         

        d_qq2[idx] = U1y[idx] - U2y[idx]
                    + cc*(q2[idx]*G2y[idx] + dv2y[idx] - u2y[idx] - 2*U2y[idx])
                    - cc*(q1[idx]*G1y[idx] + dv1y[idx] - u1y[idx] - 2*U1y[idx]);
        d_qq2[idx] /=  (2*cc);
        d_qq2[idx] = MAX(MIN(d_qq2[idx],fPenalty),-fPenalty);  

        d_qq3[idx] = U1z[idx] - U2z[idx]
                    + cc*(q2[idx]*G2z[idx] + dv2z[idx] - u2z[idx] - 2*U2z[idx])
                    - cc*(q1[idx]*G1z[idx] + dv1z[idx] - u1z[idx] - 2*U1z[idx]);
        d_qq3[idx] /=  (2*cc);
        d_qq3[idx] = MAX(MIN(d_qq3[idx],fPenalty),-fPenalty);  
    }
}


static __global__ void krnl_9(float *dv1x, float *dv1y, float *dv1z, float *dv2x, float *dv2y, float *dv2z,
        float *G1x, float *G1y, float *G1z, float *G2x, float *G2y, float *G2z,
        float *q1, float *q2, float *u1x, float *u1y, float *u1z, float *u2x, float *u2y, float *u2z, 
        float *FPS, float  cc, int iNx, int iNy, int iNz, float *d_qq1, float *d_qq2, float *d_qq3){
    
    int idx = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * (blockIdx.y + gridDim.y * blockIdx.z));
    if (idx >= iNx*iNy*iNz) return;
    
    if( ( (idx%iNy) != (iNy-1) ) &&
            ( (idx/(iNx*iNy)) < (iNz-1) ) &&
            ( ((idx/iNy)%iNx) != (iNx-1))
            ){
        
        /* update ux */
        float fp1t = cc*(dv1x[idx] + q1[idx]*G1x[idx] + d_qq1[idx]);
        float fp2t = cc*(dv2x[idx] + q2[idx]*G2x[idx] - d_qq1[idx]);
        FPS[idx] = fabsf(fp1t) + fabsf(fp2t);
        
        u1x[idx] -= fp1t;
        u2x[idx] -= fp2t;

        /* update uy */
        
        fp1t = cc*(dv1y[idx] + q1[idx]*G1y[idx] + d_qq2[idx]);
        fp2t = cc*(dv2y[idx] + q2[idx]*G2y[idx] - d_qq2[idx]);
        FPS[idx] += fabsf(fp1t) + fabsf(fp2t);
        u1y[idx]  -= fp1t;
        u2y[idx]  -= fp2t;

        /* update uz */
        
        fp1t = cc*(dv1z[idx] + q1[idx]*G1z[idx] + d_qq3[idx]);
        fp2t = cc*(dv2z[idx] + q2[idx]*G2z[idx] - d_qq3[idx]);
        FPS[idx] += fabsf(fp1t) + fabsf(fp2t);
        u1z[idx] -= fp1t;
        u2z[idx] -= fp2t;
    }
}

