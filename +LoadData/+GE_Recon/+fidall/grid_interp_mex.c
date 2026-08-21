#include "mex.h"
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<stdio.h>

/*
assume that it is called with a wrapper to ensure arguments are correct
input:

1. data complex doubles: nPoints x nreps
2. trajectory doubles:   3 x nPoints
3. [gridFOVx gridFOVy] doubles
4. [matrix_x matrix_y] ints

output: 
   complex doubles: matrix_x x matrix_y x nreps 
*/


struct lookup_kernel{
  float * C[2];
  int S_shift[2];
  int W[2];
};

struct kspace_point_2D {
  double kx, ky, weight;
};

struct trajectory2D {
  struct kspace_point_2D * sample;
  mwSize length;
};

void load_kernel_W4_alpha150(struct lookup_kernel * lk, int slot);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  /* declare all variables (welcome to strict C) */
  double * data_real, * data_imag;
  struct trajectory2D trajectory;
  double gridFOVx, gridFOVy;
  mwSize Gx_i, Gy_i;
  double * matrix_real, * matrix_imag;
  mwSize matrix_length;
  mwSize index;
  mwSize nsamples;		/* number of data samples (or trajectory points) */
  mwSize nreps, lrep;   /* number of repetitions (eg coils) */
  struct lookup_kernel lk;
  int slot;
  float * Cx, * Cy; /*, * sinLookup; */
  int Wy, S_shift_x, S_shift_y;
  double Sx_d, Sy_d;
  double SMatrixFOVx, SMatrix_2x, SMatrixFOVy, SMatrix_2y;
  int iPoint;

  /* Get data values */
  nsamples = (mwSize)mxGetM(prhs[0]);
  nreps = (mwSize)mxGetN(prhs[0]);
  data_real = (double*)mxGetPr(prhs[0]);
  data_imag = (double*)mxGetPi(prhs[0]);
  /* mexPrintf("nreps=%d nsamples=%d\n",nreps,nsamples); */
  
  /* get trajectory values */
  trajectory.sample = (struct kspace_point_2D *)mxGetPr(prhs[1]);
  trajectory.length = mxGetN(prhs[1]);
  if (trajectory.length != nsamples) {
    mexPrintf("Error: trajectory.length(=%d) != nsamples(=%d)\n",
            trajectory.length,nsamples);
    return;
  }

  /* get FOV values */
  gridFOVx = mxGetPr(prhs[2])[0];
  gridFOVy = mxGetPr(prhs[2])[1];

  /* get grid matrix size */  
  Gx_i = mxGetPr(prhs[3])[0];
  Gy_i = mxGetPr(prhs[3])[1];

  /* Make space for grid matrix */
  /* plhs[0] = mxCreateDoubleMatrix(Gx_i, Gy_i, mxCOMPLEX); */
  mwSize dims[3];
  dims[0] = Gx_i;
  dims[1] = Gy_i;
  dims[2] = nreps;
    
  plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxCOMPLEX);
  matrix_real = mxGetPr(plhs[0]);
  matrix_imag = mxGetPi(plhs[0]);
  matrix_length = Gx_i * Gy_i * nreps;
  for( index = 0; index < matrix_length; index++ ) {
    matrix_real[index] = 0.0;
    matrix_imag[index] = 0.0;
  }

  /* Initialize our kernel */
  lk.C[0] = (float *)mxMalloc(4096); /*NJS*/
  slot = 0;
  load_kernel_W4_alpha150(&lk, slot);
  Cx = lk.C[slot];
  Cy = lk.C[slot];
  Wy = lk.W[slot];
  S_shift_x = lk.S_shift[slot];
  S_shift_y = lk.S_shift[slot];

  /* igrid */
  Sx_d = 1<<S_shift_x;
  Sy_d = 1<<S_shift_y;
  
  /* X direction */
  SMatrixFOVx  = Sx_d * gridFOVx;
  SMatrix_2x = Sx_d * ((double) Gx_i) * 0.5; 
  
  /* Y direction */
  SMatrixFOVy  = Sy_d * gridFOVy;
  SMatrix_2y = Sy_d * ((double) Gy_i) * 0.5;
  
  for ( iPoint = 0; iPoint < trajectory.length; iPoint++ ) {

    /* Put trajectory point on SMatrix */
    double SMatrix_x_d = SMatrix_2x + SMatrixFOVx * trajectory.sample[iPoint].kx;
    double SMatrix_y_d = SMatrix_2y + SMatrixFOVy * trajectory.sample[iPoint].ky;

    int SMatrix_x_i = (int)SMatrix_x_d;
    int SMatrix_y_i = (int)SMatrix_y_d;
    
    /* find gridMatrix_x/y */
    int gridMatrix_x_i = ((SMatrix_x_i >> S_shift_x) - 1);
    int gridMatrix_y_i = ((SMatrix_y_i >> S_shift_y) - 1);

    /* Make sure this point and kernel are on the gridMatrix */
    double weighted_data_real, weighted_data_imag;
    
    /* Weights Calculation */
    double kerneltheta_x_d = SMatrix_x_d - floor(SMatrix_x_d);
    double kerneltheta_y_d = SMatrix_y_d - floor(SMatrix_y_d);
    int kerneloffsetx = (SMatrix_x_i & 0x3F) << 3; /* << 1 to make it skip 2 4-float vectors */
    int kerneloffsety = (SMatrix_y_i & 0x3F) << 3; /* << 1 to make it skip 2 4-float vectors */
    
    double xWeight[4];
    double yWeight[4];
    
    int ix, iy;
    int xindex[4];
    int yindex[4];
    
    for( ix = 0; ix < 4; ix++ )
        xindex[ix] = (gridMatrix_x_i + ix + Gx_i) % Gx_i;
    for( iy = 0; iy < 4; iy++ )
        yindex[iy] = Gx_i * ((gridMatrix_y_i + iy + Gy_i) % Gy_i);
    
    
    xWeight[0] = Cx[ kerneloffsetx + 4] * kerneltheta_x_d + Cx[ kerneloffsetx ];
    xWeight[1] = Cx[ kerneloffsetx + 5] * kerneltheta_x_d + Cx[ kerneloffsetx + 1 ];
    xWeight[2] = Cx[ kerneloffsetx + 6] * kerneltheta_x_d + Cx[ kerneloffsetx + 2 ];
    xWeight[3] = Cx[ kerneloffsetx + 7] * kerneltheta_x_d + Cx[ kerneloffsetx + 3 ];
    
    yWeight[0] = Cy[ kerneloffsety + 4] * kerneltheta_y_d + Cy[ kerneloffsety ];
    yWeight[1] = Cy[ kerneloffsety + 5] * kerneltheta_y_d + Cy[ kerneloffsety + 1 ];
    yWeight[2] = Cy[ kerneloffsety + 6] * kerneltheta_y_d + Cy[ kerneloffsety + 2 ];
    yWeight[3] = Cy[ kerneloffsety + 7] * kerneltheta_y_d + Cy[ kerneloffsety + 3 ];
    
    
    for (lrep=0; lrep<nreps; lrep++) {
      weighted_data_real = data_real[iPoint+lrep*nsamples] * trajectory.sample[iPoint].weight;
      weighted_data_imag = data_imag[iPoint+lrep*nsamples] * trajectory.sample[iPoint].weight;
        
      /* Perform gridding*/
      for ( iy = 0; iy < Wy; iy++ ) {
        double wdr_y = weighted_data_real * yWeight[iy];
        double wdi_y = weighted_data_imag * yWeight[iy];
        
        double * local_grid_real = matrix_real + yindex[iy] + lrep*Gx_i*Gy_i;
        double * local_grid_imag = matrix_imag + yindex[iy] + lrep*Gx_i*Gy_i;
        
        local_grid_real[xindex[0]] += wdr_y * xWeight[0];
        local_grid_imag[xindex[0]] += wdi_y * xWeight[0];
        
        local_grid_real[xindex[1]] += wdr_y * xWeight[1];
        local_grid_imag[xindex[1]] += wdi_y * xWeight[1];
        
        local_grid_real[xindex[2]] += wdr_y * xWeight[2];
        local_grid_imag[xindex[2]] += wdi_y * xWeight[2];
        
        local_grid_real[xindex[3]] += wdr_y * xWeight[3];
        local_grid_imag[xindex[3]] += wdi_y * xWeight[3];
        
      } /* for iy */
    }   /* for lrep */
  }     /* Point Loop */
   

  mxFree( lk.C[0] );
  return;
};

void load_kernel_W4_alpha150(struct lookup_kernel * lk, int slot) {
unsigned int cc[] = {
0x3EBFAB45, 0x3F800000, 0x3EBFAB45, 0x00000000, 0xBC4862F2, 0xB96BF0DA, 0x3C4AE827, 0x3B5527F8, 
0x3EB9682E, 0x3F7FF141, 0x3EC60286, 0x3B5527F8, 0xBC45BD9A, 0xBA30E2D1, 0x3C4D4BBB, 0x3A42852C, 
0x3EB33A41, 0x3F7FC508, 0x3ECC6CE4, 0x3B82E4A1, 0xBC42F9A5, 0xBA9349FE, 0x3C4F8C38, 0x3A57AB7D, 
0x3EAD2274, 0x3F7F7B63, 0x3ED2E946, 0x3B9DDA11, 0xBC40189B, 0xBACDF619, 0x3C51A82E, 0x3A6E09BB, 
0x3EA721AF, 0x3F7F1468, 0x3ED97687, 0x3BBB9B48, 0xBC3D1C0D, 0xBB043205, 0x3C539E38, 0x3A82D2B6, 
0x3EA138CE, 0x3F7E9036, 0x3EE01379, 0x3BDC4FF6, 0xBC3A058E, 0xBB21411E, 0x3C556CF9, 0x3A8F41D2, 
0x3E9B68A2, 0x3F7DEEF5, 0x3EE6BEE1, 0x3C001035, 0xBC36D6B3, 0xBB3E1FA1, 0x3C57131F, 0x3A9C547B, 
0x3E95B1EC, 0x3F7D30D5, 0x3EED777A, 0x3C139AC5, 0xBC339116, 0xBB5AC4EA, 0x3C588F63, 0x3AAA0CBA, 
0x3E901564, 0x3F7C5611, 0x3EF43BF5, 0x3C28DC5C, 0xBC303650, 0xBB77286F, 0x3C59E08B, 0x3AB86C53, 
0x3E8A93B1, 0x3F7B5EE8, 0x3EFB0AF9, 0x3C3FE9E6, 0xBC2CC7FE, 0xBB89A0E2, 0x3C5B0568, 0x3AC774C1, 
0x3E852D71, 0x3F7A4BA6, 0x3F00F192, 0x3C58D87E, 0xBC2947BA, 0xBB97844F, 0x3C5BFCD9, 0x3AD72734, 
0x3E7FC667, 0x3F791C9E, 0x3F046186, 0x3C73BD65, 0xBC25B71E, 0xBBA53A6A, 0x3C5CC5CC, 0x3AE7848F, 
0x3E756AF5, 0x3F77D229, 0x3F07D49D, 0x3C8856FB, 0xBC2217C5, 0xBBB2BF33, 0x3C5D5F3B, 0x3AF88D61, 
0x3E6B4979, 0x3F766CAA, 0x3F0B4A1A, 0x3C97DFD1, 0xBC1E6B45, 0xBBC00EC0, 0x3C5DC832, 0x3B0520F4, 
0x3E6162C4, 0x3F74EC8D, 0x3F0EC13B, 0x3CA883F0, 0xBC1AB331, 0xBBCD253E, 0x3C5DFFCD, 0x3B0E5105, 
0x3E57B791, 0x3F735242, 0x3F12393A, 0x3CBA4E11, 0xBC16F11B, 0xBBD9FEF2, 0x3C5E0537, 0x3B17D6AA, 
0x3E4E487F, 0x3F719E45, 0x3F15B14F, 0x3CCD48E6, 0xBC13268E, 0xBBE6983E, 0x3C5DD7AC, 0x3B21B17E, 
0x3E451617, 0x3F6FD114, 0x3F1928AD, 0x3CE17F16, 0xBC0F5511, 0xBBF2ED9C, 0x3C5D767C, 0x3B2BE0E9, 
0x3E3C20C6, 0x3F6DEB39, 0x3F1C9E87, 0x3CF6FB33, 0xBC0B7E25, 0xBBFEFBA8, 0x3C5CE106, 0x3B366428, 
0x3E3368E3, 0x3F6BED42, 0x3F20120B, 0x3D06E3DC, 0xBC07A346, 0xBC055F8D, 0x3C5C16BF, 0x3B413A44, 
0x3E2AEEAF, 0x3F69D7C3, 0x3F238266, 0x3D12F780, 0xBC03C5E8, 0xBC0B1A65, 0x3C5B172D, 0x3B4C6217, 
0x3E22B250, 0x3F67AB5A, 0x3F26EEC3, 0x3D1FBDA1, 0xBBFFCEEF, 0xBC10ACD9, 0x3C59E1EA, 0x3B57DA48, 
0x3E1AB3D9, 0x3F6568A6, 0x3F2A564B, 0x3D2D3B46, 0xBBF812B2, 0xBC161578, 0x3C5876A4, 0x3B63A14B, 
0x3E12F343, 0x3F631051, 0x3F2DB825, 0x3D3B755B, 0xBBF059D2, 0xBC1B52E3, 0x3C56D51D, 0x3B6FB562, 
0x3E0B7075, 0x3F60A305, 0x3F31137A, 0x3D4A70B1, 0xBBE8A6F6, 0xBC2063CC, 0x3C54FD2C, 0x3B7C149B, 
0x3E042B3D, 0x3F5E2176, 0x3F34676E, 0x3D5A31FA, 0xBBE0FCAF, 0xBC2546F8, 0x3C52EEBC, 0x3B845E68, 
0x3DFA46AF, 0x3F5B8C5A, 0x3F37B329, 0x3D6ABDC7, 0xBBD95D7F, 0xBC29FB41, 0x3C50A9CE, 0x3B8AD5D3, 
0x3DECB0D7, 0x3F58E46D, 0x3F3AF5D1, 0x3D7C1882, 0xBBD1CBCF, 0xBC2E7F93, 0x3C4E2E78, 0x3B916F48, 
0x3DDF941A, 0x3F562A6F, 0x3F3E2E8A, 0x3D872335, 0xBBCA49F4, 0xBC32D2F0, 0x3C4B7CE4, 0x3B982967, 
0x3DD2EF7B, 0x3F535F23, 0x3F415C7E, 0x3D90A5CC, 0xBBC2DA2E, 0xBC36F46C, 0x3C489554, 0x3B9F02B4, 
0x3DC6C1D8, 0x3F508351, 0x3F447ED3, 0x3D9A95F7, 0xBBBB7EA4, 0xBC3AE332, 0x3C45781D, 0x3BA5F99B, 
0x3DBB09EE, 0x3F4D97C4, 0x3F4794B4, 0x3DA4F591, 0xBBB43966, 0xBC3E9E81, 0x3C4225AC, 0x3BAD0C6E, 
0x3DAFC658, 0x3F4A9D4A, 0x3F4A9D4A, 0x3DAFC658, 0xBBAD0C6E, 0xBC4225AC, 0x3C3E9E81, 0x3BB43966, 
0x3DA4F591, 0x3F4794B4, 0x3F4D97C4, 0x3DBB09EE, 0xBBA5F99B, 0xBC45781D, 0x3C3AE332, 0x3BBB7EA4, 
0x3D9A95F7, 0x3F447ED3, 0x3F508351, 0x3DC6C1D8, 0xBB9F02B4, 0xBC489554, 0x3C36F46C, 0x3BC2DA2E, 
0x3D90A5CC, 0x3F415C7E, 0x3F535F23, 0x3DD2EF7B, 0xBB982967, 0xBC4B7CE4, 0x3C32D2F0, 0x3BCA49F4, 
0x3D872335, 0x3F3E2E8A, 0x3F562A6F, 0x3DDF941A, 0xBB916F48, 0xBC4E2E78, 0x3C2E7F93, 0x3BD1CBCF, 
0x3D7C1882, 0x3F3AF5D1, 0x3F58E46D, 0x3DECB0D7, 0xBB8AD5D3, 0xBC50A9CE, 0x3C29FB41, 0x3BD95D7F, 
0x3D6ABDC7, 0x3F37B329, 0x3F5B8C5A, 0x3DFA46AF, 0xBB845E68, 0xBC52EEBC, 0x3C2546F8, 0x3BE0FCAF, 
0x3D5A31FA, 0x3F34676E, 0x3F5E2176, 0x3E042B3D, 0xBB7C149B, 0xBC54FD2C, 0x3C2063CC, 0x3BE8A6F6, 
0x3D4A70B1, 0x3F31137A, 0x3F60A305, 0x3E0B7075, 0xBB6FB562, 0xBC56D51D, 0x3C1B52E3, 0x3BF059D2, 
0x3D3B755B, 0x3F2DB825, 0x3F631051, 0x3E12F343, 0xBB63A14B, 0xBC5876A4, 0x3C161578, 0x3BF812B2, 
0x3D2D3B46, 0x3F2A564B, 0x3F6568A6, 0x3E1AB3D9, 0xBB57DA48, 0xBC59E1EA, 0x3C10ACD9, 0x3BFFCEEF, 
0x3D1FBDA1, 0x3F26EEC3, 0x3F67AB5A, 0x3E22B250, 0xBB4C6217, 0xBC5B172D, 0x3C0B1A65, 0x3C03C5E8, 
0x3D12F780, 0x3F238266, 0x3F69D7C3, 0x3E2AEEAF, 0xBB413A44, 0xBC5C16BF, 0x3C055F8D, 0x3C07A346, 
0x3D06E3DC, 0x3F20120B, 0x3F6BED42, 0x3E3368E3, 0xBB366428, 0xBC5CE106, 0x3BFEFBA8, 0x3C0B7E25, 
0x3CF6FB33, 0x3F1C9E87, 0x3F6DEB39, 0x3E3C20C6, 0xBB2BE0E9, 0xBC5D767C, 0x3BF2ED9C, 0x3C0F5511, 
0x3CE17F16, 0x3F1928AD, 0x3F6FD114, 0x3E451617, 0xBB21B17E, 0xBC5DD7AC, 0x3BE6983E, 0x3C13268E, 
0x3CCD48E6, 0x3F15B14F, 0x3F719E45, 0x3E4E487F, 0xBB17D6AA, 0xBC5E0537, 0x3BD9FEF2, 0x3C16F11B, 
0x3CBA4E11, 0x3F12393A, 0x3F735242, 0x3E57B791, 0xBB0E5105, 0xBC5DFFCD, 0x3BCD253E, 0x3C1AB331, 
0x3CA883F0, 0x3F0EC13B, 0x3F74EC8D, 0x3E6162C4, 0xBB0520F4, 0xBC5DC832, 0x3BC00EC0, 0x3C1E6B45, 
0x3C97DFD1, 0x3F0B4A1A, 0x3F766CAA, 0x3E6B4979, 0xBAF88D61, 0xBC5D5F3B, 0x3BB2BF33, 0x3C2217C5, 
0x3C8856FB, 0x3F07D49D, 0x3F77D229, 0x3E756AF5, 0xBAE7848F, 0xBC5CC5CC, 0x3BA53A6A, 0x3C25B71E, 
0x3C73BD65, 0x3F046186, 0x3F791C9E, 0x3E7FC667, 0xBAD72734, 0xBC5BFCD9, 0x3B97844F, 0x3C2947BA, 
0x3C58D87E, 0x3F00F192, 0x3F7A4BA6, 0x3E852D71, 0xBAC774C1, 0xBC5B0568, 0x3B89A0E2, 0x3C2CC7FE, 
0x3C3FE9E6, 0x3EFB0AF9, 0x3F7B5EE8, 0x3E8A93B1, 0xBAB86C53, 0xBC59E08B, 0x3B77286F, 0x3C303650, 
0x3C28DC5C, 0x3EF43BF5, 0x3F7C5611, 0x3E901564, 0xBAAA0CBA, 0xBC588F63, 0x3B5AC4EA, 0x3C339116, 
0x3C139AC5, 0x3EED777A, 0x3F7D30D5, 0x3E95B1EC, 0xBA9C547B, 0xBC57131F, 0x3B3E1FA1, 0x3C36D6B3, 
0x3C001035, 0x3EE6BEE1, 0x3F7DEEF5, 0x3E9B68A2, 0xBA8F41D2, 0xBC556CF9, 0x3B21411E, 0x3C3A058E, 
0x3BDC4FF6, 0x3EE01379, 0x3F7E9036, 0x3EA138CE, 0xBA82D2B6, 0xBC539E38, 0x3B043205, 0x3C3D1C0D, 
0x3BBB9B48, 0x3ED97687, 0x3F7F1468, 0x3EA721AF, 0xBA6E09BB, 0xBC51A82E, 0x3ACDF619, 0x3C40189B, 
0x3B9DDA11, 0x3ED2E946, 0x3F7F7B63, 0x3EAD2274, 0xBA57AB7D, 0xBC4F8C38, 0x3A9349FE, 0x3C42F9A5, 
0x3B82E4A1, 0x3ECC6CE4, 0x3F7FC508, 0x3EB33A41, 0xBA42852C, 0xBC4D4BBB, 0x3A30E2D1, 0x3C45BD9A, 
0x3B5527F8, 0x3EC60286, 0x3F7FF141, 0x3EB9682E, 0xBB5527F8, 0xBC4AE827, 0x396BF0DA, 0x3C4862F2
};
memcpy(lk->C[slot], cc, 2048);
lk->S_shift[slot] = 6;
lk->W[slot] = 4;
};


