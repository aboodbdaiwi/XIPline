/* ********************************************
 * Spiral Generation code
 **********************************************
 * Author: Jim Pipe
 * Date: May 2011
 * Rev: Oct 2012
 **********************************************
 * A Subset of Relevant Literature
 *
 * Spiral Invented:
 * High-speed spiral-scan echo planar NMR imaging-I.
 * Ahn, C.B., Kim, J.H. & Cho, Z.H., IEEE Transactions on Medical Imaging, 5(1) 1986.
 *
 * Spiral Improved:
 * Fast Spiral Coronary Artery Imaging.
 * Meyer CH, Hu BS, Nishimura DG, Macovski A, Magnetic Resonance in Medicine, 28(2) 1992.
 *
 * Variable Density Spiral
 * Reduced aliasing artifacts using variable-density k-space sampling trajectories.
 * Tsai CM, Nishimura DG, Magnetic Resonance in Medicine, 43(3), 2000
 *
 * "SLOPPY" SPIRAL
 * Faster Imaging with Randomly Perturbed Undersampled Spirals and L_1 Reconstruction
 * M. Lustig, J.H. Lee, D.L. Donoho, J.M. Pauly, Proc. of the ISMRM '05
 *
 * FLORET
 * A new design and rationale for 3D orthogonally oversampled k-space trajectories
 * Pipe JG, Zwart NR, Aboussouan EA, Robison RK, Devaraj A, Johnson KO, Mag Res Med 66(5) 2011
 *
 * Distributed Spirals
 * Distributed Spirals: A New Class of 3D k-Space Trajectories
 * Turley D, Pipe JG, Magnetic Resonance in Medicine, in press (also proc of ISMRM '12)
 */
 /* RFS: includes */
#include "mex.h"
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define DEBUG
#ifndef MAX
#define UNDEFMAX
#define MAX(a,b) (((a)<(b))?(b):(a))
#endif

#ifndef MIN
#define UNDEFMIN
#define MIN(a,b) (((a)>(b))?(b):(a))
#endif

/* #define GRAST    0.005 * Desired gradient raster time (ms) */
/* #define subrast    5      * number of numerical cycles per gradient raster time */ 
#define GRAST    0.004 /* RFS: Desired gradient raster time (ms) */
#define subrast    4      /* RFS: number of numerical cycles per gradient raster time */ 


/* Data passed in spparams array which is spARRSIZE large */

#define spARRSIZE  20

#define spGAMMA     0
#define spGMAX      1
#define spSLEWMAX   2

#define spGTYPE     3

#define spFOVXY     4
#define spFOVZ      5
#define spRESXY     6
#define spRESZ      7
#define spARMS      8

#define spSTYPE     9
#define spUSTYPE   10
#define spUS0      11
#define spUS1      12
#define spUSR      13

#define spDWELL    14
#define spREADPTS  15

#define spSLOP_PER 16

#define spDEBUG 17

#define spGXSTART 18    /* RFS */

void bnispiralgen(double* spparams, int maxarray, 
        double *gxarray, double *gyarray, double *gzarray,
        int *spgrad_na, int *spgrad_nb, int *spgrad_nc, int *spgrad_nd) {
  /* ***********************************************************
   ************************************************************

  This function takes parameters passed in spparams array and
  returns a single spiral arm calculated numerically

  The corresponding gradient waveforms are in gxarray and gyarray
  spgrad_na reflects the number of gradient points to reach the end of k-space
  spgrad_nb = spgrad_na + the number of gradient points to ramp G to zero
  spgrad_nc = spgrad_nb + the number of gradient points to rewind k to zero
  spgrad_nd = spgrad_nc + the number of gradient points for first moment compensation

  Assignments below indicate units of input parameters
  All units input using kHz, msec, mT, and m!

  grad = gm exp(i theta) i.e. gm, theta are magnitude and angle of gradient
  kloc = kr exp(i phi)   i.e. kr, phi are magnitude and angle of k-space
  alpha = theta - phi    the angle of the gradient relative to that of k-space
                         (alpha = Pi/2, you go in a circle
                          alpha = 0, you go out radially)

  The variable rad_spacing determines the radial spacing
  in units of the Nyquist distance.
  rad_spacing = 1 gives critical sampling
  rad_spacing > 1 gives undersampling
  rad_spacing can vary throughout spiral generation to create variable density spirals

  KEY EQUATIONS:
  (1) dkr/dphi = rad_spacing*Nyquist/(2 pi)
  (2) dphi/dt = gamma gm Sin(alpha)/kr
  (3) dkr/dt = gamma gm Cos(alpha)

  Solving (1)*(2) = (3) gives
  (4) Tan(alpha) = (2*pi*kr)/(rad_spacing*Nyquist)

  *************************************************************/
  /* Initializations */
  /************************************************************/
  double rast     = GRAST / (double)(subrast);   /* calculation "raster time" in msec */

  double gamma    = spparams[spGAMMA];   /* typically 42.577 kHz/mT */
  double gmax     = spparams[spGMAX];    /* max gradient amplitude in mT/m */
  double slewmax  = spparams[spSLEWMAX]; /* max slew rate, in mT/m/msec*/
  int    gtype    = spparams[spGTYPE];   /* 0 = calculate through readout
                                            1 = include grad ramp-down
                                            2 = include rewinder to end at k=0
                                            3 = include first moment comp */
  double fovxy    = spparams[spFOVXY];   /* enter in m */
  double resxy    = spparams[spRESXY];   /* enter in m : this should be true resolution */
  double fovz     = spparams[spFOVZ];    /* enter in m */ 
  double resz     = spparams[spRESZ];    /* enter in m : this should be true resolution */
  double arms     = spparams[spARMS];    /* number of spiral interleaves*/
  int   sptype    = spparams[spSTYPE];   /* 0 = Archimedean
                                            1 = Cylinder DST 
                                            2 = Spherical DST
                                            3 = Fermat:Floret */
  /* the next 4 variables are for variable density spirals */
  /* they create a transition in the radial spacing as the k-space radius goes from 0 to 1, i.e.*/
  /*    0 < kr < us_0 : spacing = Nyquist distance */
  /* us_0 < kr < us_1 : spacing increases to us_r (affected by ustype)*/
  /* us_1 < kr < 1    : spacing = us_r*/
  int   ustype   = spparams[spUSTYPE]; /* rate of change in undersampling
                                          0 = linear
                                          1 = quadratic
                                          2 = hanning */
  double us_0    = spparams[spUS0];
  double us_1    = spparams[spUS1];
  double us_r    = spparams[spUSR];

  /* For sloppy spirals, this lets us define periodicity in units of iteration loop time */
  /* set this to zero if you do not want sloppy spirals */
  double slop_per = spparams[spSLOP_PER];
  int debug       = spparams[spDEBUG];
  double gxstart  = spparams[spGXSTART];  /* RFS: Gx starting amplitude [mT/m] */
  
  if (debug>0) {
    printf("rast = %f [ms]\t",rast);
    printf("gamma = %f [kHz/mT]\n",gamma);
    printf("gmax = %f [mT/m]\t",gmax);
    printf("slewmax  = %f [mT/m/ms]\n",slewmax);
    printf("gtype = %d (0=readout; 1=0+ramp-down; 2=0+(k->0))\n",gtype);
    printf("fovxy = %f [m]\tresxy = %f [m]\n",fovxy,resxy);
    printf("fovz = %f [m]\tresz = %f [m]\n",fovz,resz);
    printf("arms = %f\n",arms);
    printf("sptype = %d (0=Archimedean; 1=Cylinder DST; 2=Spherical DST; 3=Fermat:Floret)\n",sptype);
    printf("ustype = %d (0=linear; 1=quadratic; 2=hanning)\n",ustype);
    printf("us_0 = %f; us_1 = %f; us_r = %f\n",us_0,us_1,us_r);
    printf("slop_per = %f\n",slop_per);
    printf("gxstart = %f [mT/m]\n",gxstart);    /* RFS  */
  }
  
  double nyquist = (float)(arms)/fovxy; /* radial distance per arm to meet the Nyquist limit*/
  double gamrast = gamma*rast; /* gamrast*g = dk*/
  double dgc     = slewmax*rast; /* the most the gradients can change in 1 raster period*/
  double sub_gamrast = (double)(subrast)*gamrast;
  double sub_dgc     = (double)(subrast)*dgc;

  if (debug>1) {
    printf("nyquist = %f\n",nyquist);
    printf("gamrast = %f\n",gamrast);
    printf("dgc = %f\n",dgc);
    printf("sub_gamrast = %f\n",sub_gamrast);
    printf("sub_dgc = %f\n",sub_dgc);
  }
  
  double *kx    = NULL;
  double *ky    = NULL;
  double *kz    = NULL;
  double *gsign    = NULL;

  double kr, kmx, kmy, kmz, kmr, rnorm;
  double rad_spacing=1;
  double alpha, phi, theta;
  double ux=0,uy=0,uz=0, umag;
  double gx=0,gy=0,gz=0;
  double us_i;
  double gm=0,term;
  double gsum_ramp, gz_sum_ramp;
  double gsum, gsum0, gradtweak, gxsum, gysum, gzsum;
  double krmax, kzmax, krmax2, kzmax2;
  double krlim;
  int i, i0, i1, i_end;
  int j;

  kx    = (double*) malloc(subrast*maxarray*sizeof(double));
  ky    = (double*) malloc(subrast*maxarray*sizeof(double));
  kz    = (double*) malloc(subrast*maxarray*sizeof(double));
  gsign = (double*) malloc(subrast*maxarray*sizeof(double));

  if (kx == NULL || ky == NULL || gsign == NULL) printf ("cant allocate memory\n"); 

  for (i=0;i<subrast*maxarray;i++) gsign[i] = 1.;
  for (i=0;i<subrast*maxarray;i++) kx[i] = 0.;
  for (i=0;i<subrast*maxarray;i++) ky[i] = 0.;
  for (i=0;i<subrast*maxarray;i++) kz[i] = 0.;
  for (i=0;i<maxarray;i++) gxarray[i] = 0.;
  for (i=0;i<maxarray;i++) gyarray[i] = 0.;
  for (i=0;i<maxarray;i++) gzarray[i] = 0.;

  krmax = 0.5/resxy;
  kzmax = 0.5/resz;
  krmax2 = krmax*krmax;
  kzmax2 = kzmax*kzmax;
  krlim = krmax*(1.-(resxy/fovxy));

  /* start out spiral going radially at max slew-rate for 2 time-points */
  kx[0] = 0.;
  ky[0] = 0.;
  kx[1] = gamrast*dgc;
  ky[1] = 0.;
  kx[2] = 3.*gamrast*dgc;
  ky[2] = 0.;
  
  /* RFS: skip slewrate-limited regime */
  if ((gxstart>0.0001) || (gxstart<-0.0001)) {
    if (debug>0) printf("Skipping slewrate-limited regime\n");
    kx[1] = gxstart*gamrast;
    kx[2] = 2*kx[1];
  }
  
  /* IF SPHERE */
  if (sptype == 2) {
    kz[0] = kzmax;
    kz[1] = sqrt(kzmax2*(1.-((kx[1]*kx[1]+ky[1]*ky[1])/krmax2))); /* stay on surface of ellipsoid */
    kz[2] = sqrt(kzmax2*(1.-((kx[2]*kx[2]+ky[2]*ky[2])/krmax2))); /* stay on surface of ellipsoid */
  }

  i = 2;
  kr = kx[2];

  /******************************/
  /* LOOP UNTIL YOU HIT MAX RES */
  /******************************/
  while ((kr <= krlim) && (i < subrast*maxarray-1) ) {

    /**************************/
    /*** STEP 1:  Determine the direction (ux,uy) of the gradient at ~(i+0.5) */
    /**************************/
    /* calculate dk/rast = gamma G*/
    kmx = 1.5*kx[i] - 0.5*kx[i-1];
    kmy = 1.5*ky[i] - 0.5*ky[i-1];
    kmr = sqrt(kmx*kmx + kmy*kmy);

    /***************************
     * Start rad_spacing logic *
     ***************************/
    rnorm = 2.*resxy*kmr; /* the k-space radius, normalized to go from 0 to 1 */
    
    /* determine the undersample factor */
    if (rnorm <= us_0)
      rad_spacing = 1;
    else if (rnorm < us_1) {
      us_i = (rnorm-us_0)/(us_1 - us_0); /* goes from 0 to 1 as rnorm goes from us_0 to us_1*/
      if (ustype == 0) {
        /* linearly changing undersampling*/
        rad_spacing = 1. + (us_r - 1.)*us_i;
      } else if (ustype == 1) {
        /* quadratically changing undersampling*/
        rad_spacing = 1. + (us_r - 1.)*us_i*us_i;
        } else if (ustype == 2) {
          /* Hanning-type change in undersampling */
          rad_spacing = 1. + (us_r - 1.)*0.5*(1.-cos(us_i*M_PI));
        }
      } /* if (rnorm < us_1) */
      else {
        rad_spacing = us_r;
      } /* rnorm > us_1 */

    /* Undersample spiral for Spherical-Distributed Spiral */
    if (sptype == 2) {
      if(rnorm < 1.0)
        rad_spacing = MIN(fovz/resz, rad_spacing/sqrt(1.0 - (rnorm*rnorm)));
      else
        rad_spacing = fovz/resz;
    } /* SDST */

    /* MAKE FERMAT SPIRAL FOR FLORET*/
    if (sptype == 3 && rnorm > 0.) rad_spacing *= 1./rnorm;

    /* Sloppy Spirals - add variability to rad_spacing for reduced aliasing coherence */ 
    /* A couple different options here are commented out */
    /* Lots of ways to be sloppy */
    if (slop_per > 0) {
      /* rad_spacing = MAX(1., (rad_spacing + ((rad_spacing-1.)*sin(2.*M_PI*(double)(i)/slop_per))));
       rad_spacing += (rad_spacing-1.)*sin(2.*M_PI*slop_per*atan2(ky[i],kx[i])); */
      rad_spacing += (rad_spacing-1.)*sin(2.*M_PI*slop_per*rnorm);
    }

    /*************************
     * End rad_spacing logic *
     *************************/

    /* See the Key Equation 4 at the beginning of the code */
    alpha = atan(2.*M_PI*kmr/(rad_spacing*nyquist));
    phi = atan2(kmy,kmx);
    theta = phi + alpha;

    ux = cos(theta);
    uy = sin(theta);

    /* IF SPHERICAL DST
     * u dot km is zero if moving on a sphere (km is radial, u is tangential,
     * thus km stays on the sphere)
     * We are on an ellipsoid, but can normalize u and km by krmax and kzmax to make this work
     * The final gradient vector (ux uy uz) will be tangential to the sphere
     */
    if (sptype == 2) {
      kmz = 1.5*kz[i] - 0.5*kz[i-1];
      uz = -((ux*kmx + uy*kmy)/krmax2)*(kzmax2/kmz);
      umag = sqrt(ux*ux + uy*uy + uz*uz);
      ux = ux/umag;
      uy = uy/umag;
      uz = uz/umag;
      gz = (kz[i] - kz[i-1])/gamrast;
    }

    /**************************/
    /*** STEP 2: Find largest gradient magnitude with available slew */
    /**************************/

    /* Current gradient*/
    gx = (kx[i] - kx[i-1])/gamrast;
    gy = (ky[i] - ky[i-1])/gamrast;

    /*
     * solve for gm using the quadratic equation |gm u - g| = dgc
     * which is
     * (gm u - g)(gm u* - g*) = dgc^2
     * which gives
     * gm^2 (u u*) - gm (g u* + u g*) + g g* - dgc^2 = 0
     *
     * Replacing u u* with 1 (i.e. u is a unit vector) and
     * replacing (g u* + u g*) with 2 Real[g u*]
     * this is
     * gm^2 + gm (2 b) + c = 0
     * giving
     * gm = -b +/- Sqrt(b^2 - c)
     * The variable "term" = (b^2 - c) will be positive if we can meet the 
     * desired new gradient
     */
    term = dgc*dgc - (gx*gx + gy*gy + gz*gz) + (ux*gx + uy*gy + uz*gz)*(ux*gx + uy*gy + uz*gz);

    if (term >= 0) {
        /* Slew constraint is met! Now assign next gradient and then next k value
       NOTE gsign is +1 or -1
         if gsign is positive, we are using slew to speed up (increase gm) 
         * as much as possible
         if gsign is negative, we are using slew to slow down (decrease gm) 
         * as much as possible
         */
      gm  = MIN((ux*gx + uy*gy + uz*gz) + gsign[i]*sqrt(term),gmax);
      gx = gm*ux;
      gy = gm*uy;

      kx[i+1] = kx[i] + gx*gamrast;
      ky[i+1] = ky[i] + gy*gamrast;

      /* If SPHERE */
      if (sptype == 2)
        kz[i+1] = sqrt(kzmax2*(1.-((kx[i+1]*kx[i+1]+ky[i+1]*ky[i+1])/krmax2))); 
      /* stay on surface of ellipsoid */

      i++;
    } /* term >= 0 */
    else {
      /* We can't go further without violating the slew rate
       * This means that we've sped up too fast to turn here at the desired 
       * curvature. We are going to iteratively go back in time and slow 
       * down, rather than speed up, at max slew. Here we'll keep looking 
       * back until gsign is positive, then add another negative gsign, 
       * just far enough to make the current corner
       */
      while ((i>3) && (gsign[i-1] == -1)) {
        i--;
      }
      gsign[i-1] = -1;
      i = i-2;
      if ((i<2) && ((gxstart>0.0001)||(gxstart<-0.0001))) {  /* RFS */
        printf("Error: cannot decrease i(=%d) further\n",i);
        return;
      }
    } /* term < 0 */

    kr = sqrt(kx[i]*kx[i] + ky[i]*ky[i]);

  } /* MAIN kr loop */

  i_end = i;

  /********************************************
   * DONE LOOPING FOR SAMPLING PORTION
   * recast k to g while subsampling by subrast
   ********************************************/
  gxarray[0] = 0.;
  gyarray[0] = 0.; 
  gzarray[0] = 0.; 
  gxsum = 0.;
  gysum = 0.;
  gzsum = 0.;
  
  /* RFS: skip slewrate-limited regime */
  if ((gxstart>0.0001) || (gxstart<-0.0001)) {
    gxarray[0] = gxstart;
    gxsum = gxstart;
  }
  
  for (j=1;j<=(i_end/subrast);j++) {
    i1 = j*subrast;
    i0 = (j-1)*subrast;
    gxarray[j] = (kx[i1]-kx[i0])/sub_gamrast;
    gyarray[j] = (ky[i1]-ky[i0])/sub_gamrast;
    gzarray[j] = (kz[i1]-kz[i0])/sub_gamrast;
    gxsum = gxsum + gxarray[j];
    gysum = gysum + gyarray[j];
    gzsum = gzsum + gzarray[j];
  }
  (*spgrad_na) = j;

  /* recalculate these ending gradient points */
  gm = sqrt(gxarray[(*spgrad_na)-1]*gxarray[(*spgrad_na)-1] +
            gyarray[(*spgrad_na)-1]*gyarray[(*spgrad_na)-1] +
            gzarray[(*spgrad_na)-1]*gzarray[(*spgrad_na)-1]);
  ux = gxarray[(*spgrad_na)-1]/gm;
  uy = gyarray[(*spgrad_na)-1]/gm;
  uz = gzarray[(*spgrad_na)-1]/gm;

  /**************************************************
   * NOW, if requested via gtype, go to g=0 and k=0
   * I've tried other ways to be faster, can't find them
   **************************************************/
  
  /* first we'll ramp gradients to zero */
  /* note {ux,uy} is still pointing in the gradient direction */
  if (gtype > 0) {
    gz_sum_ramp = 0;
    while ((gm > 0) && (j < maxarray-1)) {
      gm = MAX(0,gm - sub_dgc);
      gxarray[j] = gm*ux;
      gyarray[j] = gm*uy;
      gzarray[j] = gm*uz;
      gxsum = gxsum + gxarray[j];
      gysum = gysum + gyarray[j];
      gzsum = gzsum + gzarray[j];
      gz_sum_ramp += gzarray[j];
      j++;
    }
  }
  *spgrad_nb = j;

  /* now point gradient towards the k-space origin */
  /* {ux,uy} will be a unit vector in that direction */
  if (gtype > 1) {

    /* NOTE: spherical needs a prephaser not a rewinder 
     * so just rewind x and y in that case */
    gsum = sqrt(gxsum*gxsum + gysum*gysum + gzsum*gzsum);
    if (sptype == 2 ) gsum = sqrt(gxsum*gxsum + gysum*gysum + gz_sum_ramp*gz_sum_ramp);
    gsum0 = gsum;
    ux = -gxsum/gsum;
    uy = -gysum/gsum;
    uz = -gzsum/gsum;
    if (sptype == 2) uz = -gz_sum_ramp/gsum;
    gsum_ramp = 0.5*gm*(gm/sub_dgc); 
    /* this is *roughly* how much the area changes if we ramp down the gradient NOW*/
    /* this value is zero right now (gm = 0), but it will make sense below */
    
    /* increase gm while we can */
    while ((gsum_ramp < gsum) && (j < maxarray-1)) {
      gm = MIN(gmax,gm+sub_dgc);
      gxarray[j] = gm*ux;
      gyarray[j] = gm*uy;
      gzarray[j] = gm*uz;
      gsum = gsum - gm;
      j++;
      gsum_ramp = 0.5*gm*(gm/sub_dgc); 
      /* see - now this makes sense; this tells us when to start ramping down */
    }

    /* We've overshot it by a tiny bit, but we'll fix that later */
    /* Ramp down for now */
    while ((gm > 0) && (j < maxarray-1)) {
      gm = MAX(0,gm-sub_dgc);
      gxarray[j] = gm*ux;
      gyarray[j] = gm*uy;
      gzarray[j] = gm*uz;
      gsum = gsum - gm;
      j++;
    }
    *spgrad_nc = j;

    /* OK - gm is zero, but gsum is probably not EXACTLY zero. Now scale 
     * the rewinder to make the sum exactly zero */
    gradtweak = gsum0/(gsum0-gsum);
    for (j=(*spgrad_nb); j<(*spgrad_nc); j++) {
      gxarray[j] = gradtweak*gxarray[j];
      gyarray[j] = gradtweak*gyarray[j];
      gzarray[j] = gradtweak*gzarray[j];
    }
  }
  free (kx);
  free (ky);
  free (kz);
  free (gsign);
}  /* main function bnispiralgen */


/* *************************** */
/* RFS: mex function interface */
void mexFunction(int nlhs, mxArray *plhs[], 
        int nrhs, const mxArray *prhs[]) {
    double spparams[spARRSIZE];
    double *fov=NULL;        /* fov (xy,z) [m] */
    double *mtx=NULL;        /* matrix size (xy,z) */
    int sptype=0;            /* spiral type */
    double arms=1.0;         /* #spiral interleaves */
    double gmax=33;          /* maximum gradient strength [mT/m] */
    double smax=120;         /* max. slewrate [T/m/s] */
    double gamma=42.577469;  /* gyromagnetic ration [kHz/mT] */
    double *varden=NULL;     /* variable density */
    int gtype=1;             /* add rewinder */
    double slop_per=0;       /* sloppy spiral */
    int maxarray=32766;      /* max. gradient size */
    int debug=0;             /* print info */
    double gxstart=0;        /* RFS: gx starting amplitude */
    double *gxarray = NULL;  /* Gx */
    double *gyarray = NULL;  /* Gy */
    double *gzarray = NULL;  /* Gz */
    double *grad = NULL;     /* full gradient matrix */
    int spgrad_na=0;     /* #gradient points to reach the end of k-space */
    int spgrad_nb=0;         /* spgrad_na + #gradient points to ramp G to zero */
    int spgrad_nc=0;         /* spgrad_nb + #gradient points to rewind k to zero */
    int spgrad_nd=0;         /* unused -> not actually implemented */
    int ng_tot=0;            /* total #gradient points */
    int ndim = 2;            /* #spatial dimensions (only sptype==2 is 3D) */
    double *ng = NULL;       /* spgrad_na,ng_tot */
    int l;                   /* loop count */
    
    if (nrhs<2) {
        printf("BNISPIRALGEN  Spiral Generation; Jim Pipe\n");
        printf("Calculate single spiral readout arm\n");
        printf("Usage:\n");
        printf("\t[grad,ng]=bnispiralgen(fov,mtx,sptype,arms,gmax,smax,...\n");
        printf("\t\tvarden,gtype,slop_per,gamma,maxarray,debug,gxstart);\n");
        printf("\tfov      Field of view (xy,z)[m]\n");
        printf("\tmtx      Matrix size (xy,z)\n");
        printf("\tsptype   0=Archimedean; 1=Cylinder DST; 2=Spherical DST; 3=Fermat:Floret (0)\n");
        printf("\tarms     #spiral interleaves (1)\n");
        printf("\tgmax     Max. gradient strength [mT/m] (33)\n");
        printf("\tsmax     Max. slewrate [T/m/s] (120)\n");
        printf("\tvarden   Variable density (ustype,us_0,us_1,us_r) (0,1,1,1) \n");
        printf("\tgtype    Add rewinder:  0=off; 1=ramp-down;  2=return to k=0 (1)\n");
        printf("\tslop_per Sloppy spiral (0)\n");
        printf("\tgamma    Gyrogmagnetic ratio [kHz/mT] (42.577)\n");
        printf("\tmaxarray max. #grad points\n");
        printf("\tdebug    Print debug information (0)\n");
        printf("\tgxstart  Gx starting amplitude (0) [mT/m]; skip slewrate-limited regime\n");
        printf("\tgrad     Spiral gradient waveform [mT]\n");
        printf("\tng       #waveform points: 1=spiral; 2=spiral+rewinder\n");
        printf("Source: www.ismrm.org/mri_unbound/sequence.htm\n");
        printf("\t-> spiralgen_jgp_12oct\n");
        return;
    }

    /* input variables */
    if ((mxGetM(prhs[0])>1) || (mxGetN(prhs[0])>1)) {
        fov = mxGetPr(prhs[0]);
    } else {
        fov = (double*) malloc(2*sizeof(double));
        fov[0] = mxGetScalar(prhs[0]);
        fov[1] = fov[0];
    }
    if ((mxGetM(prhs[1])>1) || (mxGetN(prhs[1])>1)) {
        mtx = mxGetPr(prhs[1]);
    } else {
        mtx = (double*) malloc(2*sizeof(double));
        mtx[0] = mxGetScalar(prhs[1]);
        mtx[1] = mtx[0];
    }
    
    if ((nrhs>2) && (mxGetM(prhs[2])>0)) sptype = mxGetScalar(prhs[2]);
    if ((nrhs>3) && (mxGetM(prhs[3])>0)) arms = mxGetScalar(prhs[3]);
    if ((nrhs>4) && (mxGetM(prhs[4])>0)) gmax = mxGetScalar(prhs[4]);
    if ((nrhs>5) && (mxGetM(prhs[5])>0)) smax = mxGetScalar(prhs[5]);
    if ((nrhs>6) && (mxGetM(prhs[6])>0)) {
        if ((mxGetM(prhs[6])==4) || (mxGetN(prhs[6])==4)) {
            varden = mxGetPr(prhs[6]);
        } else {
            printf("Error: length(varden)!=4\n");
            return;
        }
    } else {  /* assign default values */
        varden = (double*) malloc(4*sizeof(double));
        varden[0]=0; varden[1]=1; varden[2]=1; varden[3]=1;
    }
    if ((nrhs>7) && (mxGetM(prhs[7])>0)) gtype =  mxGetScalar(prhs[7]);
    if ((nrhs>8) && (mxGetM(prhs[8])>0)) slop_per = mxGetScalar(prhs[8]);
    if ((nrhs>9) && (mxGetM(prhs[9])>0)) gamma =  mxGetScalar(prhs[9]);
    if ((nrhs>10)&& (mxGetM(prhs[10])>0)) maxarray = mxGetScalar(prhs[10]);
    if ((nrhs>11)&& (mxGetM(prhs[11])>0)) debug = mxGetScalar(prhs[11]);
    if ((nrhs>12)&& (mxGetM(prhs[12])>0)) gxstart =  mxGetScalar(prhs[12]);

    /* parameter checks */
    if (arms<=0) {
        printf("Error: arms (=%f)<=0\n",arms);
        return;
    }
    
    /* assign to spparams vector */
    spparams[spGAMMA] = gamma;
    spparams[spGMAX] = gmax; 
    spparams[spSLEWMAX] = smax;
    spparams[spGTYPE] = gtype;
    spparams[spFOVXY] = fov[0];
    spparams[spRESXY] = fov[0]/(double)mtx[0];
    spparams[spFOVZ] = fov[1];
    spparams[spRESZ] = fov[1]/(double)mtx[1];
    spparams[spARMS] = arms;
    spparams[spSTYPE] = sptype;
    spparams[spUSTYPE] = varden[0];
    spparams[spUS0] = varden[1];
    spparams[spUS1] = varden[2];
    spparams[spUSR] = varden[3];
    spparams[spSLOP_PER] = slop_per;
    spparams[spDEBUG] = debug;
    spparams[spGXSTART] = gxstart;     /* RFS */
    
        
    /* initialise memory for gradient arrays */
    gxarray = (double*) malloc(maxarray*sizeof(double));
    gyarray = (double*) malloc(maxarray*sizeof(double));
    gzarray = (double*) malloc(maxarray*sizeof(double));
    
    /* actual execution of main program */
    bnispiralgen(spparams, maxarray, gxarray, gyarray, gzarray, 
                   &spgrad_na, &spgrad_nb, &spgrad_nc, &spgrad_nd);
    
    /* #grad points */
    ng_tot = spgrad_na;
    if (spparams[spGTYPE]>0) ng_tot = spgrad_nb;
    if (spparams[spGTYPE]>1) ng_tot = spgrad_nc;
    if (debug>0) {
        printf("spgrad = %d %d %d; ",spgrad_na,spgrad_nb,spgrad_nc);
        printf("ng_tot = %d\n",ng_tot);
    }
    if (ng_tot==maxarray) {
        printf("Warning: max. #gradient points (=%d) exceeded\n",maxarray);
    }
    
    /* #spatial dimension */
    if (sptype==2) ndim = 3;
    
    /* copy waveform to output matrix */
    plhs[0] = mxCreateDoubleMatrix(ng_tot, ndim, mxREAL);
    grad = (double*)mxGetData(plhs[0]);
    for (l=0; l<ng_tot; l++) {
        grad[l] = gxarray[l];
        grad[l + ng_tot] = gyarray[l];
        if (ndim>2) grad[l+2*ng_tot] = gzarray[l];
    }
    
    /* return #gradient points: pure spiral, incl. rewinder */
    plhs[1] = mxCreateDoubleMatrix(1, 2, mxREAL);
    ng = (double*)mxGetData(plhs[1]);
    ng[0] = (double)spgrad_na;
    ng[1] = (double)ng_tot;
}

/* undo common macro names */
#ifdef UNDEFMAX
#undef MAX
#endif

#ifdef UNDEFMIN
#undef MIN
#endif

