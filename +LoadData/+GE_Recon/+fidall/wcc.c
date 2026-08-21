#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "mex.h"

#define PI 3.14159265358979
#define max(a,b) ((a) < (b) ? (b) : (a)) 
#define min(a,b) ((a) < (b) ? (a) : (b)) 


int rootof(double a, double b, double c, double *s1, double *s2)  {
double descrim;
if (a == 0) {
  if (b == 0) {
    *s1 = 0;
    *s2 = 0;
    return 0;
  } else {
    *s1 = -c/b;
    *s2 = -c/b;
    return 1;
  }
} else {
  descrim = b*b-4*a*c;
  if (descrim<0) {
    *s1 = 0;
    *s2 = 0;
    return 0;
  } else {
    *s1 = (-b + sqrt(descrim))/2/a;
    *s2 = (-b - sqrt(descrim))/2/a;
    return 2;
  }
}
}


int circcirc(double d, double e, double k, double s, double type, double *x, double *y) {

double s1;
double s2;
double t;
double ade;
double tx1,tx2,ty1,ty2;
double a,b,c;
int numsols;
if (e == 0)
    ade = 100;
else
    ade = abs(d/e);

if ((d == 0) && (e==0)) {
    *x = 0;
    *y = 0;
    return 0;
} else if (ade<1) {
    t = e;
    e = d;
    d = t;
    a = 4*d*d+4*e*e;
    b = (-4*e*e*e-4*k*k*e+4*e*s*s-4*d*d*e);
    c = d*d*d*d+k*k*k*k-2*e*e*s*s+2*d*d*e*e-2*d*d*s*s+s*s*s*s-2*k*k*s*s-2*k*k*d*d+e*e*e*e+2*k*k*e*e;
    numsols = rootof(a,b,c,&s1,&s2);
    if (numsols>0) {
        tx1 = s1;
        ty1 = -0.5*(2*e*s1-k*k-d*d-e*e+s*s)/d;
        tx2 = s2;
        ty2 = -0.5*(2*e*s2-k*k-d*d-e*e+s*s)/d;
    } else {
        *y = 0;
        *x = 0;
        return 0;
    }  
    t = e;
    e = d;
    d = t;
} else {
    a = 4*d*d+4*e*e;
    b = (-4*e*e*e-4*k*k*e+4*e*s*s-4*d*d*e);
    c = d*d*d*d+k*k*k*k-2*e*e*s*s+2*d*d*e*e-2*d*d*s*s+s*s*s*s-2*k*k*s*s-2*k*k*d*d+e*e*e*e+2*k*k*e*e;
    numsols = rootof(a,b,c,&s1,&s2);
    if (numsols>0) {
        ty1 = s1;
        tx1 = -0.5*(2*e*s1-k*k-d*d-e*e+s*s)/d;
        ty2 = s2;
        tx2 = -0.5*(2*e*s2-k*k-d*d-e*e+s*s)/d;
    } else {
        *y = 0;
        *x = 0;
        return 0;
    }
}


if ((type*(-tx1*e+ty1*d))<0) {
    *x = tx2;
    *y = ty2;
} else {
    *x = tx1;
    *y = ty1;
}
    return numsols;
}


int intlinecirc(double a, double b, double d, double e, double r, double *s1, double *s2) {
return rootof((a*a+b*b),(-2*e*b-2*d*a),(d*d+e*e-r*r),s1,s2);

}

double anglebetween_new(double ar, double ai, double br, double bi) {
	double rr,ri;
	rr = ar*br+ai*bi;
	ri = -ar*bi+ai*br;
	return atan2(ri,rr);
}

double anglebetween(double ar, double ai, double br, double bi) {
	double rr,ri;
	rr = -ar*br-ai*bi;
	ri = ar*bi-ai*br;
	return atan2(ri,rr);
}

double __inline__ calcn(double k, double sintheta, double FOV1, double FOV2, double NINT, double kmax, double dcf, int nosmooth) {
  double Gtwist;
  double n;
  double c1,c2,kstart,ff,a,nt1,nt2,cc;
  if (dcf<0) {
    Gtwist = 2*PI*sintheta*k*FOV2/NINT;
    n = 2*PI*FOV2/NINT; 
  } else if (dcf>0.5) {
    Gtwist = pow(2*PI*sintheta*k*FOV2/NINT*k/kmax,2)-FOV2*FOV2/(FOV1*FOV1);
    Gtwist = sqrt((Gtwist>0)?Gtwist:0);
    n = Gtwist/(k*sintheta);
    /* n = (n>0.00001)?n:0.00001; */
  } else if (dcf>0) {
	  c1 = 4*PI*PI*sintheta*sintheta*FOV2*FOV2/(NINT*NINT);
	  c2 = (FOV2*FOV2)/(FOV1*FOV1);
	  kstart = sqrt(c1*c2)/c1; 
	  if (kstart<0.5*kmax) 
	     kstart = 2*sqrt(c1*c2)*kmax*kmax/(c1*kmax*kmax+4*c2); 
	  cc = max(0.0001,min(0.5*kmax-kstart,max(0.0001,(kmax/20.0))));
	  if (k<0.5*kmax-2*cc) {
		  ff = 1;
	  } else if (k<0.5*kmax-cc) {
		  ff = 1-(k-(0.5*kmax-2*cc))*(k-(0.5*kmax-2*cc))/(2.0*cc*cc); 
		  if (nosmooth)
			  ff = 1;
	  } else if (k<0.5*kmax) {
		  ff = (k-(0.5*kmax))*(k-(0.5*kmax))/(2.0*cc*cc); 
		  if (nosmooth)
			  ff = 1;
	  } else {
                  ff = 0;
		  if (nosmooth)
			  ff = 1;
          }
	  nt1 = sqrt(max(0,1-4*(k*k/(kmax*kmax))));
	  nt2 = 1+nt1*ff;
	  n = c1*k*k/(nt2*nt2)-c2;
	  n = max(n,0);
	  n = sqrt(n)/(k*sintheta);
  } else {
    Gtwist = pow(2*PI*sintheta*k*FOV2/NINT,2)-FOV2*FOV2/(FOV1*FOV1);
    Gtwist = sqrt((Gtwist>0)?Gtwist:0);
    n = Gtwist/(k*sintheta);
    /* n = (n>0.00001)?n:0.00001; */
  }
  return n;
}

double __inline__ calcdn(double k, double sintheta, double FOV1, double FOV2, double NINT, double kmax, double dcf, int nosmooth) {
  double Gtwist;
  double dn;
  double c1,c2,kstart,ff,dff,a,nt1,nt2,nt3,cc;
  if (dcf<0) {
    dn = 0;
  } else if (dcf>0.5) {
     Gtwist = pow(2*PI*sintheta*k*FOV2/NINT*k/kmax,2)-pow(FOV2/FOV1,2);
     Gtwist = sqrt((Gtwist>0)?Gtwist:0);
     dn = 8*k*k*sintheta*PI*PI*FOV2*FOV2/(Gtwist*NINT*NINT*kmax*kmax)-Gtwist/(k*k*sintheta);
  } else if (dcf>0)  {
     c1 = 4*PI*PI*sintheta*sintheta*FOV2*FOV2/(NINT*NINT);
     c2 = (FOV2*FOV2)/(FOV1*FOV1);
	  kstart = sqrt(c1*c2)/c1; 
	  if (kstart<0.5*kmax) 
	    kstart = 2*sqrt(c1*c2)*kmax*kmax/(c1*kmax*kmax+4*c2); 
	  cc = max(0.0001,min(0.5*kmax-kstart,max(0.0001,(kmax/20.0))));
	  if (k<0.5*kmax-2*cc) {
		  ff = 1;
		  dff = 0;
	  } else if (k<0.5*kmax-cc) {
		  ff = 1-(k-(0.5*kmax-2*cc))*(k-(0.5*kmax-2*cc))/(2.0*cc*cc); 
		  dff = -(k-0.5*kmax+2*cc)/(cc*cc);
		  if (nosmooth) {
			  ff = 1;
			  dff = 0;
		  }
	  } else if (k<0.5*kmax) {
		  ff = (k-(0.5*kmax))*(k-(0.5*kmax))/(2.0*cc*cc); 
		  dff = (k-0.5*kmax)/(cc*cc);
		  if (nosmooth) {
			  ff = 1;
			  dff = 0;
		  }
	  } else {
                  ff = 0;
		  dff = 0;
		  if (nosmooth) {
			  ff = 0;
			  dff = 0;
		  }
          }
	  nt1 = sqrt(max(1e-8,1-4*(k*k/(kmax*kmax))));
	  nt2 = 1+nt1*ff;
	  nt3 = 4*k*ff/(nt1*kmax*kmax); 
          
         dn = (c1*k/(nt2*nt2) - c1*k*k*(-4*k*ff/(nt1*kmax*kmax)+nt1*dff)/(nt2*nt2*nt2))
		 / (sqrt(c1*k*k/(nt2*nt2)-c2)*k*sintheta)
	      - (sqrt(c1*k*k/(nt2*nt2)-c2)/(k*k*sintheta));

/*
	  n = c1*k*k/(nt2*nt2)-c2;
	  n = max(n,0);
	  n = sqrt(n)/(k*sintheta);
	  
	  nt1 = sqrt(1-4*(ff*ff/(kmax*kmax)));
	  nt2 = 1+nt1;
	  nt3 = c1*k*k/(nt2*nt2)-c2;
	  nt3 = max(nt3,0);
	  nt3 = sqrt(nt3);
	  dn  = 0.5*(2*c1*k/(nt2*nt2)+8*c1*k*k*ff*dff/(nt2*nt2*nt2*((nt1==0)?1:nt1)*kmax*kmax))/(nt3*k*sintheta)-(nt3/(k*k*sintheta));
	  */
  } else {
     Gtwist = pow(2*PI*sintheta*k*FOV2/NINT,2)-pow(FOV2/FOV1,2);
     Gtwist = sqrt((Gtwist>0)?Gtwist:0);
     dn = 4*sintheta*PI*PI*FOV2*FOV2/(Gtwist*NINT*NINT)-Gtwist/(k*k*sintheta);
  }
  return (dn); 

}

int rto(double Gc, double kdes, double Smax, double Gmax, double gam, double Ts, double *g, double *k) {
double ta,tb,tc,tba,ttot, karea, dkarea, ddkarea;
int numc,numb,numa,numd,ntot;
double gend,kreq,sumg,greq;
int ai;

  /* ta:  ramp from Gc down to 0. */
  /* tb:  triangle above Gc */
  /* tc:  plateau at Gmax */

  ta = Gc/(Smax/Ts);
  karea = gam/2*(Smax/Ts)*ta*ta; /* Area achievable with an immediate ramp from Gc down to 0*/
  dkarea = kdes-karea;
  if (dkarea<0) {
	  return -1;
  }
  tb = (-gam*Gc+sqrt((gam*gam*Gc*Gc)-(gam*Smax/Ts)*(-dkarea)))/(gam*Smax/Ts)*2;
  if (((ta+tb/2)*Smax/Ts)>Gmax) {
	  tb = 2*(Gmax/(Smax/Ts)-ta);
          ddkarea = kdes-(karea+gam*Smax/Ts*tb*tb/4.0+gam*Gc*tb);
	  tc = ddkarea/(gam*Gmax);
  } else {
	  tc = 0;
  }

  tba = ta+tb/2;
  double tb_leftover;
  double tc_leftover;

  numc = tc/Ts;
  if (g!=0) {
  k[0] = kdes/gam;
  /* Start in the k-domain */
  for (ai = 1; ai<=floor(tb/Ts/2); ai++) {
	  k[ai] = k[0] - Gc*ai*Ts - (ai*Ts*ai*Ts)/2*(Smax/Ts); 
  }
  tb_leftover = tb/Ts/2-floor(tb/Ts/2);

  if ((tc/Ts+tb_leftover)<1) {
	  tc_leftover = tc/Ts+tb_leftover;
  } else { 
    for (ai = floor(tb/Ts/2)+1; ai<=floor(tb/Ts/2+tc/Ts); ai++) {
	    k[ai] = k[0] - Gc*tb/2 - tb*tb/8*Smax/Ts - (ai*Ts - tb/2)*Gmax;
    }
    tc_leftover = tb/Ts/2+tc/Ts-floor(tb/Ts/2+tc/Ts);
  }

  if ((tc_leftover+tba/Ts)<1) {
	  k[(int)ceil(tb/Ts/2+tc/Ts+tba/Ts)] = 0;
  } else {
    for (ai = floor(tb/Ts/2+tc/Ts)+1; ai<=floor(tb/Ts/2+tc/Ts+tba/Ts); ai++) {
	    k[ai] = k[0] - Gc*tb/2 - tb*tb/8*Smax/Ts - (tc*Gmax) - (ai*Ts - tb/2 - tc)*(Gc+tb/2*Smax/Ts) + (ai*Ts-tb/2-tc)*(ai*Ts-tb/2-tc)/2*Smax/Ts;
    }
    k[(int)ceil(tb/Ts/2+tc/Ts+tba/Ts)] = 0;
  }

  for (ai = ceil(tb/Ts/2+tc/Ts+tba/Ts); ai>=1; ai--) {
     g[(int)ceil(tb/Ts/2+tc/Ts+tba/Ts)-ai+1] = (k[ai-1]-k[ai])/Ts;
  }
     g[0] = 0;
}

  return ceil(tb/Ts/2+tc/Ts+tba/Ts)+1;

}


int rtz(double gcx, double gcy, double gcz, double kcx, double kcy, double kcz, double theta, int maxai, double Gmax_xy, double Gmax_z, double Smax_xy, double Smax_z, double Ts, double * pgx, double *pgy, double *pgz, double *pkx, double *pky, double *pkz, double FOV1, double FOV2, double NINT, double kmax, double dcf) {

  double dang,dangold;
  double kxy, gxy;
  double kxymin, kxymax;
  int ai;
  int bi;
  int done;
  int fail;
  int lenai;
  int ttot;
  double s1,s2;
  double s1min, s2min;
  double s1max, s2max;
  double s1l, s2l;
  double Smax_c;
  double Gmax_c;
      double kr_next;
      double kr_cur;
      double kr_dangle;
      double n_des;
      double n_act;
  int numsols,numsolsmin,numsoll;
  double t1,tmin,tmax,tl;
  int numsolsline, numsolsmax;
      double range_min;
      double range_max;
      double rang;
      double kang,tt,tkx,tky,tkz;
      double s1kxy,s2kxy;
      int numsolskxy;
      double t1xy,t2xy;
      int numsolstmax;
      double s1tmax, s2tmax;
      double ttmax;
      int nosmooth;
      nosmooth = 1;
  pgx[0] = gcx;
  pgy[0] = gcy;
  pgz[0] = gcz;
  pkx[0] = kcx/4258/Ts;
  pky[0] = kcy/4258/Ts;
  pkz[0] = kcz/4258/Ts;
  ai = 0;
  done = 0;
  fail = 0;
  lenai = 1;
  kxy = sqrt(pkx[ai]*pkx[ai]+pky[ai]*pky[ai]);
  gxy = sqrt(pgx[ai]*pgx[ai]+pgy[ai]*pgy[ai]);
  dang = anglebetween(pgx[ai],pgy[ai],pkx[ai],pky[ai]); 
  while ((done == 0) && (fail == 0) && (ai<(maxai))) {

      kxymin = (pkz[ai]+max(-Gmax_z,(pgz[ai]-Smax_z)))*tan(theta);
      kxymax = (pkz[ai]+min(Gmax_z,(pgz[ai]+Smax_z)))*tan(theta);


      /* Find intersection of SMAX with Gc */
      numsols = circcirc(pgx[ai],pgy[ai],gxy,Smax_xy,-1,&s1,&s2);
      pgx[ai+1] = s1;
      pgy[ai+1] = s2;
      s1 = s1+pkx[ai];
      s2 = s2+pky[ai];
      /* Find intersection of Gc with xy min (Z slew and amp constraints) */
      numsolsmin = circcirc(pkx[ai],pky[ai],kxymin,gxy,-1,&s1min,&s2min);
      numsolsmax = circcirc(pkx[ai],pky[ai],kxymax,gxy,-1,&s1max,&s2max);
      /* Find intersection of Gc with line back to origin*/
      numsolsline = intlinecirc(pkx[ai]/kxy,pky[ai]/kxy,pkx[ai],pky[ai],gxy,&s1l,&s2l);
      s1l = s2l*pkx[ai]/kxy;
      s2l = s2l*pky[ai]/kxy;

      /* Find intersection of kxymin and SMAX */
      numsolskxy = circcirc(pkx[ai]+pgx[ai],pky[ai]+pgy[ai],kxymin,Smax_xy,1,&s1kxy,&s2kxy);

      /* Find intersection of kxy and Gc */
      numsolstmax = circcirc(pkx[ai],pky[ai],kxy,gxy,-1,&s1tmax,&s2tmax);
      ttmax = anglebetween_new(s1tmax-pkx[ai],s2tmax-pky[ai],-pkx[ai],-pky[ai]);
      /* Convert to parametric form */
      if (numsols>0) {
      t1 = anglebetween_new(s1-pkx[ai],s2-pky[ai],-pkx[ai],-pky[ai]);
      } else {
      t1 = 0;
      }

      if (numsolskxy>0) {
      t1xy = anglebetween_new(pkx[ai],pky[ai],s1kxy,s2kxy);
      } else {
      t1xy = 0;
      }
      if (numsolsmin>0) {
      t2xy = anglebetween_new(pkx[ai],pky[ai],s1min,s2min);
      } else {
      t2xy = 0;
      }

      if (t1xy<0) {
      t1xy = 0; 
      }

      tmin = anglebetween_new(s1min-pkx[ai],s2min-pky[ai],-pkx[ai],-pky[ai]);
      tmax = anglebetween_new(s1max-pkx[ai],s2max-pky[ai],-pkx[ai],-pky[ai]);
      tl = anglebetween_new(pgx[ai],pgy[ai],-pkx[ai],-pky[ai]);

      range_min = max(0,t1);
      range_min = (numsolsmin>0)?max(tmin,range_min):range_min; 
      range_max = ((numsolsmax>0)&&(tmax>0))?min(tmax,tl):tl; 
      range_max = min(range_max,ttmax);

      

      kang = atan2(-pky[ai],-pkx[ai]); 

      /* Ensure a point in the range exists. If not, fail!*/
      if (range_max-range_min>=0) {
      
      /* 1. Check at range_max (ensure that at least one point has n_act>=n_des.  If not, fail!) */
         tt = range_max;
         tkx = pkx[ai]+gxy*(cos(kang)*cos(tt)-sin(kang)*sin(tt));
         tky = pky[ai]+gxy*(cos(kang)*sin(tt)+sin(kang)*cos(tt));
	 tkz = sqrt(tkx*tkx+tky*tky)/tan(theta);
         kr_cur =  4258*Ts*sqrt(pkx[ai]*pkx[ai]+pky[ai]*pky[ai]+pkz[ai]*pkz[ai]);
         kr_next =  4258*Ts*sqrt(tkx*tkx+tky*tky+tkz*tkz);
         n_des = calcn(kr_next, sin(theta), FOV1, FOV2, NINT, kmax, dcf,nosmooth); 
	 n_act = anglebetween_new(pkx[ai],pky[ai],tkx,tky)/(kr_cur-kr_next);
	 pkx[ai+1] = tkx;
	 pky[ai+1] = tky;
	 pkz[ai+1] = tkz;

	 if (n_act>=n_des) {
      /* 2. Check at range_min (if range_min has n_act>=n_des, don't need to search) */
         tt = range_min;
         tkx = pkx[ai]+gxy*(cos(kang)*cos(tt)-sin(kang)*sin(tt));
         tky = pky[ai]+gxy*(cos(kang)*sin(tt)+sin(kang)*cos(tt));
	 tkz = sqrt(tkx*tkx+tky*tky)/tan(theta);
         kr_cur =  4258*Ts*sqrt(pkx[ai]*pkx[ai]+pky[ai]*pky[ai]+pkz[ai]*pkz[ai]);
         kr_next =  4258*Ts*sqrt(tkx*tkx+tky*tky+tkz*tkz);
         n_des = calcn(kr_next, sin(theta), FOV1, FOV2, NINT, kmax, dcf,nosmooth); 
	 n_act = anglebetween_new(pkx[ai],pky[ai],tkx,tky)/(kr_cur-kr_next);
	 if (n_act>=n_des) {
	   pkx[ai+1] = tkx;
	   pky[ai+1] = tky;
	   pkz[ai+1] = tkz;
	 }

	 if (n_act<n_des) {  /*if range_min has n_act<n_des, search for the crossover point */
      /* 3. Do bissection 10 or so times to find the crossover point */
		 for (bi = 0; bi<10; bi++) {
			 tt = range_min+(range_max-range_min)/2.0;
                         tkx = pkx[ai]+gxy*(cos(kang)*cos(tt)-sin(kang)*sin(tt));
                         tky = pky[ai]+gxy*(cos(kang)*sin(tt)+sin(kang)*cos(tt));
	                 tkz = sqrt(tkx*tkx+tky*tky)/tan(theta);
                         kr_cur =  4258*Ts*sqrt(pkx[ai]*pkx[ai]+pky[ai]*pky[ai]+pkz[ai]*pkz[ai]);
                         kr_next =  4258*Ts*sqrt(tkx*tkx+tky*tky+tkz*tkz);
                         n_des = calcn(kr_next, sin(theta), FOV1, FOV2, NINT, kmax, dcf,nosmooth); 
	                 n_act = anglebetween_new(pkx[ai],pky[ai],tkx,tky)/(kr_cur-kr_next);
			 if (n_act>=n_des) {
	                    pkx[ai+1] = tkx;
	                    pky[ai+1] = tky;
	                    pkz[ai+1] = tkz;
                            range_max = tt;
			 } else {
			    range_min = tt;
			 }
		 }
	 } else { /*range_min has n_act>=n_des*/
	    if ((numsolsmin>0) && (max(0,tmin)>max(0,t1))) {
		 /* Check at t1xy */
		 tkx = (pkx[ai]*cos(t1xy)+pky[ai]*sin(t1xy))/kxy*kxymin;
		 tky = (-pkx[ai]*sin(t1xy)+pky[ai]*cos(t1xy))/kxy*kxymin;
	         tkz = sqrt(tkx*tkx+tky*tky)/tan(theta);
                 kr_cur =  4258*Ts*sqrt(pkx[ai]*pkx[ai]+pky[ai]*pky[ai]+pkz[ai]*pkz[ai]);
                 kr_next =  4258*Ts*sqrt(tkx*tkx+tky*tky+tkz*tkz);
                 n_des = calcn(kr_next, sin(theta), FOV1, FOV2, NINT, kmax, dcf,nosmooth); 
	         n_act = anglebetween_new(pkx[ai],pky[ai],tkx,tky)/(kr_cur-kr_next);
		 if (n_act>=n_des) {
	           pkx[ai+1] = tkx;
	           pky[ai+1] = tky;
	           pkz[ai+1] = tkz;
		 }
		 if (n_act<n_des) { /* find crossover point */
			 for (bi=0; bi<10; bi++) {
                            tt = t1xy+(t2xy-t1xy)/2.0; 
		            tkx = (pkx[ai]*cos(tt)+pky[ai]*sin(tt))/kxy*kxymin;
		            tky = (-pkx[ai]*sin(tt)+pky[ai]*cos(tt))/kxy*kxymin;
	                    tkz = sqrt(tkx*tkx+tky*tky)/tan(theta);
                            kr_cur =  4258*Ts*sqrt(pkx[ai]*pkx[ai]+pky[ai]*pky[ai]+pkz[ai]*pkz[ai]);
                            kr_next =  4258*Ts*sqrt(tkx*tkx+tky*tky+tkz*tkz);
                            n_des = calcn(kr_next, sin(theta), FOV1, FOV2, NINT, kmax, dcf,nosmooth); 
	                    n_act = anglebetween_new(pkx[ai],pky[ai],tkx,tky)/(kr_cur-kr_next);
			    if (n_act>=n_des) {
	                       pkx[ai+1] = tkx;
	                       pky[ai+1] = tky;
	                       pkz[ai+1] = tkz;
                               t2xy = tt;
			    } else {
			       t1xy = tt;
			    }
			 }
		 }


	    }
	 }
             pgz[ai+1] = pkz[ai+1]-pkz[ai];
             pgy[ai+1] = pky[ai+1]-pky[ai];
             pgx[ai+1] = pkx[ai+1]-pkx[ai];

	 } else {
		 fail = 1; /* Due to there being no points in the range which have n_act>n_des. */
	 }
      } else {
	      fail = 1;  /* Due to no points existing */
      }
      
      if (n_act<1e-8) {
               done = 1;
      }

      kxy = sqrt(pkx[ai+1]*pkx[ai+1]+pky[ai+1]*pky[ai+1]);
      gxy = sqrt(pgx[ai+1]*pgx[ai+1]+pgy[ai+1]*pgy[ai+1]);
      Smax_c = min(Smax_z/cos(theta),Smax_xy/sin(theta));
      Gmax_c = min(Gmax_z/cos(theta),Gmax_xy/sin(theta));
      ttot = rto(sqrt(gxy*gxy+pgz[ai+1]*pgz[ai+1]),sqrt(kxy*kxy+pkz[ai+1]*pkz[ai+1])*4258*Ts,Smax_c,Gmax_c,4258,Ts,0,0);

      if (ttot==-1) {
	      fail = 1;
      }

      ai = ai+1;

  }
  if ((done == 1) && (fail == 0))
    return ai;
  else 
  {
    return -1;
  }

}
      
       
double fd2k(double k, double dk, double n, double dn, double SMS, double cmpsz) {
  double a,b,c,d,descrim,out;
  a = -dk*dk*n*n*k;
  b = 1;
  c = dk*dk*(2*n+k*dn);
  d = k*n;
  descrim = 2*d*c*b*a-b*b*c*c+b*b*SMS-d*d*a*a+d*d*SMS;
  descrim = (descrim>0)?descrim:0;
  out = -((d*c+b*a)-sqrt(descrim))/(b*b+d*d);
  out = min(out,cmpsz);
  out = max(out,-cmpsz); 
  return out;
}

double apod(double k, double kmax,double falloff, double alpha) {
	double tt;
 tt = (0.5+0.5*cos(PI*max(0,1/falloff*(k/kmax-1)+1)))*exp(-(alpha*alpha*k*k/(kmax*kmax))); 
 return max(1e-10,tt);
}

int whirlcone(  double theta,
                double FOV1,
                double FOV2,
		double FOV3,
		double FOV4,
                double dcf,
                double kmax,
                double NINT,
                int numpt,
                double Ts,
                double Smax_xy,
                double Smax_z,
                double Gmax_xy,
                double Gmax_z,
		double dens_min,
                double *pg,
                double *pk,
                double *pai,
		double *pn,
		double falloff,
		double alpha,
		int numptx) {

double *pkx,*pky,*pkz,*pgx,*pgy,*pgz;
double *pkfx,*pkfy,*pkfz,*pgfx,*pgfy,*pgfz;
double *prkx,*prky,*prkz,*prgx,*prgy,*prgz;
double *pkr;
double *pcsndk;
double kstart;
double n;
double dn;
double dk;
double d2k;
double cmp;
double ggmax;
double a,b,c,d;
double descrim;
double SMS;
double gxyend;
double kxyend;
double gxyzend;
double kxyzend;
double mg;
double k;
double maxs;
double maxg;
double Smax_c;
double Gmax_c;
double gxy_in;
double gxy_old;
double sintheta;
double costheta;
   double n_kmax,kmaxg,Gc,n_cur,phi_high,phi_low,dphi,gzt;
   double Gc_top, Gc_bot;
   double gxt,gyt,kxt,kyt,kzt;
   double rgxt,rgyt,rgzt;
int rtzlen;
int rtolen;
int dlen;
int ai,bi,ci;
int rangelow,rangehigh;
int maxrangehigh;
double dphimax;
double q1,q2,q3,q4,r1,r2,r3,r4,qq,rr;
double rtzlenhigh, rtzlenlow;
double c1, c2;
double dorungekutta;
double nosmooth;
nosmooth = 0; /*tohere*/
dorungekutta = 0;
  
pkr = mxCalloc((numptx),sizeof(double));
pcsndk = mxCalloc((numptx),sizeof(double));
prgx = mxCalloc((numptx),sizeof(double));
prgy = mxCalloc((numptx),sizeof(double));
prgz = mxCalloc((numptx),sizeof(double));
prkx = mxCalloc((numptx),sizeof(double));
prky = mxCalloc((numptx),sizeof(double));
prkz = mxCalloc((numptx),sizeof(double));
pgx = mxCalloc((numptx),sizeof(double));
pgy = mxCalloc((numptx),sizeof(double));
pgz = mxCalloc((numptx),sizeof(double));
pkx = mxCalloc((numptx),sizeof(double));
pky = mxCalloc((numptx),sizeof(double));
pkz = mxCalloc((numptx),sizeof(double));


pkfx = pk; pkfy = pk+numptx; pkfz = pk+2*numptx;
pgfx = pg; pgfy = pg+numptx; pgfz = pg+2*numptx;

  sintheta = sin(theta);
  costheta = cos(theta);


/* Find the point that the twist begins */
  if (dcf<0) {
    kstart = 0;
  } else if (dcf==0.5)  {
     c1 = 4*PI*PI*sintheta*sintheta*FOV2*FOV2/(NINT*NINT);
     c2 = (FOV2*FOV2)/(FOV1*FOV1);
     kstart = sqrt(c1*c2)/c1;
     if (kstart<0.5*kmax)
      kstart = 2*sqrt(c1*c2)*kmax*kmax/(c1*kmax*kmax+4*c2); 
  } else if (dcf>0) {
  kstart = sqrt(NINT*kmax/(2*PI*sin(theta)*FOV1));
} else {
  kstart = NINT/(2*PI*sin(theta)*FOV1);
}
SMS = Smax_xy*Smax_xy*4258*4258/(sin(theta)*sin(theta)*Ts*Ts);
	
if (kstart>(kmax-0.0001)) {
    mg = 0;
    ai = 0;
    pgfx[ai] = 0;
    pgfy[ai] = 0;
    pgfz[ai] = 0;
    pkfx[ai] = 0;
    pkfy[ai] = 0;
    pkfz[ai] = 0;
    if (theta<PI/4.0) {
	maxs = Smax_z*sqrt(1+tan(theta)*tan(theta));
	maxg = Gmax_z*sqrt(1+tan(theta)*tan(theta));
    } else {
	maxs = Smax_xy*sqrt(1+1/(tan(theta)*tan(theta)));
	maxg = Gmax_xy*sqrt(1+1/(tan(theta)*tan(theta)));
    }
    while ((k<kmax) && (ai<(numpt-3))) {
	    mg = min(maxg,mg+maxs);
	    k = k+mg*4258*Ts;
	    ai++;
	    pgfx[ai] = mg*sin(theta);
	    pgfy[ai] = 0;
	    pgfz[ai] = mg*cos(theta);
	    pkfx[ai] = pkfx[ai-1]+pgfx[ai]*4258*Ts;
	    pkfy[ai] = pkfy[ai-1]+pgfy[ai]*4258*Ts;
	    pkfz[ai] = pkfz[ai-1]+pgfz[ai]*4258*Ts;
    }
    *pai = ai;

} else {
	if (dcf<0) {
   pkr[0] = kstart;
} else {
   pkr[0] = kstart+0.001;
}
   pkx[0] = pkr[0]*sin(theta);
   pky[0] = 0;
   pkz[0] = pkr[0]*cos(theta);
   pgx[0] = 0;
   pgy[0] = 0;
   pgz[0] = 0;
   pcsndk[0] = 0;
   dk = 0;
   ai = 0;

   while ((ai<(numpt-1)) && (pkr[ai]<kmax)) {
if (0) {
      n = calcn(pkr[ai],sintheta,FOV1,FOV2,NINT,kmax,dcf,0);
      n = (n>0.00001)?n:(0.00001);
      dn = calcdn(pkr[ai],sintheta,FOV1,FOV2,NINT,kmax,dcf,0);
      /* Solve differential equation for Smax in x-y */
      a = -dk*dk*n*n*pkr[ai];
      b = 1;
      c = dk*dk*(2*n+pkr[ai]*dn);
      d = pkr[ai]*n;
      descrim = (2*d*c*b*a-b*b*c*c+b*b*SMS-d*d*a*a+d*d*SMS);
      descrim = (descrim>0)?descrim:0;
      d2k = -((d*c+b*a)-sqrt(descrim))/(b*b+d*d);
      /* Apply Smax constraint in z */
      cmp = Smax_z*4258/(costheta*Ts);
      d2k = (d2k<cmp)?d2k:cmp;
      d2k = (d2k>(-cmp))?d2k:(-cmp);
      dk = dk+d2k*Ts;
}

cmp = Smax_z*4258/(costheta*Ts);
q1 = Ts*dk;
r1 = Ts*fd2k(pkr[ai],dk,calcn(pkr[ai],sintheta,FOV1,FOV2,NINT,kmax,dcf,nosmooth),
		        calcdn(pkr[ai],sintheta,FOV1,FOV2,NINT,kmax,dcf,nosmooth),SMS,cmp);

if (dorungekutta) {

q2 = Ts*(dk+r1/2.0);
r2 = Ts*fd2k(pkr[ai]+q1/2.0,dk+r1/2.0,calcn(pkr[ai]+q1/2.0,sintheta,FOV1,FOV2,NINT,kmax,dcf,0),
		        calcdn(pkr[ai]+q1/2.0,sintheta,FOV1,FOV2,NINT,kmax,dcf,0),SMS,cmp);

q3 = Ts*(dk+r2/2.0);
r3 = Ts*fd2k(pkr[ai]+q2/2.0,dk+r2/2.0,calcn(pkr[ai]+q2/2.0,sintheta,FOV1,FOV2,NINT,kmax,dcf,0),
		        calcdn(pkr[ai]+q2/2.0,sintheta,FOV1,FOV2,NINT,kmax,dcf,0),SMS,cmp);

q4 = Ts*(dk+r3);
r4 = Ts*fd2k(pkr[ai]+q3,dk+r3,calcn(pkr[ai]+q3,sintheta,FOV1,FOV2,NINT,kmax,dcf,0),
		        calcdn(pkr[ai]+q3,sintheta,FOV1,FOV2,NINT,kmax,dcf,0),SMS,cmp);

qq = (q1+2*q2+2*q3+q4)/6;
rr = (r1+2*r2+2*r3+r4)/6;
} else {
qq = q1;
rr = r1;
}

 d2k = rr/Ts;
 dk = dk + rr;

      n = calcn(pkr[ai],sintheta,FOV1,FOV2,NINT,kmax,dcf,nosmooth);
      n = (n>0.00001)?n:(0.00001);
      dn = calcdn(pkr[ai],sintheta,FOV1,FOV2,NINT,kmax,dcf,nosmooth);
      /*mexPrintf("ai(%d) curk(%g) n(%g) dn(%g)\n",ai,pkr[ai],n,dn);*/

      /* Apply Gmax constraint in x-y */
      cmp = Gmax_xy*4258/(sqrt(1+pkr[ai]*pkr[ai]*n*n)*sintheta);
      dk = (dk<cmp)?dk:cmp;
      qq = min(qq,cmp/Ts);
      /* Apply Gmax constraint in z */
      cmp = Gmax_z*4258/costheta;
      dk = (dk<cmp)?dk:cmp;
      qq = min(qq,cmp/Ts);
      /* Apply constant density constraint */
      if (dens_min>0) {
      ggmax = NINT*sqrt((n*n*pkr[ai]*pkr[ai]*sintheta*sintheta)+1)/(2*PI*pkr[ai]*pkr[ai]*sintheta*dens_min*apod(pkr[ai],kmax,falloff,alpha)); 
      cmp = ggmax/sqrt(1+pkr[ai]*pkr[ai]*n*n*sintheta*sintheta)*4258;
      dk = (dk<cmp)?dk:cmp;
      qq = min(qq,cmp/Ts);
      }

      pkr[ai+1] = pkr[ai]+qq;/* dk*Ts;*/
/*      pcsndk[ai+1] = pcsndk[ai]+dk*n;*/
      pcsndk[ai+1] = pcsndk[ai]+dk*n;
      ai++;
   }
   dlen = ai;
   for (ai = 1; ai<dlen; ai++) {
	   pkx[ai] = pkr[ai]*cos(pcsndk[ai]*Ts)*sintheta;
	   pky[ai] = pkr[ai]*sin(pcsndk[ai]*Ts)*sintheta;
	   pkz[ai] = pkr[ai]*costheta;
	   pgx[ai] = (pkx[ai]-pkx[ai-1])/(Ts*4258);
	   pgy[ai] = (pky[ai]-pky[ai-1])/(Ts*4258);
	   pgz[ai] = (pkz[ai]-pkz[ai-1])/(Ts*4258);
   }
   
   rangelow = 1;
   rangehigh = 1;
   ai = rangehigh;
   rtzlen = 0;
   gxy_old = 0;
   rtzlenhigh = 100000;
   rtzlenlow = 100000;

   ai = 1;
   gxy_in =  sqrt(pgx[ai+1]*pgx[ai+1]+pgy[ai+1]*pgy[ai+1]);
   while ((gxy_in>(gxy_old+Smax_xy/1000)) && (ai<dlen-1)) {
	gxy_old = gxy_in;
	ai++;
	gxy_in =  sqrt(pgx[ai+1]*pgx[ai+1]+pgy[ai+1]*pgy[ai+1]);
   }
   maxrangehigh = ai-1;

   while ((rtzlen!=-1) && (rangehigh<maxrangehigh)) {
        rangehigh = min(maxrangehigh,rangehigh*2);
	ai = rangehigh;
	gxy_in =  sqrt(pgx[ai+1]*pgx[ai+1]+pgy[ai+1]*pgy[ai+1]);
        rtzlen =  rtz(-pgx[ai+1],-pgy[ai+1],-pgz[ai+1],pkx[ai],pky[ai],pkz[ai],theta,numpt,Gmax_xy,Gmax_z,Smax_xy,Smax_z,Ts,prgx,prgy,prgz,prkx,prky,prkz,FOV1,FOV2,NINT,kmax,dcf);
   }
   if ((rangehigh == dlen-2) && (rtzlen>0)) {
           n_kmax = calcn(kmax,sin(theta),FOV1,FOV2,NINT,kmax,dcf,0);
	   kmaxg = kmax/(4258*Ts)*sin(theta);
	   /* We reached the end of the pre-canned waveform, which means that we */
	   /* should search for the right Gc value to start with at kmax */

           /* Search over all Gc */
	   Gc_top = Gmax_xy;
	   Gc_bot = 0;
	   for (ci = 0; ci<10; ci++) {
		   Gc = Gc_bot+(Gc_top-Gc_bot)/2;

	   /* Start by finding an appropriate starting seed */
	       circcirc(kmaxg,0,kmaxg,Gc,1,&kxt,&kyt);
	       dphimax = atan2(kyt,kxt-kmaxg);
		   gxt = Gc*cos(dphimax);
		   gyt = Gc*sin(dphimax);
		   kxt = kmaxg+gxt; 
		   kyt = gyt; 
		   kzt = sqrt(kxt*kxt+kyt*kyt)/tan(theta);
	           gzt = kzt-kmaxg/tan(theta);
	           n_cur = atan2(kyt,kxt) / (4258*Ts*max(1e-8,(sqrt(kxt*kxt+kyt*kyt+kzt*kzt)-kmaxg/sin(theta))));
	   if (n_cur>n_kmax) {
	     phi_high = dphimax;
	     phi_low = 0;
	     for (bi = 0; bi < 20; bi++) {
		   dphi = phi_low+(phi_high-phi_low)/2;
		   gxt = Gc*cos(dphi);
		   gyt = Gc*sin(dphi);
		   kxt = kmaxg+gxt; 
		   kyt = gyt; 
		   kzt = sqrt(kxt*kxt+kyt*kyt)/tan(theta);
	           gzt = kzt-kmaxg/tan(theta);
	           n_cur = atan2(kyt,kxt) / (4258*Ts*max(1e-8,(sqrt(kxt*kxt+kyt*kyt+kzt*kzt)-kmaxg/sin(theta))));
		   if (n_cur>n_kmax) {
                     phi_high = dphi; 
		   } else {
                     phi_low = dphi;
		   }
	     }
	     dphi = phi_high; 
		   gxt = Gc*cos(dphi);
		   gyt = Gc*sin(dphi);
		   kxt = kmaxg+gxt; 
		   kyt = gyt; 
		   kzt = sqrt(kxt*kxt+kyt*kyt)/tan(theta);
	           gzt = kzt-kmaxg/tan(theta);

             rtzlen =  rtz(-gxt,-gyt,-gzt,kmaxg*4258*Ts,0*4258*Ts,kmaxg/tan(theta)*4258*Ts,theta,numpt,Gmax_xy,Gmax_z,Smax_xy,Smax_z,Ts,prgx,prgy,prgz,prkx,prky,prkz,FOV1,FOV2,NINT,kmax,dcf);
	   } else {
             rtzlen = -1;
	   }

	   if (rtzlen>0) {
               Gc_bot = Gc;
	       rgxt = -gxt;
	       rgyt = -gyt;
	       rgzt = -gzt;
	   } else {
	       Gc_top = Gc;
	   }

	}
	   if (Gc_bot>0) {
        rtzlen =  rtz(rgxt,rgyt,rgzt,kmaxg*4258*Ts,0*4258*Ts,kmaxg/tan(theta)*4258*Ts,theta,numpt,Gmax_xy,Gmax_z,Smax_xy,Smax_z,Ts,prgx,prgy,prgz,prkx,prky,prkz,FOV1,FOV2,NINT,kmax,dcf);
	   } else {
            rtzlen = -1;
	   }
	   ai++;

   } else {
   while ((rangehigh-rangelow)>0) {
	ai = ceil((rangehigh-rangelow)/2.0)+rangelow;
        rtzlen =  rtz(-pgx[ai+1],-pgy[ai+1],-pgz[ai+1],pkx[ai],pky[ai],pkz[ai],theta,numpt,Gmax_xy,Gmax_z,Smax_xy,Smax_z,Ts,prgx,prgy,prgz,prkx,prky,prkz,FOV1,FOV2,NINT,kmax,dcf);
	if ((rtzlen==-1) /*|| (rtzlen>rtzlenlow)*/) {
		rangehigh = ai-1;
		
	} else {
		rangelow = ai;
		rtzlenlow = rtzlen;
	}
   }
    ai = rangelow;

   rtzlen =  rtz(-pgx[ai+1],-pgy[ai+1],-pgz[ai+1],pkx[ai],pky[ai],pkz[ai],theta,numpt,Gmax_xy,Gmax_z,Smax_xy,Smax_z,Ts,prgx,prgy,prgz,prkx,prky,prkz,FOV1,FOV2,NINT,kmax,dcf);
   }

   /* COPY OVER SPOKE TRAJECTORY */
   if ((rtzlen>0)) {
      Smax_c = min(Smax_z/cos(theta),Smax_xy/sin(theta));
      Gmax_c = min(Gmax_z/cos(theta),Gmax_xy/sin(theta));

   gxyend = sqrt(prgx[rtzlen]*prgx[rtzlen]+prgy[rtzlen]*prgy[rtzlen]);
   kxyend = sqrt(prkx[rtzlen]*prkx[rtzlen]+prky[rtzlen]*prky[rtzlen]);
   gxyzend = sqrt(prgx[rtzlen]*prgx[rtzlen]+prgy[rtzlen]*prgy[rtzlen]+prgz[rtzlen]*prgz[rtzlen]);
   kxyzend = sqrt(prkx[rtzlen]*prkx[rtzlen]+prky[rtzlen]*prky[rtzlen]+prkz[rtzlen]*prkz[rtzlen]);
   rtolen = rto(gxyzend,kxyzend*4258*Ts,Smax_c,Gmax_c,4258,Ts,pgfx,pkfx);
        for (bi = 0; bi<rtolen; bi++) {
            pgfz[bi] = pgfx[bi]*cos(theta);
            pgfy[bi] = pgfx[bi]*sin(theta)*prky[rtzlen]/kxyend;
            pgfx[bi] = pgfx[bi]*sin(theta)*prkx[rtzlen]/kxyend;
            pkfx[bi] = ((bi>0)?pkfx[bi-1]:0) + pgfx[bi]*4258*Ts;
            pkfy[bi] = ((bi>0)?pkfy[bi-1]:0) + pgfy[bi]*4258*Ts;
            pkfz[bi] = ((bi>0)?pkfz[bi-1]:0) + pgfz[bi]*4258*Ts;
	}

   /* COPY OVER CIRCULAR TRAJECTORY */
   for (bi = 0; bi<min(rtzlen,numpt-rtolen); bi++) {
       pgfx[bi+rtolen] = -prgx[rtzlen-bi];
       pgfy[bi+rtolen] = -prgy[rtzlen-bi];
       pgfz[bi+rtolen] = -prgz[rtzlen-bi];
       pkfx[bi+rtolen] = prkx[rtzlen-bi-1]*4258*Ts;
       pkfy[bi+rtolen] = prky[rtzlen-bi-1]*4258*Ts;
       pkfz[bi+rtolen] = prkz[rtzlen-bi-1]*4258*Ts;
       
   }
   /* COPY OVER DIFFERENTIAL TRAJECTORY */
   for (bi = ai+1; bi<min(dlen,numpt-rtolen-rtzlen+(ai+1)); bi++) {
       pgfx[bi-(ai+1)+rtolen+rtzlen] = pgx[bi];
       pgfy[bi-(ai+1)+rtolen+rtzlen] = pgy[bi];
       pgfz[bi-(ai+1)+rtolen+rtzlen] = pgz[bi];
       pkfx[bi-(ai+1)+rtolen+rtzlen] = pkx[bi];
       pkfy[bi-(ai+1)+rtolen+rtzlen] = pky[bi];
       pkfz[bi-(ai+1)+rtolen+rtzlen] = pkz[bi];
   }
   
    
*pai = min(numpt,dlen-ai+rtolen+rtzlen)-1;
   } else {

   Smax_c = min(Smax_z/cos(theta),Smax_xy/sin(theta));
   Gmax_c = min(Gmax_z/cos(theta),Gmax_xy/sin(theta));
   kxyend = sqrt(pkx[0]*pkx[0]+pky[0]*pky[0]);
   kxyzend = sqrt(pkx[0]*pkx[0]+pky[0]*pky[0]+pkz[0]*pkz[0]);
   rtolen = rto(0,kxyzend,Smax_c,Gmax_c,4258,Ts,pgfx,pkfx);
        for (bi = 0; bi<min(rtolen,numpt); bi++) {
            pgfz[bi] = pgfx[bi]*cos(theta);
            pgfy[bi] = pgfx[bi]*sin(theta)*pky[0]/kxyend;
            pgfx[bi] = pgfx[bi]*sin(theta)*pkx[0]/kxyend;
            pkfx[bi] = ((bi>0)?pkfx[bi-1]:0) + pgfx[bi]*4258*Ts;
            pkfy[bi] = ((bi>0)?pkfy[bi-1]:0) + pgfy[bi]*4258*Ts;
            pkfz[bi] = ((bi>0)?pkfz[bi-1]:0) + pgfz[bi]*4258*Ts;
	}

   for (bi = 0; bi<min(dlen,numpt-rtolen-1); bi++) {
       pgfx[rtolen+bi] = pgx[bi];
       pgfy[rtolen+bi] = pgy[bi];
       pgfz[rtolen+bi] = pgz[bi];
       pkfx[rtolen+bi] = pkx[bi];
       pkfy[rtolen+bi] = pky[bi];
       pkfz[rtolen+bi] = pkz[bi];
   }
   *pai = min(dlen,numpt-rtolen-1)+min(rtolen,numpt);
   }

}

return 1;
          
}


/* mexFunction is the gateway routine for the MEX-file. */ 
void
mexFunction( int nlhs, mxArray *plhs[],
             int nrhs, const mxArray *prhs[] )
{
  int mrows;
  int ncols;
  double * pr;
  double theta;
  double FOV1;
  double FOV2;
  double FOV3;
  double FOV4;
  double dcf;
  double kmax;
  double NINT;
  double numpt;
  double Ts;
  double Smax_xy;
  double Smax_z;
  double Gmax_xy;
  double Gmax_z;
  double dens_min;
  double falloff,alpha;
  double *pg;
  double *pk;
  double *pai;
  double *pn;
  int numptx;

  pr = mxGetPr(prhs[0]);
  theta = *pr;  
  if (theta<=0) {
	  theta = 0.00001;
  } else if (theta>= PI/2) {
	  theta = PI/2-0.00001;

  }
  

  pr = mxGetPr(prhs[1]);
  FOV1 = *pr++;  
  if (mxGetNumberOfElements(prhs[1]) < 2)
     FOV2 = FOV1;
  else
     FOV2 = *pr++;

  if (mxGetNumberOfElements(prhs[1]) < 4) {
     FOV3 = FOV1;
     FOV4 = FOV2;
  } else {
     FOV3 = *pr++;
     FOV4 = *pr;
  }

  pr = mxGetPr(prhs[2]);
  dcf = *pr;  

  pr = mxGetPr(prhs[3]);
  kmax = *pr;  

  pr = mxGetPr(prhs[4]);
  NINT = *pr;  

  if (nrhs<6)
     numpt = 2000;
  else {
     pr = mxGetPr(prhs[5]);
     numpt = *pr;  
     if (numpt>100000)
       numpt = 100000;
  }

  if (nrhs<7)
     Ts = 0.000002;
  else {
     pr = mxGetPr(prhs[6]);
     Ts = *pr;  
  }
  
  if (nrhs<8) {
     Smax_xy = Ts*15000;
     Smax_z = Ts*15000;
  } else {
     pr = mxGetPr(prhs[7]);
     if (mxGetNumberOfElements(prhs[7]) < 2)  {
       Smax_xy = *pr;  
       Smax_z = *pr;  
     } else {
       Smax_xy = *pr++;  
       Smax_z = *pr;  
     }
  }

  if (nrhs<9) {
     Gmax_xy = 3.89;
     Gmax_z = 3.89;
  } else {
      pr = mxGetPr(prhs[8]);
     if (mxGetNumberOfElements(prhs[8]) < 2)  {
       Gmax_xy = *pr;  
       Gmax_z = *pr;  
     } else {
       Gmax_xy = *pr++;  
       Gmax_z = *pr;  
     }
  }

  if (nrhs<10)
     dens_min = 0.0000001;
  else {
      pr = mxGetPr(prhs[9]);
      dens_min = *pr;  
  }

  if (nrhs<11) {
	  alpha = 0.0;
	  falloff = 1;
  } else {
	  pr = mxGetPr(prhs[10]);
	  alpha = *pr++;
          if (mxGetNumberOfElements(prhs[8]) < 2)  {
           falloff = 1;
	  } else {
           falloff = *pr;
	  }
  }

/*mexPrintf("alpha(%g), falloff(%g)\n",alpha,falloff);*/
   mrows = 4;
   ncols = 3;
 /* Look at each input (right-hand-side) argument. */
   numptx = numpt+500;

  plhs[0] = mxCreateDoubleMatrix( (int)(numptx), 3, mxREAL);
  plhs[1] = mxCreateDoubleMatrix( (int)(numptx), 3, mxREAL);
  plhs[2] = mxCreateDoubleMatrix( 1,1,mxREAL);
  
  pg = mxGetPr(plhs[0]);
  pk = mxGetPr(plhs[1]);
  pai = mxGetPr(plhs[2]);
  whirlcone(theta,FOV1,FOV2,FOV3,FOV4,dcf,kmax,NINT,numpt,Ts,Smax_xy,Smax_z,Gmax_xy,Gmax_z,dens_min,pg,pk,pai,pn,falloff,alpha,numptx);

}
