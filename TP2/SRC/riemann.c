#include "riemann.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

double g = 9.81;


/*
int main(void){

  double hL=2;
  double uL=0;

  double hR=1;
  double uR=0;

  double wL[2]={hL, hL * uL};		// variables conservatives
  double wR[2]={hR, hR * uR};

  double xi = 0;
  double w[2];
  
  riem_stvenant(wL, wR, xi, w);
  plot_riem(wL,wR);

  return 0;
}
*/

// INUTILE POUR NOTRE CAS QUI EST 1D
void flux_riem_2d(double *wL, double *wR, double *vnorm, double *flux){

  double qnL = wL[1] * vnorm[0] +  wL[2] * vnorm[1]; 
  double qnR = wR[1] * vnorm[0] +  wR[2] * vnorm[1];

  double qtL = -wL[1] * vnorm[1] +  wL[2] * vnorm[0];
  double qtR = -wR[1] * vnorm[1] +  wR[2] * vnorm[0];

  double vL[2] = {wL[0], qnL};
  double vR[2] = {wR[0], qnR};

  double v[2];
  double xi = 0;

  riem_stvenant(vL, vR, xi, v);

  double un = v[1] / v[0];

  double ut;

  if (un > 0)
    ut = qtL / wL[0];
  else
    ut = qtR / wR[0];

  double qn = v[1];
  double qt = ut * v[0];

  double w[3];
  w[0] = v[0];
  w[1] = qn * vnorm[0] - qt * vnorm[1];
  w[2] = qn * vnorm[1] + qt * vnorm[0];

}

void riem_stvenant(double *wL,
		   double *wR,
		   double xi,
		   double *w){

  double hL = wL[0];
  double uL = wL[1]/wL[0];
  double hR = wR[0];
  double uR = wR[1]/wR[0];

  double hs = 1e-6;
  int itermax = 10;
  
  // printf("NOUVEAU PB DE RIEMANN\n");
  
  for(int it = 0; it < itermax; it++){
    double f = uL - (hs - hL) * Z(hs, hL) -
      uR - (hs - hR) * Z(hs, hR);
    double df = -(hs - hL) * dZ(hs, hL) -
      Z(hs, hL) -
      (hs - hR) * dZ(hs, hR) -
      Z(hs, hR);
    double dhs = -f / df;

    hs += dhs;

    // printf("it=%d f=%e df=%e hs=%e dhs=%e\n",it,f,df,hs,dhs);
    
  }

  double us = uL - (hs - hL) * Z(hs, hL);

  double v1m, v1p, v2m, v2p;

  // 1-onde
  if (hs < hL){
    v1m = uL - sqrt(g * hL);
    v1p = us - sqrt(g * hs);
  } else {
    double a = sqrt(hs) / (sqrt(hs) + sqrt(hL));
    double u = a * us + (1 - a) * uL;
    double h = (hs + hL) / 2;
    v1m = u - sqrt(g * h);
    v1p = v1m;
  }

  // 2 onde
  if (hs < hR){
    v2m = us + sqrt(g * hs);
    v2p = uR + sqrt(g * hR);
  } else {
    double a = sqrt(hs) / (sqrt(hs) + sqrt(hR));
    double u = a * us + (1 - a) * uR;
    double h = (hs + hR) / 2;
    v2m = u + sqrt(g * h);
    v2p = v2m;
  }

  //printf("v=%f %f %f %f\n hs=%f us=%f\n", v1m,v1p,v2m,v2p, hs,us);
  
  if (xi < v1m) {
    w[0] = wL[0];
    w[1] = wL[1];
  } else if (xi < v1p){
    double u = (uL + 2 * xi + 2 *sqrt(g * hL)) / 3;
    double h = (u - xi) * (u - xi) / g;
    w[0] = h;
    w[1] = h * u;
  } else if (xi < v2m){
    w[0] = hs;
    w[1] = hs * us;
  } else if (xi < v2p){
    double u = (uR + 2 * xi - 2 *sqrt(g * hR)) / 3;
    double h = (u - xi) * (u - xi) / g;
    w[0] = h;
    w[1] = h * u;
  } else {
    w[0] = wR[0];
    w[1] = wR[1];
  }

}

void plot_riem(double *wL,
	       double *wR){

  FILE *plotfile;

  double xmin = -5;
  double xmax = 5;

  int n = 1000;
  double dx = (xmax - xmin)/n;
  plotfile = fopen("riem.dat", "w");
  for(int i = 0; i < n; i++){
    double xi = xmin+i*dx;
    double w[2];
    riem_stvenant(wL,wR,xi,w);
    double h = w[0];
    double u = w[1]/w[0];
    fprintf(plotfile, "%f %f %f\n",
	    xi, h, u);
  }
  fclose(plotfile);

  system("gnuplot riemcom");

}

double Heaviside(double x){
  if (x > 0)
    return 1;
  else
    return 0;
}

double Dirac(double x){
    return 0;
}


double Z(double hs,
	 double h){

double t0 = 2.0*sqrt(g)/(sqrt(hs)+sqrt(h))*Heaviside(h-hs)+sqrt(2.0)*sqrt(g)*sqrt(h+hs)/sqrt(h*hs)/2.0-sqrt(2.0)*sqrt(g)*sqrt(h+hs)/sqrt(h*hs)*Heaviside(h-hs)/2.0;

 return t0;
  
}

double dZ(double hs,
	  double h){

   double t0 = -sqrt(g)/pow(sqrt(hs)+sqrt(h),2.0)*Heaviside(h-hs)/sqrt(hs)-2.0*sqrt
(g)/(sqrt(hs)+sqrt(h))*Dirac(-h+hs)+sqrt(2.0)*sqrt(g)/sqrt(h+hs)/sqrt(h*hs)/4.0
-sqrt(2.0)*sqrt(g)*sqrt(h+hs)/sqrt(h*h*h*hs*hs*hs)*h/4.0-sqrt(2.0)*sqrt(g)/sqrt
(h+hs)/sqrt(h*hs)*Heaviside(h-hs)/4.0+sqrt(2.0)*sqrt(g)*sqrt(h+hs)/sqrt(h*h*h*
hs*hs*hs)*Heaviside(h-hs)*h/4.0+sqrt(2.0)*sqrt(g)*sqrt(h+hs)/sqrt(h*hs)*Dirac(-
h+hs)/2.0;

   return t0;

}
