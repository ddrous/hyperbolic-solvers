#ifndef _RIEM_H
#define _RIEM_H


void riem_stvenant(double *wL,
		   double *wR,
		   double xi,
		   double *w);

double Z(double h1,
	 double h2);

double dZ(double h1,
	  double h2);

void plot_riem(double *wL,
	       double *wR);

void flux_riem(double *wL, double *wR,double *flux);

#endif
