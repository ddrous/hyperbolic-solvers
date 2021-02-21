// QUELQUES COMMANDES DE COMPILATTION:
// clang godunov.c -lm -g -fsanitize=address
// ou
// clang godunov.c -lm -O

#include <assert.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "godunov.h"
#include "riemann.h"

// Gravité
#define _G 9.81

int main(void){

    bool etudeCvg = false;     

    if (!etudeCvg){              // On résout un seul problème et on l'affiche
        godunov gd = {0};

        /* CAS TEST */
        // double xmin = -1;       // Ne pas commencer en 0 car le problleme de Rieman en 0
        // double xmax = 1;
        // double cfl = 0.5;
        // int m = 2;
        // int N = 50;
        // double tmax = 0.1;

        double xmin = -10;       // Ne pas commencer en 0 car le problleme de Rieman en 0
        double xmax = 10;
        double cfl = 0.5;
        int m = 2;
        int N = 1000;
        double tmax = 1;      // mettre 1.366671, puis 1.808131, puis 2.048256

        // 
        // DIFFCULTË APPARIT EN (0,1) avec si on prend une onde sonique (voir page 3)

        // Essayer aussi MUCL, ca devrait etre meilleur 

        godunov_init(&gd, xmin, xmax, cfl, m, N);
        godunov_solve(&gd, tmax);

        double err[2];
        godunov_error(&gd, err);
        // printf("N: %d  \ndx: %f \nErreur sur h: %f  \nErreur sur u: %f \n", gd.N, gd.dx, err[0], err[1]);
        printf("N: %d  \ndx: %f \nErreur : %f \n", gd.N, gd.dx, err[0]+err[1]);

        godunov_plot(&gd, true);
        godunov_free(&gd);
    } else                  // On résout plusieur problemes
    {  
        // Calcul des erreurs L1 pour plusieurs problemes et plot
        int NMin = 100;
        int NMax = 500;
        int nbSimuMax = 8;          // Nombre maximal de simu a effectuer

        printEstimatedOrder(NMax/100);
        plotError(NMin, NMax, nbSimuMax);
    }
    
    return 0;
}


/*******************************************************************************
 * POUR L'EQUATION DE ST VENANT
 ******************************************************************************/

void solexacte(double x, double t, double* w) {

    /* CAS TEST DE BASE */
    // double hL=2;
    // double uL=0;

    // double hR=1;
    // double uR=0;

    /* VFROE AVEC POINT SONIQUE */ 
    /************************************** */
    // double dv = 1 - 1.56; // on décale la vitesse a droite pour enlever le point sonique. Quel est le plus decallalge qu'on peut faire pour faire sortie le point sonique de 0.
    // double dv = 1.0; // 
    double dv = 0.0; // 
    double hL=1;
    double uL=-1 + dv;

    double hR=0.25;
    double uR= uL + 2*sqrt(_G)*(sqrt(hL)-sqrt(hR));
    /************************************** */

    // CORRECTION ENTROPIQUE
    // double lambda_1L = uL - sqrt(_G*hL);
    // double lambda_2L = uL + sqrt(_G*hL);
    // double lambda_1R = uR - sqrt(_G*hR);
    // double lambda_2R = uR + sqrt(_G*hR);
    // printf("%f, %f, %f, %f \n", lambda_1L, lambda_1R, lambda_2L, lambda_2R);


    double wL[2]={hL, hL * uL};		    // variables conservatives
    double wR[2]={hR, hR * uR};

    double z = x/ (t + 1e-12);           // A cause des initialisations

    //   riem_stvenant(wL, wR, z, w);
    riem_stvenant(wL, wR, z, w);
}

void fluxphy(double* w, double* flux) { 
    flux[0] = w[1]; 
    flux[1] = (w[1]*w[1]/w[0]) + (_G*w[0]*w[0]/2); 
}

double lambda_max(double* w) {
    double h = w[0]; 
    double u = w[1] / w[0]; 
    return fabs(u + sqrt(_G*h));
}

/* VFROE */ // apres on peu faire MUSCL + VFROE pou eassayer
// Plus precis que Rusanov, bcp de temps de calcul par rapport a volume finis
void riem_vfroe(double *wL, double *wR, double z, double *w){
    double hL = wL[0];
    double uL = wL[1] / wL[0];

    double hR = wR[0];
    double uR = wR[1] / wR[0];

    double hBar = (hL + hR) /2.0;
    double uBar = (uL + uR) /2.0;
    double cBar = sqrt(_G * hBar);

    double lambda_1 = uBar - cBar;
    double lambda_2 = uBar + cBar;
    
    // CORRECTION ENTROPIQUE
    // double lambda_1L = uL - sqrt(_G*hL);
    // double lambda_2L = uL + sqrt(_G*hL);
    // double lambda_1R = uR - sqrt(_G*hR);
    // double lambda_2R = uR + sqrt(_G*hR);

    double h, u;
    if ((lambda_1 > 0) && (lambda_2 > 0)){
        h = hL;
        u = uL;
    } else if (lambda_1 < 0 && lambda_2 < 0){
        h = hR;
        u = uR;
    // } else if (lambda_1 < 0 && lambda_2 > 0){
    } else {
        h = hBar - (hBar*(uR-uL))/(2*sqrt(_G*hBar));
        u = uBar - (_G*(hR-hL))/(2*sqrt(_G*hBar));
    }

    w[0] = h;
    w[1] = h*u;
}



/* VFROE */ // apres on peu faire MUSCL + VFROE pou eassayer
// Plus precis que Rusanov, bcp de temps de calcul par rapport a volume finis
void flux_rusanov(double *wL, double *wR, double *flux){
    double hL = wL[0];
    double uL = wL[1] / wL[0];

    double hR = wR[0];
    double uR = wR[1] / wR[0];

    double lambdaL = fabs(uL) + sqrt(_G*hL);
    double lambdaR = fabs(uR) + sqrt(_G*hR);
    
    double lambda = lambdaL>lambdaR ? lambdaL : lambdaR;

    double fL[2]; 
    double fR[2]; 
    fluxphy(wL, fL);
    fluxphy(wR, fR);

    flux[0] = (fL[0]+fR[0])/2. - lambda*(wR[0]-wL[0])/2.;
    flux[1] = (fL[1]+fR[1])/2. - lambda*(wR[1]-wL[1])/2.;
}


/*******************************************************************************
 * PARTIE COMMUNE A TOUS LES PROBLEMES
 ******************************************************************************/

void fluxnum(double* a, double* b, double* flux) {
    double w[2];       

    // //-------------------------------------------------------------
    double *wL = a;
    double *wR = b;    

    double hL = wL[0];
    double uL = wL[1] / wL[0];

    double hR = wR[0];
    double uR = wR[1] / wR[0];

    
    // CORRECTION ENTROPIQUE
    double lambda_1L = uL - sqrt(_G*hL);
    double lambda_2L = uL + sqrt(_G*hL);
    double lambda_1R = uR - sqrt(_G*hR);
    double lambda_2R = uR + sqrt(_G*hR);

    // /* */
    // // --------------------------DEBUT CORRECTION ENTROPIQUE FACILE
    double epsL = 0;
    double epsR = 0;
    if (lambda_1L < 0 && lambda_1R > 0){
        epsL = fmin(-lambda_1L, lambda_1R);     // plus pttit des deux decallalges
        // printf("lambda 1L, %f")
    }
    if (lambda_2L < 0 && lambda_2R > 0)
        epsR = fmin(-lambda_2L, lambda_2R);
    // double eps = (epsR>epsL)?epsR:epsL;
    double eps = fmax(epsL,epsR);
    //-------------------------------------------------------------

    // riem_stvenant(a, b, 0., w);
    riem_vfroe(a, b, 0., w);
    fluxphy(w, flux);


    // Arangeons le flux-----COmmenter ca pour enlever la correction entropique
    for (int i = 0; i < 2; i++)
    {
        flux[i] -= eps/2.0  * (wR[i]-wL[i]); /// rajouter ls viscosité
    }

    // // flux_rusanov(a, b, flux);
    // // -------------------------- FIN CORRECTION ENTROPIQUE FACILE


    //---------Question 8--------------------------

    // if ((lambda_1L<0 && lambda_1R>0) || (lambda_2L<0 && lambda_2R>0))
    // // if (fabs(lambda_1R-lambda_1L)<tol || fabs(lambda_2R-lambda_2L)<tol)
    //     flux_rusanov(a, b, flux);
    // else{
    //     // printf("I'M IN HERE\n");
    //     riem_vfroe(a, b, 0., w);
    //     fluxphy(w, flux);
    // }

    // riem_stvenant(a, b, 0., w);
    // fluxphy(w, flux);

    //-------------------- FIn question 8----------


}

void godunov_init(godunov* gd, double xmin, double xmax, double cfl, int m,
                  int N) {
    gd->xmin = xmin;
    gd->xmax = xmax;
    gd->m = m;
    gd->N = N;
    gd->cfl = cfl;
    gd->dx = (xmax - xmin) / N;
    gd->dt = 0;
    gd->tfin = 0;

    gd->xi = malloc((N + 2) * sizeof(double));
    gd->un = malloc((N + 2) * sizeof(double) * m);
    gd->unp1 = malloc((N + 2) * sizeof(double) * m);

    for (int i = 0; i < N + 2; i++) {
        gd->xi[i] = xmin + gd->dx / 2 + (i - 1) * gd->dx;
        double t = 0;
        solexacte(gd->xi[i], t, gd->un + i * m);        // La solution initiale est la solution exacte a t=0
    }
}

void godunov_free(godunov* gd) {
    free(gd->xi);
    free(gd->un);
    free(gd->unp1);
}



double minmod(double a, double b, double c)
{
    double res = 0.0;
    if ((a > 0) && (b > 0) && (c > 0))
        res = a < b ? (a < c ? a : c) : (b < c ? b : c);    // Pas sur
    else if ((a < 0) && (b < 0) && (c < 0))
        res = a > b ? (a > c ? a : c) : (b > c ? b : c);    // Pas sur
    return res;
}


void godunov_solve(godunov* gd, double tmax) {
    double tnow = 0;
    int m = gd->m;

    while (tnow < tmax) {
        double vmax = 0;
        // calcul de la vitesse max
        // #pragma omp parallel for
        for (int i = 0; i < gd->N + 2; i++) {
            double vloc = lambda_max(gd->un + m * i);
            vmax = vmax > vloc ? vmax : vloc;
        }

        // Calcul des pentes pour MUSCL----------------------------------
        double si[m*(gd->N+2)], ri[m*(gd->N+2)];
        for (int i = 1; i < gd->N + 1; i++) {
            for (int iv = 0; iv < m; iv++) {
                double alpha = (gd->un[i * m + iv] - gd->un[(i-1) * m + iv])/gd->dx;
                double beta = (gd->un[(i+1) * m + iv] - gd->un[i * m + iv])/gd->dx;
                double gamma = (gd->un[(i+1) * m + iv] - gd->un[(i-1) * m + iv])/(2.0*gd->dx);
                si[i * m + iv] = minmod(alpha, beta, gamma);
                // si[i * m + iv] = 0;
            }   

            // ri[i * m + 0] = -si[i * m +1];
            // ri[i * m + 1] = -(_G*gd->un[i*m] - gd->un[i*m+1]*gd->un[i*m+1]/(gd->un[i*m]*gd->un[i*m]))*si[i*m+0] - 2*si[i*m+1]*gd->un[i*m+1]/gd->un[i*m];
        }
        // ----------------------------------------------------------------


        gd->dt = gd->cfl * gd->dx / vmax;

        // #pragma omp parallel for
        for (int i = 1; i < gd->N + 1; i++) {            
            // fluxnum(gd->un + i * m, gd->un + (i + 1) * m, flux);

            double uLplus[2], uRplus[2], uLmoins[2], uRmoins[2];
            // Application de MUSCL--------------------------------------------
            for (int iv = 0; iv < m; iv++) {
                uLplus[iv] = gd->un[i*m+iv] + si[i*m+iv]*gd->dx/2.0;
                uRplus[iv] = gd->un[(i+1)*m+iv] - si[(i+1)*m+iv]*gd->dx/2.0;

                uLmoins[iv] = gd->un[(i-1)*m+iv] + si[(i-1)*m+iv]*gd->dx/2.0;
                uRmoins[iv] = gd->un[i*m+iv] - si[i*m+iv]*gd->dx/2.0;

            }
            // printf("tnow = %f w-L = %0.9f w-R = %0.9f\n", tnow, uLmoins[0], uRmoins[0]);
            double fluxPlus[m];
            double fluxMoins[m];
            fluxnum(uLplus, uRplus, fluxPlus);
            fluxnum(uLmoins, uRmoins, fluxMoins);
            // printf("tnow = %f flux+ = %0.9f flux- = %0.9f\n", tnow, fluxPlus[0], fluxMoins[0]);

            double k1[m]; 
            for (int iv = 0; iv < m; iv++) {
                k1[iv] = - gd->dt * (fluxPlus[iv] - fluxMoins[iv])/gd->dx;
            }
            // printf("k1 = %f h[0] = %f v[0] = %f\n", k1[0], gd->unp1[i*m], gd->unp1[3]);

            double uLplusPrime[m], uRplusPrime[m], uLmoinsPrime[m], uRmoinsPrime[m];
            for (int iv = 0; iv < m; iv++) {
                uLplusPrime[iv] = uLplus[iv] + k1[iv];
                uRplusPrime[iv] = uRplus[iv] + k1[iv];

                uLmoinsPrime[iv] = uLmoins[iv] + k1[iv];
                uRmoinsPrime[iv] = uRmoins[iv] + k1[iv];
            }
            double fluxPlusPrime[m];
            double fluxMoinsPrime[m];
            fluxnum(uLplusPrime, uRplusPrime, fluxPlusPrime);
            fluxnum(uLmoinsPrime, uRmoinsPrime, fluxMoinsPrime);
            double k2[m]; 
            for (int iv = 0; iv < m; iv++) {
                k2[iv] = - gd->dt * (fluxPlusPrime[iv] - fluxMoinsPrime[iv])/gd->dx;

            }

            // ----------------------------------------------------------------
            for (int iv = 0; iv < m; iv++) {
                gd->unp1[i * m + iv] =  gd->un[i*m+iv] + (k1[iv] + k2[iv]) / 2.0;
            }
        }
        // mise à jour
        tnow += gd->dt;
        // printf("tnow = %f vmax = %f tmax = %f\n", tnow, vmax, tmax);
        // printf("tnow = %f k1 = %f k2 = %f\n", tnow, k1[i * m], k2[i * m]);

        // //**********************************************************************
        // // CONDITIONS AUX LIMITES TRICHÉES
        // //**********************************************************************
        int i = 0;
        solexacte(gd->xi[i], tnow, gd->unp1 + i * m);
        i = gd->N + 1;
        solexacte(gd->xi[i], tnow, gd->unp1 + i * m);
        // //**********************************************************************

        //**********************************************************************
        // PISCINE BORNÉE
        //**********************************************************************
        // int i = 0;
        // gd->unp1[i*m] = gd->unp1[(i+1)*m];
        // gd->unp1[i*m+1] = -gd->unp1[(i+1)*m+1];
        // i = gd->N + 1;
        // gd->unp1[i*m] = gd->unp1[(i-1)*m];
        // gd->unp1[i*m+1] = -gd->unp1[(i-1)*m+1];
        //**********************************************************************

        // // **********************************************************************
        // // PISCINE INFINIE
        // // **********************************************************************
        // int i = 0;
        // gd->unp1[i*m] = gd->unp1[(i+1)*m];
        // gd->unp1[i*m+1] = gd->unp1[(i+1)*m+1];
        // i = gd->N + 1;
        // gd->unp1[i*m] = gd->unp1[(i-1)*m];
        // gd->unp1[i*m+1] = gd->unp1[(i-1)*m+1];
        // //**********************************************************************

        memcpy(gd->un, gd->unp1, (gd->N + 2) * m * sizeof(double));
    }
    gd->tfin = tnow;
}

void godunov_plot(godunov *gd, bool visu)
{
    FILE *fic = fopen("godu.dat", "w");

    for (int i = 1; i < gd->N + 1; i++)
    {
        double uex[gd->m];
        // double u0[gd->m]; // Solution initiale
        solexacte(gd->xi[i], gd->tfin, uex);
        // solexacte(gd->xi[i], 0, u0);

        fprintf(fic, "%f %f %f %f %f\n", gd->xi[i], uex[0], gd->un[i*gd->m], uex[1]/uex[0], gd->un[i*gd->m+1]/gd->un[i*gd->m]);
    }

    fclose(fic);

    if (visu)
    {
        int status = system("gnuplot plotcom");        // plot la hauteur
        assert(status == EXIT_SUCCESS);
    }
}

void godunov_error(godunov* gd, double *err) {
    double erreur1 = 0;
    double erreur2 = 0;

    for (int i = 1; i < gd->N + 1; i++) {
        double uex[gd->m];
        solexacte(gd->xi[i], gd->tfin, uex);
        int iv = 0;     // Erreur sur h
        erreur1 += gd->dx * fabs(uex[iv] - gd->un[i * gd->m + iv]);

        iv = 1;         // Erreur sur h
        erreur2 += gd->dx * fabs(uex[iv]/uex[iv-1] - gd->un[i*gd->m+iv]/gd->un[i*gd->m+iv-1]);
    }

    err[0] = erreur1;
    err[1] = erreur2;
}

void plotError(int NMin, int NMax, int nbSimu){
    FILE * fic = fopen("erreur.dat", "w");

    /* Remplissage automatique des tailles de maillages */
    int pbSizes[nbSimu];

    int step = (NMax - NMin) / nbSimu;      
    int ctr = 0;
    // for (int i = NMin; i < NMax; i+step)     // Pour un maillage uniforme
    for (int i = NMin; i < NMax && ctr < nbSimu; i*=2)
    {
        pbSizes[ctr] = i;
        ctr ++;
    }
    nbSimu = ctr;       // Vrai nombre de simulation a effectuer

    /* -------------------------------------------------- */

    // Remplissage à la main des tailles de maillages
    // int pbSizes[] = {10, 50, 100, 500, 1000, 5000, 10000, 20000, 20000};
    // --------------------------------------------------
    
    /* Resolution des differents problemes */
    godunov gd = {0};

    // Constantes d'un prbleme a l'autre
    double xmin = -10;
    double xmax = 10;
    double cfl = 0.5;
    int m = 2;
    double tmax = 1;

    for (int i = 0; i < nbSimu; i++) {      // On n'utilise que les nbSimu premieres valeurs
        int N = pbSizes[i];

        godunov_init(&gd, xmin, xmax, cfl, m, N);
        godunov_solve(&gd, tmax);

        double err[2];
        godunov_error(&gd, err);
        double errtot = err[0]+err[1];
        printf("N = %d,    dx = %f,    erreur L1 = %f\n", gd.N, gd.dx, errtot);
        fprintf(fic, "%f %f\n", gd.dx, errtot);

        godunov_free(&gd);
    }

    fclose(fic);
    system("gnuplot plotcomError");
}

/**
 * Calcule l'ordre de convergence estimé
 */
void printEstimatedOrder(int N){
    /* Resolution des differents problemes */
    godunov gd = {0};

    // Constantes d'un prbleme a l'autre
    double xmin = -1;
    double xmax = 1;
    double cfl = 0.5;
    int m = 1;
    double tmax = 0.7;

    godunov_init(&gd, xmin, xmax, cfl, m, N);
    godunov_solve(&gd, tmax);
    double err1[2];
    godunov_error(&gd, err1);

    godunov_init(&gd, xmin, xmax, cfl, m, 2*N);
    godunov_solve(&gd, tmax);
    double err2[2];
    godunov_error(&gd, err2);
    godunov_free(&gd);

    double errtot1 = err1[0]+err1[0];
    double errtot2 = err1[1]+err1[1];
    double order = log(errtot1 / errtot2) / log(2);
    printf("\nOrdre de convergence estimé: %f\n\n", order);
}
