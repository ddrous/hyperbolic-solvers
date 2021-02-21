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

// Vitesse de transport
#define _C 1

int main(void){

    bool etudeCvg = true;       

    if (!etudeCvg){              // On résout un seul problème et on l'affiche
        godunov gd = {0};

        double xmin = -1;
        double xmax = 2;
        double cfl = 0.5;
        int m = 1;
        int N = 2500;
        double tmax = 2;

        godunov_init(&gd, xmin, xmax, cfl, m, N);
        godunov_solve(&gd, tmax);
        godunov_plot(&gd, true);
        // double err = godunov_error(&gd);
        // printf("erreur L1=%f\n", err);
        godunov_free(&gd);
    } else                  // On résout plusieur problemes
    {  
        // Calcul des erreurs L1 pour plusieurs problemes et plot
        int NMin = 10;
        int NMax = 10000;
        int nbSimuMax = 9;          // Nombre maximal de simu a effectuer

        printEstimatedOrder(NMax/100);
        plotError(NMin, NMax, nbSimuMax);
    }
    
    return 0;
}


/*******************************************************************************
 * POUR L'EQUATION DE BURGERS
 ******************************************************************************/

void fluxphy(double* w, double* flux) { flux[0] = w[0] * w[0] / 2; }

double lambda_max(double* u) { return fabs(u[0]); }

void riemann(double* a, double* b, double z, double* w) {
    if (a[0] > b[0]) {
        // choc
        double s = (a[0] + b[0]) / 2;
        if (z < s)
            w[0] = a[0];
        else
            w[0] = b[0];
    } else {
        // onde de détente
        if (z < a[0]){
            w[0] = a[0];
        } else if (z > b[0]){
            w[0] = b[0];
        } else
            w[0] = z;
    }
}


void solexacte(double x, double t, double* w) {
    if (x < t){
        w[0] = 1;
    } else if(x > 1) {
        w[0] = 0;
    } else { // if ( 0 <= (x-t)/(1-t) && (x-t)/(1-t) <=1.)
        w[0] = (1-x)/(1-t);
    }
}


/*******************************************************************************
 * PARTIE COMMUNE AUX DEUX PROBLEMES
 ******************************************************************************/

void fluxnum(double* a, double* b, double* flux) {
    double w[10];
    riemann(a, b, 0., w);
    fluxphy(w, flux);
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
        double si[gd->N+1], ri[gd->N+1];
        for (int i = 1; i < gd->N + 1; i++) {
            double alpha = (gd->un[i] - gd->un[i-1])/gd->dx;
            double beta = (gd->un[i+1] - gd->un[i])/gd->dx;
            double gamma = (gd->un[i+1] - gd->un[i-1])/(2.0*gd->dx);
            si[i] = minmod(alpha, beta, gamma);
            ri[i] = - gd->un[i] * si[i];
        }
        // ----------------------------------------------------------------

        gd->dt = gd->cfl * gd->dx / vmax;

        // #pragma omp parallel for
        for (int i = 1; i < gd->N + 1; i++) {

            double flux[m];

            // fluxnum(gd->un + i * m, gd->un + (i + 1) * m, flux);
            
            // Application de MUSCL--------------------------------------------
            double uL[1] = {gd->un[i] + si[i]*gd->dx/2.0 + ri[i]*gd->dt/2.0};          // de gauche
            double uR[1] = {gd->un[i+1] - si[i+1]*gd->dx/2.0 + ri[i+1]*gd->dt/2.0};    // de droite
            fluxnum(uL, uR, flux);
            // ----------------------------------------------------------------

            for (int iv = 0; iv < m; iv++) {
                gd->unp1[i * m + iv] =
                    gd->un[i * m + iv] - gd->dt / gd->dx * flux[iv];
            }       // unp1 est calculée en deux temps, mais 1 seul suffit aussi. 
            
            // fluxnum(gd->un + (i - 1) * m, gd->un + i * m, flux);

            // Application de MUSCL--------------------------------------------
            uL[0] = gd->un[i-1] + si[i-1]*gd->dx/2.0 + ri[i-1]*gd->dt/2.0;          // de gauche
            uR[0] = gd->un[i] - si[i]*gd->dx/2.0 + ri[i]*gd->dt/2.0;    // de droite
            fluxnum(uL, uR, flux);
            // ----------------------------------------------------------------
            for (int iv = 0; iv < m; iv++) {
                gd->unp1[i * m + iv] += gd->dt / gd->dx * flux[iv];
            }
        }
        // mise à jour
        tnow += gd->dt;
        // printf("tnow = %f vmax = %f tmax = %f\n", tnow, vmax, tmax);

        // conditions aux limites
        int i = 0;
        solexacte(gd->xi[i], tnow, gd->unp1 + i * m);
        //**********************************************************************
        i = gd->N + 1;
        solexacte(gd->xi[i], tnow, gd->unp1 + i * m);   // Petite tricherie, ca ne change rien vu que le shcema converge !
        //**********************************************************************

        memcpy(gd->un, gd->unp1, (gd->N + 2) * m * sizeof(double));
    }
    gd->tfin = tnow;
}

void godunov_plot(godunov* gd, bool visu) {
    FILE* fic = fopen("godu.dat", "w");

    for (int i = 1; i < gd->N + 1; i++) {
        double uex[gd->m];
        double u0[gd->m];     // Solution initiale
        // solexacte(gd->xi[i], gd->tfin, uex);
        solexacte(gd->xi[i], gd->tfin, uex);
        solexacte(gd->xi[i], 0, u0);
        int iv = 0;

        fprintf(fic, "%f %f %f %f\n", gd->xi[i], uex[iv], gd->un[i * gd->m + iv], u0[iv]);
    }

    fclose(fic);

    if (visu) {
        int status = system("gnuplot plotcom");
        assert(status == EXIT_SUCCESS);
    }
}

double godunov_error(godunov* gd) {
    double erreur = 0;
    for (int i = 1; i < gd->N + 1; i++) {
        double uex[gd->m];
        solexacte(gd->xi[i], gd->tfin, uex);
        int iv = 0;

        erreur += gd->dx * fabs(uex[iv] - gd->un[i * gd->m + iv]);
    }
    return erreur;
}

double erreurL1(godunov *gd){
    double erreur = 0;
    for (int i = 1; i < gd->N+1; i++){
        double uex[gd->m];
        solexacte(gd->xi[i], gd->tfin, uex);

        int iv = 0;     // Indice pour parcourir les m equations de transport
        erreur += fabs(uex[iv] - gd->un[i*gd->m]) * gd->dx;
    }
    return erreur;
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
    double xmin = -1;
    double xmax = 2;
    double cfl = 0.5;
    int m = 1;
    double tmax = 1;

    for (int i = 0; i < nbSimu; i++) {      // On n'utilise que les nbSimu premieres valeurs
        int N = pbSizes[i];

        godunov_init(&gd, xmin, xmax, cfl, m, N);
        godunov_solve(&gd, tmax);

        // double err = godunov_error(&gd);
        double err = erreurL1(&gd);
        printf("N = %d,    dx = %f,    erreur L1 = %f\n", gd.N, gd.dx, err);
        // fprintf(fic, "%f %f\n", log(gd.dx), log(err));
        fprintf(fic, "%f %f\n", gd.dx, err);

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
    double xmax = 2;
    double cfl = 0.5;
    int m = 1;
    double tmax = 0.5;

    godunov_init(&gd, xmin, xmax, cfl, m, N);
    godunov_solve(&gd, tmax);
    double err1 = erreurL1(&gd);

    godunov_init(&gd, xmin, xmax, cfl, m, 2*N);
    godunov_solve(&gd, tmax);
    double err2 = erreurL1(&gd);
    godunov_free(&gd);

    double order = log(err1 / err2) / log(2);
    printf("\nOrdre de convergence estimé: %f\n\n", order);
}