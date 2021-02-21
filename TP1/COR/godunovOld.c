// compilation avec:
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

// vitesse de transport
#define _C 1

int main(void) {
    godunov gd = {0};

    double xmin = -1;
    double xmax = 1;
    double cfl = 0.5;
    int m = 1;
    int N = 5000;
    double tmax = 0.7;

    godunov_init(&gd, xmin, xmax, cfl, m, N);
    godunov_solve(&gd, tmax);
    bool visu = false;
    godunov_plot(&gd, visu);
    // double err = godunov_error(&gd);
    // printf("erreur L1=%f\n", err);
    godunov_free(&gd);
    
    // Calcul des erreurs L1 pour plusieurs problemes et plot de 
    plotError(8);
}

void fluxphy(double* w, double* flux) { flux[0] = _C * w[0]; }

void riemann(double* a, double* b, double z, double* w) {
    if (z < _C) {       // z = x/t
        w[0] = a[0];
    } else {
        w[0] = b[0];
    }
}

void fluxnum(double* a, double* b, double* flux) {
    double w[10];
    riemann(a, b, 0., w);
    fluxphy(w, flux);
}

void solexacte(double x, double t, double* w) {
    // double uL = 1;
    // double uR = 0;

    double uR = 0;
    double uL = exp(-(t - x / _C));

    if (x < _C * t) {
        w[0] = uL;
    } else {
        w[0] = uR;
    }

    // w[0] = cos(x - _C * t);
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
        solexacte(gd->xi[i], t, gd->un + i * m);
    }
}

void godunov_free(godunov* gd) {
    free(gd->xi);
    free(gd->un);
    free(gd->unp1);
}

double lambda_max(double* u) { return _C; }

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

        gd->dt = gd->cfl * gd->dx / vmax;

        // #pragma omp parallel for
        for (int i = 1; i < gd->N + 1; i++) {
            double flux[m];
            fluxnum(gd->un + i * m, gd->un + (i + 1) * m, flux);
            for (int iv = 0; iv < m; iv++) {
                gd->unp1[i * m + iv] =
                    gd->un[i * m + iv] - gd->dt / gd->dx * flux[iv];
            }
            fluxnum(gd->un + (i - 1) * m, gd->un + i * m, flux);
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
        i = gd->N + 1;
        solexacte(gd->xi[i], tnow, gd->unp1 + i * m);

        memcpy(gd->un, gd->unp1, (gd->N + 2) * m * sizeof(double));
    }
    gd->tfin = tnow;
}

void godunov_plot(godunov* gd, bool visu) {
    FILE* fic = fopen("godu.dat", "w");

    for (int i = 0; i < gd->N + 2; i++) {
        double uex[gd->m];
        solexacte(gd->xi[i], gd->tfin, uex);
        int iv = 0;

        fprintf(fic, "%f %f %f\n", gd->xi[i], uex[iv], gd->un[i * gd->m + iv]);
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

void plotError(int nbSimu){
    FILE * fic = fopen("erreur.dat", "w");

    int pbSizes[] = {10, 50, 100, 500, 1000, 5000, 10000, 20000, 20000};
    godunov gd = {0};

    // Constantes d'un prbleme a láutre
    double xmin = -1;
    double xmax = 1;
    double cfl = 0.5;
    int m = 1;
    double tmax = 0.7;

    // Pour l'estimation de l'ordre
    double err1;
    double err2;

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

        // Estimation de l'erreur
        if (N == 5000)
            err1 = err;
        if (N == 2*5000)
            err2 = err;
    }

    double order = log(err1 / err2) / log(2);
    printf("Ordre estimé: %f\n", order);

    fclose(fic);

    system("gnuplot plotcomError");
}
