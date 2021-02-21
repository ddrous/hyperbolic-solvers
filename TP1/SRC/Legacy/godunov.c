
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "godunov.h"

// vitesse de transport
#define _C 1

int main() {

    godunov gd = {0};

    double xmin = -1;
    double xmax = 1;
    double cfl = 0.99;

    int m = 1;
    int N = 1000;

    double tsol = 0.25;         // Temps final

    godunov_init(&gd, xmin, xmax, cfl, m, N);

    godunov_solve(&gd, tsol);

    godunov_plot(&gd, tsol);


    // printf("log(dx) = %f \n", log(gd.dx));
    // printf("log(erreurL1) = %f", log(erreurL1(&gd, tsol)));
    // plotError(&gd, tsol);

    printf("\n");
    return 0;
}



void fluxnum(double *a,double *b, double *flux){
    flux[0] = _C * a[0];

}

// Il faut arranger ca car w est supose contenir m flotants
void solexacte(double x, double t, double *w){
    // double uL = 1; Ancienne valeur
    double uL = exp(-t + x/_C);
    double uR = 0;

    if (x < _C * t){
        w[0] = uL;
    } else
    {
        w[0] = uR;
    }
    
}


void godunov_init(godunov *gd, double xmin, double xmax, double cfl, int m, int N){
    gd->xmin = xmin;
    gd->xmax = xmax;
    gd->m = m;
    gd->N = N;
    gd->cfl = cfl;
    gd->dx = (xmax - xmin) / N;
    gd->dt = 0;

    gd->xi = malloc((N+2) * sizeof(double));
    gd->un = malloc((N+2) * sizeof(double) * m);    // Pour un systeme de Friedrich, contenant m equations
    gd->unp1 = malloc((N+2) * sizeof(double) * m);

    for (int i = 0; i < N+2; i++){
        gd->xi[i] = xmin + gd->dx / 2 + (i-1)*gd->dx;
        double t = 0;
        // Condition initiale
        solexacte(gd->xi[i], t, gd->un+i*m);        // pointeur decallé. On prend les truc de m en m
        // solexacte(gd->xi[i], t, gd->un+i + (m-1)*(N+2));        // je pense que c'est mieux comme ca, mais FAUX !!!!!!! Car on stocke dábord tous les [m] avant d'attacquer le prochain i 
    }

}

double riemann(double uL, double uR, double xSurt){
    if (xSurt < _C)
        return uL;
    else
        return uR;
}

// Fonction sur laquelle on travaille
double f(double u){
    return _C*u;
}

double fluxGodunov(double a, double b){
    return f(riemann(a, b, 0));
}

// void updateUnp1(godunov *gd, double *unp1, double *un){
//     for (int j = 0; j < gd->m; j++){
//         unp1[j] = un[j] - (gd->dt / gd->dx) * (fluxGodunov(gd->un[(i)*m], gd->un[(i+1)*m]) - fluxGodunov(gd->un[(i-1)*m], gd->un[(i)*m]));
//     }
// }


double lambda_max(double *u){
    return _C;          // tres simple car il s'agit de l'equation de transport
}


void godunov_solve(godunov *gd, double tmax){
    double tnow = 0;
    int m = gd->m;

    while (tnow < tmax){
        double vmax = 0;
        // calcul de la vitesse max
        for (int i = 0; i < gd->N+2; i++){
            double vloc = lambda_max(gd->un + gd->m * i);
            vmax = vmax > vloc? vmax : vloc;
        }
        
        gd->dt = gd->cfl * gd->dx / vmax;
        // gd->dt = gd->cfl * gd->dx / (2*vmax);       // D'apres le cours

        for (int i = 0; i < gd->N+2; i++){
            gd->unp1[i*m] = gd->un[i*m] - (gd->dt / gd->dx) * (fluxGodunov(gd->un[(i)*m], gd->un[(i+1)*m]) - fluxGodunov(gd->un[(i-1)*m], gd->un[(i)*m]));
        }

    // printf("\n");
    // for (int i = 0; i < gd->N+2; i++){
    //     printf("sol = %f \n", gd->un[i*m]);
    // }

        // Condition de Dirichlet
        solexacte(gd->xi[0], tnow, gd->un);        // Pour aller plus vite, on sait que c'est egale a la solution exacte a gauche 

        for (int i = 1; i < (gd->N+2)*m; i++){
            gd->un[i] = gd->unp1[i];
        }
        tnow += gd->dt;
        // printf("Nouvelle iteration");

    }
    
}


void godunov_plot(godunov *gd, double tsol){
    FILE * fic = fopen("godu.dat", "w");

    for (int i = 1; i < gd->N+1; i++){
        double uex[gd->m];
        solexacte(gd->xi[i], tsol, uex);
        int iv = 0;

        fprintf(fic, "%f %f %f\n", gd->xi[i], uex[iv], gd->un[i]);
    }

    fclose(fic);

    system("gnuplot plotcom");
    
}

double erreurL1(godunov *gd, double tsol){
    double erreur = 0;
    for (int i = 1; i < gd->N+1; i++){
        double uex[gd->m];
        solexacte(gd->xi[i], tsol, uex);

        int iv = 0;     // Indice pour parcourir les m equations de transport
        erreur += fabs(uex[iv] - gd->un[i*gd->m]);
    }

    return erreur;
}

void plotError(godunov *gd, double tsol){
    FILE * fic = fopen("erreur.dat", "a");

    fprintf(fic, "%f %f\n", log(gd->dx), log(erreurL1(gd, tsol)));

    fclose(fic);

    system("gnuplot plotcomError");
    
}
