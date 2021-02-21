#ifndef _GODUNOV_H
#define _GODUNOV_H
#include <stdbool.h>

typedef struct godunov {

    int m; // nombre d'equations de transport a resoudre
    int N; // nombre de cellules  (qu'on va calculer a partir de L et dx)
    double dt, dx;
    double tfin; // temps final
    double cfl; // rppport vmax dt / dx
    double xmin, xmax; // nbornes de l'intervalle

    double *xi; // tableau des centres des mileiux des cellules

    double *un; // solution a l'instant n
    double *unp1; // solution a l'instant n+1
} godunov;

/// flux numerique
void fluxnum(double *a,double *b, double *flux);        // Le resultats est passe dans flux )d'ou le void). a et b sont les parametres, 

// une fonction qui calcule la vitesse max
double lambda_max(double *u);

// initialisation de la structure godunov
void godunov_init(godunov *gd, double xmin, double xmax, double cfl, int m, int N);     // gd c'est la resolution


// resolution de
void godunov_solve(godunov *gd, double tmax); // temps de resolution

// tracer la solution
void godunov_plot(godunov *gd, bool visu); // on trace a tsol, a un instant donne 

// solution exacte de l'equation de transport
void solexacte(double x, double t, double *w);

// Caclcul de l'erreur
void godunov_error(godunov* gd, double *err);

void godunov_free(godunov* gd);

// Mon erreur
// double erreurL1(godunov *gd, double tsol);

// Mon erreur L1
void plotError(int NMin, int NMax, int nbSimu);

void printEstimatedOrder(int N);


#endif