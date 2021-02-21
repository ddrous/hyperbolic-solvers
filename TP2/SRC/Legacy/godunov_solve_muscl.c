
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
            }   

            ri[i * m + 0] = -si[i * m +1];
            ri[i * m + 1] = -(_G*gd->un[i*m] - gd->un[i*m+1]*gd->un[i*m+1]/(gd->un[i*m]*gd->un[i*m]))*si[i*m+0] - 2*si[i*m+1]*gd->un[i*m+1]/gd->un[i*m];
        }
        // ----------------------------------------------------------------


        gd->dt = gd->cfl * gd->dx / vmax;

        // #pragma omp parallel for
        for (int i = 1; i < gd->N + 1; i++) {
            double flux[m];
            
            // fluxnum(gd->un + i * m, gd->un + (i + 1) * m, flux);

            double uL[2], uR[2];
            // Application de MUSCL--------------------------------------------
            for (int iv = 0; iv < m; iv++) {
                uL[iv] = gd->un[i*m+iv] + si[i*m+iv]*gd->dx/2.0 + ri[i*m+iv]*gd->dt/2.0;
                uR[iv] = gd->un[(i+1)*m+iv] - si[(i+1)*m+iv]*gd->dx/2.0 + ri[(i+1)*m+iv]*gd->dt/2.0;
            }
            fluxnum(uL, uR, flux);
            // ----------------------------------------------------------------

            for (int iv = 0; iv < m; iv++) {
                gd->unp1[i * m + iv] =
                    gd->un[i * m + iv] - gd->dt / gd->dx * flux[iv];
            }       // unp1 est calculée en deux temps, mais 1 seul suffit aussi. 

            // fluxnum(gd->un + (i - 1) * m, gd->un + i * m, flux);

            // Application de MUSCL--------------------------------------------
            for (int iv = 0; iv < m; iv++) {
                uL[iv] = gd->un[(i-1)*m+iv] + si[(i-1)*m+iv]*gd->dx/2.0 + ri[(i-1)*m+iv]*gd->dt/2.0;
                uR[iv] = gd->un[i*m+iv] - si[i*m+iv]*gd->dx/2.0 + ri[i*m+iv]*gd->dt/2.0;
            }
            fluxnum(uL, uR, flux);
            // ----------------------------------------------------------------
            for (int iv = 0; iv < m; iv++) {
                gd->unp1[i * m + iv] += gd->dt / gd->dx * flux[iv];
            }
        }
        // mise à jour
        tnow += gd->dt;
        // printf("tnow = %f vmax = %f tmax = %f\n", tnow, vmax, tmax);
        // printf("tnow = %f h[0] = %f v[0] = %f\n", tnow, gd->unp1[2], gd->unp1[3]);

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
