/*
    E2_func.c 
*/
#include <math.h>

void calc_acc(double *a, double *u, double alpha, double kappa, double m, int nbr_particles){
    /* Calculating the acceleration on the boundaries */
    // The boundary displacements for u[-1] and u[nbr_particles] are both 0 due to our BCs.
    a[0] = (kappa*(u[1]-2*u[0]) + alpha*pow(u[1]-u[0], 2.0) - alpha*pow(u[0]-0, 2.0))/m;  // u_-1 = 0
    a[nbr_particles-1] = (kappa*(u[nbr_particles-2]-2*u[nbr_particles-1]) + alpha*pow(0-u[nbr_particles-1],2.0) - alpha*pow(u[nbr_particles-1]-u[nbr_particles-2],2.0) )/m;  //u_N+1 = 0

    for(int i=1; i<nbr_particles-1; i++){
        // Note that we only update particles 1,...,N-1, not 0 and N+1.
        a[i] = (kappa*(u[i+1]+u[i-1]-2*u[i]) + alpha*pow(u[i+1]-u[i],2.0) - alpha*pow(u[i]-u[i-1],2.0))/m;
    }
}