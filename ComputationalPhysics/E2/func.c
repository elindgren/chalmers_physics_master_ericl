/*
    E2_func.c 
*/

void calc_acc(double *a, double *u, double alpha, int nbr_particles){
    /* Calculating the acceleration on the boundaries */
    // The boundary displacements for u[-1] and u[nbr_particles] are both 0 due to our BCs.
    a[0] = (u[1]-u[0]) + alpha*(u[1]-u[0]);
    a[nbr_particles-1] = (-u[nbr_particles-1]) + alpha*(-u[nbr_particles-1]);

    for(int i=1; i<nbr_particles-1; i++){
        // Note that we only update particles 1,...,N-1, not 0 and N+1.
        a[i] = (u[i+1]+-u[i]) + alpha*(u[i+1]-u[i]);
    }
}