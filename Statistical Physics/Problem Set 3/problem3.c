#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <string.h>
#include <gsl/gsl_rng.h>


double magnetization(int N, double lat[][N]){
    /* The average magnetization is the sum of all spins, divided by the number of particles*/
    double m = 0;

    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            m += lat[i][j];
        }
    }
    printf("Magnetization: %.2f\n", m);
    return m / (double)(N*N);
}

double energy(int N, double lat[][N]){
    /* Calculation the energy of the current lattice - implementing periodic boundary conditions */
    double E = 0; 
    double neigh_1;
    double neigh_2;
    double neigh_3;
    double neigh_4;


    /* Iterate over whole lattice, check neighbors */
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            /* Check neighbors */
            /* Periodic boundary conditions */
            if(i==N-1){
                neigh_1 = lat[i][j]*lat[0][j] + 1.0;  // -> Convert to +2 or 0 from 1 or -1
                neigh_2 = lat[i][j]*lat[i-1][j] + 1.0;  // <-
            }else if(i==0){
                neigh_1 = lat[i][j]*lat[i][j] + 1.0;  
                neigh_2 = lat[i][j]*lat[N-1][j] + 1.0; 
            }else{
                neigh_1 = lat[i][j]*lat[i+1][j] + 1.0;  
                neigh_2 = lat[i][j]*lat[i-1][j] + 1.0; 
            }

            if(j==N-1){
                neigh_3 = lat[i][j]*lat[i][0] + 1.0;  // ^
                neigh_4 = lat[i][j]*lat[i][j-1] + 1.0;  // down
            }else if(j==0){
                neigh_3 = lat[i][j]*lat[i][j+1] + 1.0;  
                neigh_4 = lat[i][j]*lat[i][N-1] + 1.0; 
            }else{
                neigh_3 = lat[i][j]*lat[i][j+1] + 1.0;  
                neigh_4 = lat[i][j]*lat[i][j-1] + 1.0;  
            }
            
            E += (neigh_1 + neigh_2 + neigh_3 + neigh_4) / 2.0;  // Transform back by dividing by two
        }
    }
    return E;
}


void metropolis(int N, double beta, double J, double tol, double lat[][N], int M, double m[]){
    /* Initialize variables */
    double (*prop_lat)[N] = malloc(sizeof(double [N][N]));  // lattice with one proposal spin flip
    double prev_E;  // Previous energy in deltaE
    double new_E;  // New energy in deltaE
    double delta_E;
    double prev_m;
    double new_m;
    gsl_rng *r = gsl_rng_alloc (gsl_rng_taus);  // Init random generator
    gsl_rng_set(r, time(NULL));
    int i;
    int j;
    double eps;
    int steps = 1;
    int acc_steps = 0;

    /* Calculate energy of current lattice */
    prev_E = energy(N, lat);
    prev_m = -DBL_MAX;
    new_m = magnetization(N, lat);

    while(fabs(new_m - prev_m) > tol && steps < 10){
        if(steps%1000 == 0){
            printf("Iterations: %d \n", steps);
        }
        /* Randomize indices i and j to change */
        i = gsl_rng_uniform_int(r, N);
        j = gsl_rng_uniform_int(r, N);

        /* Flip spin at position i, j */
        printf("Spin flipped at i=%d, j=%d \n", i, j);
        lat[i][j] = -1.0 * lat[i][j]; 
    
        /* Calculate new energy */
        new_E = energy(N, lat);

        /* Calculate energy difference */
        delta_E = new_E - prev_E;
        

        if(delta_E < 0){
            /* Accept change */
            prev_E = new_E;
            acc_steps++;
            /* Calculate magnetization */
            prev_m = new_m;  // save old value
            new_m = magnetization(N, lat);

        }else{
            /* Accept change if e^{beta*\deltaE} > eps*/
            eps = gsl_rng_uniform(r);
            if (exp( -delta_E * beta) > eps){
                /* Accept change */
                prev_E = new_E;
                acc_steps++;
                /* Calculate magnetization */
                prev_m = new_m;  // save old value
                new_m = magnetization(N, lat);
            }else{
                /* Revert lattice */
                lat[i][j] = -1 * lat[i][j]; 
            }
        }

        m[steps] = new_m;
        // printf("magnetization: %.2f \n", new_m);
        steps++;

        /* Print lattice - for debug */
        // for(int k=0; k<N; k++){
        //     for(int l=0; l<N; l++){
        //         printf("%.1f \t", lat[k][j] );
        //     }
        //     printf("\n");
        // }

        if(steps==M){
            /* Reallocate more space in m-array*/
            printf("Reallcate array length \n");
            M *= 2;
            m = realloc(m, M * sizeof(double));  // Doubles length of array
        }
    }
    
    /* Report acceptance ratio */
    printf("Acceptance ration: %.4f \n", (double)acc_steps/steps);
    

    /* Free variables */
    free(prop_lat); prop_lat=NULL;
    gsl_rng_free(r);
}

void simulate_lattice(int N, double tol, double T,  double kB,  double J){
    /* RNG */
    srand(time(NULL)); // Set the seed for rand
    double rand_disp;  // random displacement

    /* Variables */
    double beta = kB*T;
    double (*lat)[N] = malloc(sizeof(double [N][N]));

    int M = 1000;
    double *m = malloc(M * sizeof(double));  // Average magnetizations for each iteration



    FILE *f;  // Output file
    char filename[50] = "";
    char dir[100] = "problem3_output/";
    char Temperature[50];
    char output[50];
    gcvt(T, 10, Temperature);
    strcat(dir, Temperature);

    /* Initialize lattice randomly */
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            rand_disp = round((double)rand() / (double)RAND_MAX);  // -1 or +1
            if(rand_disp == 0.0){
                rand_disp = -1.0;
            }
            lat[i][j] = rand_disp;
        }
    }

    // for(int k=0; k<N; k++){
    //     for(int l=0; l<N; l++){
    //         printf("%.1f \t", lat[k][l] );
    //     }
    //     printf("\n");
    // }

    // Run metropolis to tol
    metropolis(N, beta, J, tol, lat, M, m);

    // Save results
    filename[0] = '\0';  // Empty filename string
    strcat(filename, dir);
    strcat(filename, "_lattice.dat");
    f = fopen(filename, "w");
    for (int i = 0; i < N; i++)
    {
        for (int j=0; j<N; j++){
            fprintf(f, "%.1f \t", lat[i][j]);
        }
        fprintf(f, "\n");
    }
    fclose(f);

    filename[0] = '\0';  // Empty filename string
    strcat(filename, dir);
    strcat(filename, "_magnetization.dat");
    f = fopen(filename, "w");
    for (int i = 0; i < M; i++)
    {
        fprintf(f, "%.16f \n", m[i]);
    }
    fclose(f);




    // Free variables
    free(lat); lat=NULL;
}


int main(){
    // Define constants: J, kB, tol
    double J = 1;
    double kB = 0.00008617; 
    double tol = 0.00001;  // Tolerance for convergence

    // Define temperatures
    double beta[] = {10, 10, 20};

    // Define sizes N
    int N[] = {2, 4, 6};

    /* Iterate over all temperatures and all sizes */
    simulate_lattice(4, tol, 1, kB, J);




}