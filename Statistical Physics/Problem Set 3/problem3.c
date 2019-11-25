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
    return m / (double)(N*N);
}

double energy(int N, double lat[][N], double J){
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
    return -J*E;
}

double calc_E_site(int N, double lat[][N], double J, int i, int j){
    /* Calculation of the change in energy of the current lattice - implementing periodic boundary conditions */
    double E = 0; 
    double neigh_1;
    double neigh_2;
    double neigh_3;
    double neigh_4;

    clock_t begin = clock();
    
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
    
    E = (neigh_1 + neigh_2 + neigh_3 + neigh_4) / 2.0;  // Transform back by dividing by two
    
    return -J*E;
}


double metropolis(int N, double beta, double J, double tol, double lat[][N], int M, double m[]){
    
    /* Initialize variables */
    double (*prop_lat)[N] = malloc(sizeof(double [N][N]));  // lattice with one proposal spin flip
    double prev_E_site;  // Previous energy in deltaE
    double new_E_site;  // New energy in deltaE
    double delta_E;
    double prev_m;
    double new_m;
    /* Loop variables */
    long int max_iters = 1e12;
    int run_loop = 1;
    /* Average m values*/
    double old_avg_m = DBL_MAX;
    double curr_avg_m = 1;
    double new_avg_m = 0;
    /* Save spins */
    double old_spin;
    double new_spin;
    /* Initialize RNG */
    gsl_rng *r = gsl_rng_alloc (gsl_rng_taus);  // Init random generator
    gsl_rng_set(r, time(NULL));
    int i;
    int j;
    double eps;
    int steps = 1;
    int acc_steps = 0;

    /* Calculate energy of current lattice */
    double delta_m;  // Change in m when flipping spin
    double curr_m = magnetization(N, lat);
    m[0] = curr_m;  // Save current m

    while(run_loop && steps < max_iters){
        i = gsl_rng_uniform_int(r, N);
        j = gsl_rng_uniform_int(r, N);
        
        /* Calculate energy at randomized site */
        prev_E_site = calc_E_site(N, lat, J, i, j);
        
        /* Flip spin at position i, j */
        old_spin = lat[i][j];
        lat[i][j] = -1.0 * lat[i][j];
        
        
        /* Calculate new energy */
        new_E_site = calc_E_site(N, lat, J, i, j);
        
        /* Calculate energy difference */
        delta_E = new_E_site-prev_E_site;
    
        
        eps = gsl_rng_uniform(r);
        if(delta_E < 0 || exp( -1.0 * delta_E * beta) > eps){
            /* Accept change */
            acc_steps++;            
        }else{
            /* Revert lattice */
            lat[i][j] = -1 * lat[i][j];
        }
        /* Save current magnetization */
        new_spin = lat[i][j];
        delta_m = (new_spin - old_spin)/(double)(N*N);
        curr_m += delta_m;
        m[steps%M] = curr_m;
        // printf("magnetization: %.2f \n", new_m);
        steps++;

        if(steps%M==0 && steps != 0){
            /* Reallocate more space in m-array*/
            // printf("\nCalculate new average value of m array length \n");
            new_avg_m = 0;
            for(int k=0; k<M; k++){
                new_avg_m += m[k];
            }
            new_avg_m /= M;
            // printf("Old m_avg: %.6f \t Curr m_avg: %.6f \t new m_avg %.6f \n", old_avg_m, curr_avg_m, new_avg_m);
            old_avg_m = curr_avg_m;
            curr_avg_m = new_avg_m;
            /* Calculate break condition */
            if(fabs((curr_avg_m - old_avg_m)/old_avg_m) > tol){
                /* Continue loop*/
                run_loop = 1;
                memset(m, 0, sizeof(double) *M);
            }else{
                /* Stop loop */
                run_loop = 0;
                
            }
        }
    }
    
    /* Report acceptance ratio */
    printf("Acceptance ration: %.4f - ", (double)acc_steps/steps);
    if(steps-max_iters != 0){
        printf("Converged in %d steps - ", steps);
    }else{
        printf("Timed out - ");
    }

    /* Free variables */
    free(prop_lat); prop_lat=NULL;
    gsl_rng_free(r);

    /* Return last magnetization */
    return new_avg_m;
}

void simulate_lattice(int N, double tol, double T[], int Ts, double kB,  double J){
    /* RNG */
    srand(time(NULL)); // Set the seed for rand
    double rand_disp;  // random displacement

    /* Variables */
    // double beta = kB*T;
    double (*lat)[N] = malloc(sizeof(double [N][N]));
    double *m = malloc(Ts * sizeof(double));  // Final magnetization for each temperature
    int M = 10000;
    double *m_final = malloc(M * sizeof(double));  // Save last 10000 magnetization values

    FILE *f;  // Output file
    char filename[200] = "";
    char dir[100] = "problem3_output/";
    char subfolder[50] = "N=";
    char size[50];
    char Temperature[50];
    char output[50];

    gcvt(N,3, size);
    strcat(subfolder, size);
    strcat(subfolder, "/");
    strcat(dir, subfolder);
    // printf("Dir: %s \n", dir);

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
    clock_t begin;
    clock_t end; 
    double time_spent;

    for(int i=0; i<Ts; i++){
        double beta = 1.0 / (T[i] * kB);
        printf("\tTemperature: %.2f K - ", T[i]);

        begin=clock();

        // Run metropolis
        m[i] = metropolis(N, beta, J, tol, lat, M, m_final);
        
        end = clock();
        time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
        printf("Time: %.2f s \n", time_spent);
        
        /* Save results*/
        gcvt(T[i],3, Temperature);
        char T_lattice[50] = "lattices/T=";
        char T_traces[50] = "traces/T=";

        strcat(T_lattice, Temperature);
        strcat(T_traces, Temperature);

        // Save lattice for this temperature
        filename[0] = '\0';  // Empty filename string
        strcat(filename, dir);
        strcat(filename, T_lattice);
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
        
        // Save trace over last 10000 timesteps
        filename[0] = '\0';  // Empty filename string
        strcat(filename, dir);
        strcat(filename, T_traces);
        strcat(filename, "_magnetization_trace.dat");
        f = fopen(filename, "w");
        for (int i = 0; i < M; i++)
        {
            fprintf(f, "%.16f \n", m_final[i]);
        }
        fclose(f);
        
    }
    // Save results
    filename[0] = '\0';  // Empty filename string
    strcat(filename, dir);
    strcat(filename, "magnetization.dat");
    f = fopen(filename, "w");
    for (int i = 0; i < Ts; i++)
    {
        fprintf(f, "%.16f \t %.16f \n", T[i], m[i]);
    }
    fclose(f);
    
    // Free variables
    free(lat); lat=NULL;
}


int main(){
    // Define constants: J, kB, tol
    double J = 1;
    // double kB = 0.00008617; 
    double kB = 1;
    double tol = 1e-5;  // Tolerance for convergence

    // Define temperatures
    int Ts = 100;  // Nbr of temperatures
    double Tmin = 0;
    double Tmax = 3;
    double dT = (Tmax-Tmin)/(double)Ts;  // Temperature step
    double *T = malloc(Ts * sizeof(double));
    for(int i=0; i<Ts; i++){
        T[i] = i*dT; 
    }

    // Define sizes N
    int Ns = 1;
    int N[] = {512};

    /* Iterate over all temperatures and all sizes */
    clock_t begin = clock();

    for(int i=0; i<Ns; i++){
        printf("N=%d \n", N[i]);
        simulate_lattice(N[i], tol, T, Ts, kB, J);
    }
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Total execution time: %.2f \n s", time_spent);
}