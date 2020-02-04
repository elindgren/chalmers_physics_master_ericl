#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#define PI 3.141592653589793238

double calc_a(double x, double f0){
    return -x*(2*PI*f0)*(2*PI*f0);
}

void brownianDynamics(int N, int N_walkers, double dt, double relTime, double pos[][N_walkers], double vel[][N_walkers], int col){

    /* Variable declarations */
    double x = 0;  // Current position
    double v = 0;  // Current velocity
    double a = 0;   // Current acceleration

    /* Trap and particle parameters */
    double eta = 1.0/relTime;  // Viscosity - relTime in ms
    double f0 = 3e-3;  // Natural frequency of the trap - 3 kHz - 1 Hz = 1/s = 1/(1000 ms) => 3 kHz = 3000/(1000 ms) = 3/ms
    double m = 1;  // Mass of silica particle in our units - M = 1g
    double kBT = 1.700968e-8;  // Boltzmann constant in our units
    double vth = sqrt(kBT/m);
    double c0 = exp(-eta*dt);

    /* RNG */
	const gsl_rng_type *Ty;
	gsl_rng *q;
	gsl_rng_env_setup();
	Ty = gsl_rng_default;
	q = gsl_rng_alloc(Ty);
	gsl_rng_set(q,col);
    double g;  // Gaussian random variable

    /* Initialize system */
    x = 0.1;
    v = 2.0e-3;
    a = calc_a(x, f0);
    
    
    pos[0][col] = x; 
    vel[0][col] = v;
    // printf("Time = %.4ld \n", time(NULL));
    // printf("vth = %.8f \n", vth);
    
    
    g = gsl_ran_gaussian(q, 1.0);

    /* Implement algorithm */
    for(int i=1; i<N; i++){
        /* Calculate intermediate velocity */
        g = gsl_ran_gaussian(q, 1.0);  // Gaussian random number with unit variance
        // printf("last term = %.8f \n", 1/vth*sqrt(1-c0)*g);
        v = 0.5*a*dt + sqrt(c0)*v + vth*sqrt(1-c0)*g;

        /* Calculate new position */
        x = x + v*dt;

        /* Calculate new external acceleration */
        a = calc_a(x, f0);

        /* Calculate new velocity */ 
        g = gsl_ran_gaussian(q, 1.0);
        v = 0.5*sqrt(c0)*a*dt + sqrt(c0)*v + vth*sqrt(1-c0)*g;

        /* Save velocities */
        pos[i][col] = x; 
        vel[i][col] = v;
    }

    /* Write to file */
    gsl_rng_free(q);
    
}

int main(){
    /* Variable declarations */
    int N = 12000;  // Number of algorithm steps
    int N_walkers = 5000;  // Number of walkers - times the algorithm will be started 
    double dt = 0.1;
    int N_times = 5;
    int N_completed = 0;
    int steps_between_save = N/N_times;

    /* RNG */
	const gsl_rng_type *Ty;
	gsl_rng *q;
	gsl_rng_env_setup();
	Ty = gsl_rng_default;
	q = gsl_rng_alloc(Ty);
	gsl_rng_set(q,time(NULL));

    /* File IO */
    FILE *f;
    
    double (*pos)[N_walkers] = malloc(sizeof(double[N][N_walkers]));
    double (*vel)[N_walkers] = malloc(sizeof(double[N][N_walkers]));
    double (*vel_many)[N_walkers] = malloc(sizeof(double[N_times][N_walkers])); 
    double (*pos_many)[N_walkers] = malloc(sizeof(double[N_times][N_walkers])); 

    /* Task 3 - run Brownian Dynamics simulation */
    double tau = 147.3; 
    #pragma omp parallel for
    for(int i=0; i<N_walkers; i++){
        if(i%100 == 0){
            printf("Progress: %.2f \n", (double)N_completed/N_walkers * 100);
        }
        brownianDynamics(N, N_walkers, dt, tau, pos, vel, i);
        N_completed++;
    }   

    for(int i = 0; i<N_times; i++){
        for(int j = 0; j<N_walkers; j++){
            // printf("Time: %.d \n", i*steps_between_save);
            vel_many[i][j] = vel[i*steps_between_save][j];
            pos_many[i][j] = pos[i*steps_between_save][j];
        }
    }
    
    /* Write to file */
    // f = fopen("pos.dat", "w");
    // for(int i=0; i<N; i++){
    //     fprintf(f, "%.8f \t", dt*i);
    //     for(int j=0; j<N_walkers; j++){
    //         fprintf(f, "%.8f \t", pos[i][j]);
    //     }
    //     fprintf(f, "\n");
    // }
    // fclose(f); 
    // f = fopen("vel.dat", "w");
    // for(int i=0; i<N; i++){
    //     fprintf(f, "%.8f \t", dt*i);
    //     for(int j=0; j<N_walkers; j++){
    //         fprintf(f, "%.8f \t", vel[i][j]);
    //     }
    //     fprintf(f, "\n");
    // }
    // fclose(f); 
    f = fopen("pos_many.dat", "w");
    for(int i=0; i<N_times; i++){
        fprintf(f, "%.8f \t", dt*i*steps_between_save);
        for(int j=0; j<N_walkers; j++){
            fprintf(f, "%.8f \t", pos_many[i][j]);
        }
        fprintf(f, "\n");
    }
    fclose(f); 
    f = fopen("vel_many.dat", "w");
    for(int i=0; i<N_times; i++){
        fprintf(f, "%.8f \t", dt*i*steps_between_save);
        for(int j=0; j<N_walkers; j++){
            fprintf(f, "%.8f \t", vel_many[i][j]);
        }
        fprintf(f, "\n");
    }
    fclose(f); 

    free(pos); pos=NULL;
    free(vel); vel=NULL;
}