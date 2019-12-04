#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#define PI 3.141592653589793238

double Fext(double x, double m, double f0){
    return -x*m*(2*PI*f0)*(2*PI*f0);
}

void brownianDynamics(int N, int N_walkers, double dt, double relTime, double pos[][N_walkers], double vel[][N_walkers], int col){

    /* Variable declarations */
    double x = 0;  // Current position
    double v = 0;  // Current velocity
    double a = 0;   // Current acceleration
    double vt = 0;  // Intermediate velocity

    /* Trap and particle parameters */
    double eta = 1/relTime;  // Viscosity - relTime in ms
    double f0 = 3;  // Natural frequency of the trap - 3 kHz - 1 Hz = 1/s = 1/(1000 ms) => 3 kHz = 3000/(1000 ms) = 3/ms
    double m = 3.013e-11;  // Mass of silica particle in our units - M = 1g
    double T = 297;  // Temperature in Kelvin
    double kB = 1.381e-20;  // Boltzmann constant in our units
    double vth = sqrt(kB*T/m);
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
    a = Fext(x, m, f0)/m;
    
    
    pos[0][col] = x; 
    vel[0][col] = v;
    // printf("Time = %.4ld \n", time(NULL));
    // printf("x0 = %.4f \n", x);
    
    
    g = gsl_ran_gaussian(q, 1.0);

    /* Implement algorithm */
    for(int i=1; i<N; i++){
        /* Calculate intermediate velocity */
        g = gsl_ran_gaussian(q, 1.0);  // Gaussian random number with unit variance
        v = 0.5*a*dt + sqrt(c0)*v + vth*sqrt(1-c0)*g;

        /* Calculate new position */
        x = x + v*dt;

        /* Calculate new external acceleration */
        a = Fext(x, m, f0)/m;

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
    int N = 200;  // Number of algorithm steps
    int N_walkers = 5;  // Number of walkers - times the algorithm will be started 
    double dt = 0.01;

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

    /* Task 3 - run Brownian Dynamics simulation */
    double tau = 48.5; 
    for(int i=0; i<N_walkers; i++){
        brownianDynamics(N, N_walkers, dt, tau, pos, vel, i);
        // printf("pos[0][%d]: %.4f \n", i, pos[0][i]);
    }
    
    /* Write to file */
    f = fopen("pos.dat", "w");
    for(int i=0; i<N; i++){
        fprintf(f, "%.8f \t", dt*i);
        for(int j=0; j<N_walkers; j++){
            fprintf(f, "%.8f \t", pos[i][j]);
        }
        fprintf(f, "\n");
    }
    fclose(f); 
    f = fopen("vel.dat", "w");
    for(int i=0; i<N; i++){
        fprintf(f, "%.8f \t", dt*i);
        for(int j=0; j<N_walkers; j++){
            fprintf(f, "%.8f \t", vel[i][j]);
        }
        fprintf(f, "\n");
    }
    fclose(f); 

    free(pos); pos=NULL;
    free(vel); vel=NULL;
}