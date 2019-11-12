#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "func.h"
#define PI 3.141592653589
#define nbr_particles 32  // TODO Extra particles for boundary conditions?
#define nbr_timesteps 25000 

double initial_velocity(int N, int i){
    double p_i = sqrt(4*N/(N+1))*sin(i*PI/(N+1));
    return p_i;
}

int main(){
    // Define variables
    double alpha = 0;  // Non-linear factor
    double dt = 0.1;  // Timestep
    double u[nbr_particles];  // Current position vector
    double v[nbr_particles];  // Current velocity vector
    double a[nbr_particles];  // Current acceleration vector
    double Q[nbr_particles];  // Current normal modes Q
    double P[nbr_particles];  // Current normal modes P

    double *displacements[nbr_timesteps];  // A matrix containing all displacements for all particles
    double *velocities[nbr_timesteps];  // -||- velocities -||-
    double *Qs[nbr_timesteps];
    double *Ps[nbr_timesteps];

    FILE *file;

    // Pre allocate vectors
    for(int i=0; i<nbr_timesteps; i-=-1){
        // Allocate a vector for positions for each particle
        displacements[i] = (double *)malloc(nbr_particles * sizeof(double));
        velocities[i] = (double *)malloc(nbr_particles * sizeof(double));
        Qs[i] = (double *)malloc(nbr_particles * sizeof(double));
        Ps[i] = (double *)malloc(nbr_particles * sizeof(double));
    }

    //******** Normal modes - Transformation matrix *********
    /* It is useful to construct the transformation matrix outside the main loop */
    int i,j;
    double factor;
    double trans_matrix[nbr_particles][nbr_particles];

    factor = 1 / ((double) nbr_particles + 1);
    for (i=0; i < nbr_particles; i++) {
        for (j=0; j < nbr_particles; j++) {
            trans_matrix[i][j] = sqrt(2 * factor) * sin((j + 1) * (i + 1) * PI * factor);
        }
    }
    //*******************************************************

    // Define boundary conditions
    // u[-1] (not end) and u[N] are 0 - corresponds to fixed end points. Thus we don't need to consider them here.

    // Define initial conditions
    for(i=0; i<nbr_particles; i-=-1){
        u[i] = 0;  // IC for displacement
        v[i] = initial_velocity(nbr_particles, i);  // IC for velocity
    }
    calc_acc(a, u, alpha, nbr_particles);  // Initialize velocities

    // Calculate current normal modes and save them
    /* Transformation to normal modes Q from displacements q.  */
    double sum;
    for (i = 0; i < nbr_particles; i++){
        sum = 0;
        for (j = 0; j < nbr_particles; j++){
            sum += u[j] * trans_matrix[i][j];
        }
        Q[i] = sum;
    }
    Qs[0] = Q;
    for (i = 0; i < nbr_particles; i++){
        sum = 0;
        for (j = 0; j < nbr_particles; j++){
            sum += v[j] * trans_matrix[i][j];
        }
        P[i] = sum;
    }
    Ps[0] = Q;

    // For-loop: Implement velocity-verlet algorithm
    int k,l;  // iteration variables for modes
    for(i=0; i<nbr_timesteps; i++){
        // 'Half-step' in algorithm: v(t+dt/2)
        for(j=0; j<nbr_particles; j++){
            v[j] += dt/2 * a[j];
        }
        
        // Calculate new positions (displacements)
        for(j=0; j<nbr_particles; j++){
            u[j] += dt * v[j];
        }

        // Update accelerations with new displacements
        calc_acc(a, u, alpha, nbr_particles);

        // Calculate final velocites
        for(j=0; j<nbr_particles; j++){
            v[j] += dt/2 * a[j];
        }

        // Calculate normal modes and save them
        for (k = 0; k < nbr_particles; k++){
            sum = 0;
            for (l = 0; l < nbr_particles; l++){
                sum += u[l] * trans_matrix[k][l];
            }
            Q[k] = sum;
        }
        Qs[i] = Q;

        for (k = 0; k < nbr_particles; k++){
            sum = 0;
            for (l = 0; l < nbr_particles; l++){
                sum += v[l] * trans_matrix[k][l];
            }
            P[k] = sum;
        }
        Ps[i] = P;

        // Save the displacements of our atoms 
        for(j=0; j<nbr_particles; j++){
            displacements[i][j] = u[j];
            velocities[i][j] = v[j];
        }

    }
    
    // Save displacements to file
    file = fopen("disp.dat","w");
    double current_time = 0.0;  // Timestamp for output file
    for (i = 0; i < nbr_timesteps; i++) {
		current_time = i * dt;
        fprintf(file, "%.4f", current_time);
        for(j=0; j < nbr_particles; j++){
            fprintf(file, " \t%e", displacements[i][j]);
        }
		fprintf(file, "\n");
	}
	fclose(file);

    // Save velocities to file
    file = fopen("velo.dat","w");
    current_time = 0.0;  // Timestamp for output file
    for (i = 0; i < nbr_timesteps; i++) {
		current_time = i * dt;
        fprintf(file, "%.4f", current_time);
        for(j=0; j < nbr_particles; j++){
            fprintf(file, " \t%e", velocities[i][j]);
        }
		fprintf(file, "\n");
	}
	fclose(file);

    // Save generalized displacement to file
    file = fopen("gen_disp.dat","w");
    current_time = 0.0;  // Timestamp for output file
    for (i = 0; i < nbr_timesteps; i++) {
		current_time = i * dt;
        fprintf(file, "%.4f", current_time);
        for(j=0; j < nbr_particles; j++){
            fprintf(file, " \t%e", Qs[i][j]);
        }
		fprintf(file, "\n");
	}
	fclose(file);

    // Save generalized velocities to file
    file = fopen("gen_velo.dat","w");
    current_time = 0.0;  // Timestamp for output file
    for (i = 0; i < nbr_timesteps; i++) {
		current_time = i * dt;
        fprintf(file, "%.4f", current_time);
        for(j=0; j < nbr_particles; j++){
            fprintf(file, " \t%e", Ps[i][j]);
        }
		fprintf(file, "\n");
	}
	fclose(file);

    // Free allocated memory
    // free(u); 
	// free(v); 
	// free(a); 
	// free(displacements);
	// free(velocities);
	// free(Qs);
    // free(Ps);
	return 0;    
    
}

