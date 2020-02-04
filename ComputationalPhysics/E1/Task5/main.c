/*
E1_main.c
 
Created by AL on 2013-10-24.
Further developed by Martin Gren on 2014-10-20.
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "func.h"
#include "fft_func.h"
#define PI 3.141592653589
#define nbr_of_timesteps 32767 /* nbr_of_timesteps+1 = power of 2, for best speed */
#define nbr_of_particles 3 /* The number of particles is 3 */

/* Main program */
int main()
{
	/* Declartion of variables */
	double timestep;
	int i,j;
	double timestep_sq,current_time;
	double m;
	double kappa;

	/* declare file variable */
	FILE *file;

	/* displacement, velocity and acceleration */
	double q[nbr_of_particles];
	double v[nbr_of_particles];
	double a[nbr_of_particles]; 


	/* Allocating memory for large vectors */
	/* displacements for writing to file */
	double *q_1 = malloc((nbr_of_timesteps+1) * sizeof (double));
	double *q_2 = malloc((nbr_of_timesteps+1) * sizeof (double));
	double *q_3 = malloc((nbr_of_timesteps+1) * sizeof (double));
	double *E_pot = malloc((nbr_of_timesteps+1) * sizeof(double));
	double *E_kin = malloc((nbr_of_timesteps+1) * sizeof(double));
	double *E_tot = malloc((nbr_of_timesteps+1) * sizeof(double));

	/* Set variables */
	timestep = 0.01;
	m = 1.0;
	kappa = 1.0;
	timestep_sq = timestep * timestep;


	
	/* Initial conditions */
	/* Set initial displacements and velocites */
	q[0] = 0.1;
	v[0] = 0;

	for (i = 1; i < nbr_of_particles; i++) {
		q[i] = 0;
		v[i] = 0;
	}
	q_1[0] = q[0];
	q_2[0] = q[1];
	q_3[0] = q[2];
	
	// Calculate initial energies for all particles and save them
	double u_particles[] = {q_1[0], q_2[0], q_3[0]};
	double v_particles[] = {0.0, 0.0, 0.0};
	E_pot[0] = calc_pe(u_particles, kappa, nbr_of_particles);
	E_kin[0] = calc_ke(v_particles, nbr_of_particles, m);
	E_tot[0] = E_pot[0] + E_kin[0];


	/* Calculate initial accelerations based on initial displacements */
	calc_acc(a, q, m, kappa, nbr_of_particles);

	/* timesteps according to velocity Verlet algorithm */
	for (i = 1; i < nbr_of_timesteps + 1; i++) {
		/* v(t+dt/2) */
		for (j = 0; j < nbr_of_particles; j++) {
		    v[j] += timestep * 0.5 * a[j];
		} 

		/* q(t+dt) */
		for (j = 0; j < nbr_of_particles; j++) {
		    q[j] += timestep * v[j];
		}

		/* a(t+dt) */
		calc_acc(a, q, m, kappa, nbr_of_particles);

		/* v(t+dt) */
		for (j = 0; j < nbr_of_particles; j++) {
		    v[j] += timestep * 0.5 * a[j];
		} 

		/* Save the displacement of the three atoms */
		q_1[i] = q[0];
		q_2[i] = q[1];
		q_3[i] = q[2];

		// Calculate energies and save them
		double u_particles[] = {q_1[i], q_2[i], q_3[i]};
		double v_particles[] = {v[0], v[1], v[2]};
		E_pot[i] = calc_pe(u_particles, kappa, nbr_of_particles);
		E_kin[i] = calc_ke(v_particles, nbr_of_particles, m);
		E_tot[i] = E_pot[0] + E_kin[0];
	}

	// Calculate powerspectrum
	double *powspect = malloc((nbr_of_timesteps+1)*sizeof(double));
	double *freq = malloc((nbr_of_timesteps+1)*sizeof(double));
	/* make FFT (powerspectrum) */
	powerspectrum(q_1, powspect, nbr_of_timesteps);
	powerspectrum_shift(powspect,nbr_of_timesteps);
	fft_freq_shift(freq, timestep, nbr_of_timesteps);

	/* Print displacement data to output file */
	file = fopen("disp.dat","w");

	for (i = 0; i < nbr_of_timesteps + 1; i++) {
		current_time = i * timestep;
		fprintf(file, "%.4f \t %e \t %e \t %e", current_time, q_1[i], q_2[i], q_3[i] );	
		fprintf(file, "\n");
	}
	fclose(file);

	// Save energies
	file = fopen("energies.dat","w");
	for (i = 0; i < nbr_of_timesteps + 1; i++) {
		current_time = i * timestep;
		fprintf(file, "%.4f \t %e \t %e \t %e", current_time, E_pot[i], E_kin[i], E_tot[i] );	
		fprintf(file, "\n");
	}

	/*Save powerspectrum data in file */
  	file = fopen("powerspectrum.dat","w");
	for (i = 0; i < nbr_of_timesteps; i++)	{
		fprintf (file,"%e \t %e\n", freq[i], powspect[i]);
	}

	/* Free allocated memory */ 
	free(q_1); q_1 = NULL;
	free(q_2); q_2 = NULL;
	free(q_3); q_3 = NULL;
	free(E_kin); E_kin = NULL;
	free(E_pot); E_pot = NULL;
	free(E_tot); E_tot = NULL;
	return 0;    
}
