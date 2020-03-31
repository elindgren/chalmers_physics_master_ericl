/* 	
	Program for use on discrete Fourier transforms in E1 of Computaional Physics.
	Created by Martin Gren on 2015-10-13.
*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define PI 3.141592653589
#define n 250 /*number of timesteps*/

int main ()
{
	/* timestep dt */
	double dt = 0.1;

	/* Declare and set iteration parameter, data arrays, and frequency */
	int i;
	double data[n];
	double f = 1;
	
	/* declare file variables */
	FILE *file1;	

	/* calculate data points */
	for (i = 0; i < n; i++) {
		data[i] = cos(2*PI*f*dt*i + PI/2);
	}

	/*Save function values in file*/
	file1 = fopen("function_phi=90.dat","w");

	for (i = 0; i < n; i++) {
		fprintf (file1,"%e \t %e \n", i*dt, data[i]);
	}
	fclose(file1);

	return 0;
}
