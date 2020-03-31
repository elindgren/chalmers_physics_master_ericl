/*
 E2.c
 
 Created by Anders Lindman on 2014-11-04.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define PI 3.141592653589

/* Main program */
int main()
{
    
    /* It is useful to construct the transformation matrix outside the main loop */
    int i,j;
    int nbr_of_particles;
    double factor;
    double trans_matrix[nbr_of_particles][nbr_of_particles];
    
    factor = 1 / ((double) nbr_of_particles + 1);
    for (i=0; i < nbr_of_particles; i++) {
        for (j=0; j < nbr_of_particles; j++) {
            trans_matrix[i][j] = sqrt(2 * factor) * sin((j + 1) * (i + 1) * PI * factor);
        }
    }
    
    /* Transformation to normal modes Q from displacements q.  */
    double sum;
    for (i = 0; i < nbr_of_particles; i++){
        sum = 0;
        for (j = 0; j < nbr_of_particles; j++){
            sum += q[j] * trans_matrix[i][j];
        }
        Q[i] = sum;
    }
    
    
}
