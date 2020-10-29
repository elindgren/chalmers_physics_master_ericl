#include <math.h>
#include <stdio.h>


/*
 * Calculates the average temperature over the whole grid.
 * Note that nx and any runs over grid size -2 to skip boundaries.
 */
float calc_average(float* h, int nx, int ny){
  float average = 0.;
  for(int i=1; i<nx-1; i++){
    for(int j=1; j<ny-1; j++){
      int idx = i*ny + j;
      average += h[idx];
    }
  }
  return average/((nx-2)*(ny-2));
}

/* 
 * Calculate the average absolute distance from the average
 */

float calc_distance_from_average(float* h, int nx, int ny, float average){
  float average_distance = 0.;
  for(int i=1; i<nx-1; i++){
    for(int j=1; j<ny-1; j++){
      int idx = i*ny + j;
      average_distance += fabs( h[idx] - average );
    }
  }
  return average_distance/((nx-2)*(ny-2));
}
