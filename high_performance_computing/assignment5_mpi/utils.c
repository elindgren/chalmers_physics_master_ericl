#include <math.h>
#include<stdio.h>

/*
 * Calculates the average temperature over the whole grid.
 * Note that nx and any runs over grid size -2 to skip boundaries.
 */
float calc_average(float* h, int nx, int ny,int up){
  float average = 0;
  for(int i=up; i<ny-up; i++){
    for(int j=1; j<nx-1; j++){
      average += h[i*nx+j];
    }
  }
  return average/((nx-2)*(ny-2*up));
}

/* 
 * Calculate the average absolute distance from the average
 */
float calc_distance_from_average(float* h, int nx, int ny, float average, int up){
  float average_distance = 0;
  for(int i=up; i<ny-up; i++){
    for(int j=1; j<nx-1; j++){     
      average_distance += fabs(h[i*nx+j] - average);
    }
  }
  return average_distance/((nx-2)*(ny-2*up));
}
