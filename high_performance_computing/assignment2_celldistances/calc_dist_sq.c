#include<stdio.h>

long int calc_dist_sq(int* a, int* b){
  // Computes the distance squared between 3D points a and b.
  // Note that it is assumed that a and b are vectors of length 3.
  int d0 = a[0] - b[0];
  int d1 = a[1] - b[1];
  int d2 = a[2] - b[2];
  return d0*d0 + d1*d1 + d2*d2;
}
