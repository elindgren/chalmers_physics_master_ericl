#include<math.h>
#include<stdio.h>
/* 
 * Propagates the grid at time t to time t+1.
 *
 * Note that this implementation assumes 
 * the grid to be initialized with zeros 
 * on the border for the boundary conditions,
 * and thus an n+1 grid.
 * This is to simplify the code here and 
 * to avoid branch misses, at a slight penalty 
 * to memory.
 *
 * Params:
 *      ht: Grid at time t
 *      htp1: Grid at time t+1
 *      nx: Grid cols+2
 *      ny: Grid rows+2
 *      d: Diffusion coefficient
 * */

void propagateOneTimeStep(float* ht, float* htp, int nx, int ny, float d, float c){
  for(int i=1; i<ny-1; i++){
    for(int j=1; j<nx-1; j+=2){
      int idx = i*nx + j;
      int idx1 = idx+1;
      double tmp = ht[idx];
      double tmp1 = ht[idx1];
      htp[idx - nx] = c*tmp + d * (ht[idx-nx] + ht[idx+nx] + ht[idx-1] + tmp1);
      htp[idx1 - nx] = c*tmp1 + d * (ht[idx1-nx] + ht[idx1+nx] + tmp + ht[idx1+1]);
      //float term1 = c * ht[idx];
      //float term2 = d * (ht[idx-nx] + ht[idx+nx] + ht[idx-1] + ht[idx+1]);
      //htp[idx-nx] = term1 + term2;
    }
  }
}


/* Propagate one time step: Case when htp has border */
void propagateOneTimeStep1(float* ht, float* htp, int nx, int ny, float d, float c){
  for(int i=1; i<ny-1; i++){
    for(int j=1; j<nx-1; j+=2){
      int idx = i*nx + j;
      int idx1 = idx+1;
      double tmp = ht[idx];
      double tmp1 = ht[idx1];
      htp[idx] = c*tmp + d * (ht[idx-nx] + ht[idx+nx] + ht[idx-1] + tmp1);
      htp[idx1] = c*tmp1 + d * (ht[idx1-nx] + ht[idx1+nx] + tmp + ht[idx1+1]);
    }
  }
}

/* Propagate one time step: Case when htp has no border, and thus the rows are offset */
void propagateOneTimeStep2(float* ht, float* htp, int nx, int ny, float d, float c, int up){
  for(int i=0; i<ny; i++){
    for(int j=1; j<nx-1; j+=2){
      int idx = (i+up)*nx + j;
      int idx1 = idx+1;
      double tmp = ht[idx];
      double tmp1 = ht[idx1];
      htp[i*nx+j] = c*tmp + d * (ht[idx-nx] + ht[idx+nx] + ht[idx-1] + tmp1);
      htp[i*nx+j+1] = c*tmp1 + d * (ht[idx1-nx] + ht[idx1+nx] + tmp + ht[idx1+1]);
    }
  }
}

