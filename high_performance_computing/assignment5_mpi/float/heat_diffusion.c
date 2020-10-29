#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
#include"headerfile.h"
#define bl_size 100

void print_grid(float* h_old, int nx, int ny){
  for(int i=0; i<ny; i++){
    for(int j=0; j<nx; j++)
      printf("%010.2f    ", h_old[i*nx + j]);
    printf("\n");
  }
}
/*
static inline
void propagateOneTimeStep1(float* ht, float* htp, int nx, int ny, float d, float c){
  for(int bli = 0; bli<(ny-2)/bl_size; bli+=bl_size){
    for(int blj = 0; blj<(nx-2)/bl_size; blj+=bl_size){
      for(int i=bli+1; i<bli+bl_size+1; i++){
        for(int j=blj+1; j<blj+bl_size+1; j+=2){
          //printf("bli blj i j %d %d %d %d \n", bli, blj, i, j);
          int idx = i*nx + j;
          int idx1 = idx+1;
          float tmp = ht[idx];
          float tmp1 = ht[idx1];
          htp[idx] = c*tmp + d * (ht[idx-nx] + ht[idx+nx] + ht[idx-1] + tmp1);
          htp[idx1] = c*tmp1 + d * (ht[idx1-nx] + ht[idx1+nx] + tmp + ht[idx1+1]);
        }
      }
    }
  }
}

static inline
void propagateOneTimeStep2(float* ht, float* htp, int nx, int ny, float d, float c,int up){
  for(int bli = 0; bli<ny; bli+=bl_size){
    for(int blj = 0; blj<(nx-2); blj+=bl_size){
      for(int i=bli; i<bli+bl_size; i++){
        for(int j=blj+1; j<blj+bl_size+1; j+=2){
          //printf("bli blj i j %d %d %d %d \n", bli, blj, i, j);
          int idx = (i+up)*nx + j;
          int idx1 = idx+1;
          int idxret = i*nx+j;
          float tmp = ht[idx];
          float tmp1 = ht[idx1];
          htp[idxret] = c*tmp + d * (ht[idx-nx] + ht[idx+nx] + ht[idx-1] + tmp1);
          htp[idxret+1] = c*tmp1 + d * (ht[idx1-nx] + ht[idx1+nx] + tmp + ht[idx1+1]);
        }
      }
    }
  }
}
*/

static inline
void propagateOneTimeStep1(float* ht, float* htp, int nx, int ny, float d, float c){
  for(int i=1; i<ny-1; i++){
    for(int j=1; j<nx-1; j+=2){
      int idx = i*nx + j;
      int idx1 = idx+1;
      float tmp = ht[idx];
      float tmp1 = ht[idx1];
      htp[idx] = c*tmp + d * (ht[idx-nx] + ht[idx+nx] + ht[idx-1] + tmp1);
      htp[idx1] = c*tmp1 + d * (ht[idx1-nx] + ht[idx1+nx] + tmp + ht[idx1+1]);
    }
  }
}

static inline
void propagateOneTimeStep2(float* ht, float* htp, int nx, int ny, float d, float c,int up){
  for(int i=0; i<ny; i++){
    for(int j=1; j<nx-1; j+=2){
      int idx = (i+up)*nx + j;
      int idx1 = idx+1;
      int idxret = i*nx+j;
      float tmp = ht[idx];
      float tmp1 = ht[idx1];
      htp[idxret] = c*tmp + d * (ht[idx-nx] + ht[idx+nx] + ht[idx-1] + tmp1);
      htp[idxret+1] = c*tmp1 + d * (ht[idx1-nx] + ht[idx1+nx] + tmp + ht[idx1+1]);
    }
  }
}


int main(int argc, char** argv){
  MPI_Init(&argc,&argv);
  int nbr_mpi_proc, mpi_rank;
  MPI_Comm_size(MPI_COMM_WORLD, &nbr_mpi_proc);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);  
  int scatter_root = 0;
  int up = 25; // Amount of extra rows that need to be included for unrolling

  int nbr_iterations = 7;
  float d = 1/30.;
  float* h_old;
  float* h_new;
  int nx, ny;
  int nx_red, ny_red;
  if(mpi_rank == 0){
    parsing_cmdl_args(argc,argv,&nbr_iterations,&d);
    FILE *init_file;
    char filename[] = "./diffusion";
    init_file = fopen(filename, "r");
    if(init_file==NULL){
      printf("Init file could not be opened.");
      exit(1);
    }
    read_header(init_file, &ny_red ,&nx_red);
    ny = ny_red + 2*up;
    nx = nx_red + 2;
    
    printf("** Rows = %d Cols: = %d \n\n", ny,nx);
    h_old = (float*) malloc(sizeof(float)*nx*ny);
    h_new = (float*) malloc(sizeof(float)*nx*ny);
    for(int i = 0; i<ny; i++){
      for(int j = 0; j<nx ;j++){
	h_old[i*nx + j] = 0;
	h_new[i*nx + j] = 0;
      }
    }
    read_data(init_file, nx, up,  h_old);
    fclose(init_file);
    //print_grid(h_old,nx,ny);
  }
  // setup scatter on all processes
  
  MPI_Bcast(&nx, 1, MPI_INT, scatter_root, MPI_COMM_WORLD);
  MPI_Bcast(&ny, 1, MPI_INT, scatter_root, MPI_COMM_WORLD);
  MPI_Bcast(&ny_red, 1, MPI_INT, scatter_root, MPI_COMM_WORLD);
  MPI_Bcast(&nx_red, 1, MPI_INT, scatter_root, MPI_COMM_WORLD);
  MPI_Bcast(&nbr_iterations, 1, MPI_INT, scatter_root, MPI_COMM_WORLD);
  MPI_Bcast(&d, 1, MPI_DOUBLE, scatter_root, MPI_COMM_WORLD);
  int nbr_row = (ny_red - 1)/nbr_mpi_proc + 1;
  int* lens = (int*) malloc(sizeof(int)*nbr_mpi_proc);
  int* poss = (int*) malloc(sizeof(int)*nbr_mpi_proc);
  int* lens_ret = (int*) malloc(sizeof(int)*nbr_mpi_proc);
  int* poss_ret = (int*) malloc(sizeof(int)*nbr_mpi_proc);
  for(int ix = 0; ix<nbr_mpi_proc; ix++){
    lens[ix] = ix < nbr_mpi_proc-1 ? (nbr_row + 2*up)*nx : (ny_red - nbr_row*ix+2*up)*nx;
    poss[ix] = nbr_row*nx*ix;
    lens_ret[ix] = lens[ix] - 2*up*nx;
    poss_ret[ix] = poss[ix] + up*nx;
    
    if(mpi_rank==0){
      printf("lens lend_ret  och poss poss_ret %d %d %d %d \n",lens[ix], lens_ret[ix],poss[ix],poss_ret[ix]);
      //printf("lens[%d] = %d",ix, lens[ix]);
    }
    
  }
  //printf("lens: %d nx: %d  - lens/nx: %d \n", lens[mpi_rank], nx, lens[mpi_rank]/nx);
  float* h_old_sub = (float*) malloc(sizeof(float)*lens[mpi_rank]);
  float* h_tmp_sub = (float*) malloc(sizeof(float)*lens[mpi_rank]);
  float* h_new_sub = (float*) malloc(sizeof(float)*lens_ret[mpi_rank]);
  for(int i = 0; i<lens[mpi_rank]; i++){
    h_old_sub[i] = 0;
    h_tmp_sub[i] = 0;
  }
  for(int i = 0; i<lens_ret[mpi_rank];i++){
    h_new_sub[i] = 0;
  }
  float c = 1 - d;
  float d_red = d*0.25;
  //MPI_Scatterv(h_old, lens, poss, MPI_FLOAT, h_old_sub, lens[mpi_rank], MPI_FLOAT,0,MPI_COMM_WORLD);
  for(int ix = 0; ix<nbr_iterations;ix+=up){
    MPI_Scatterv(h_old, lens, poss, MPI_FLOAT, h_old_sub, lens[mpi_rank], MPI_FLOAT,0,MPI_COMM_WORLD);
    propagateOneTimeStep1(h_old_sub, h_tmp_sub, nx, lens[mpi_rank]/nx, d_red, c);
    propagateOneTimeStep1(h_tmp_sub, h_old_sub, nx, lens[mpi_rank]/nx, d_red, c);
    propagateOneTimeStep1(h_old_sub, h_tmp_sub, nx, lens[mpi_rank]/nx, d_red, c);
    propagateOneTimeStep1(h_tmp_sub, h_old_sub, nx, lens[mpi_rank]/nx, d_red, c);
    propagateOneTimeStep1(h_old_sub, h_tmp_sub, nx, lens[mpi_rank]/nx, d_red, c);
    propagateOneTimeStep1(h_tmp_sub, h_old_sub, nx, lens[mpi_rank]/nx, d_red, c);
    propagateOneTimeStep1(h_old_sub, h_tmp_sub, nx, lens[mpi_rank]/nx, d_red, c);
    propagateOneTimeStep1(h_tmp_sub, h_old_sub, nx, lens[mpi_rank]/nx, d_red, c);
    propagateOneTimeStep1(h_old_sub, h_tmp_sub, nx, lens[mpi_rank]/nx, d_red, c);
    propagateOneTimeStep1(h_tmp_sub, h_old_sub, nx, lens[mpi_rank]/nx, d_red, c);
    propagateOneTimeStep1(h_old_sub, h_tmp_sub, nx, lens[mpi_rank]/nx, d_red, c);
    propagateOneTimeStep1(h_tmp_sub, h_old_sub, nx, lens[mpi_rank]/nx, d_red, c);
    propagateOneTimeStep1(h_old_sub, h_tmp_sub, nx, lens[mpi_rank]/nx, d_red, c);
    propagateOneTimeStep1(h_tmp_sub, h_old_sub, nx, lens[mpi_rank]/nx, d_red, c);
    propagateOneTimeStep1(h_old_sub, h_tmp_sub, nx, lens[mpi_rank]/nx, d_red, c);
    propagateOneTimeStep1(h_tmp_sub, h_old_sub, nx, lens[mpi_rank]/nx, d_red, c);
    propagateOneTimeStep1(h_old_sub, h_tmp_sub, nx, lens[mpi_rank]/nx, d_red, c);
    propagateOneTimeStep1(h_tmp_sub, h_old_sub, nx, lens[mpi_rank]/nx, d_red, c);
    propagateOneTimeStep1(h_old_sub, h_tmp_sub, nx, lens[mpi_rank]/nx, d_red, c);
    propagateOneTimeStep1(h_tmp_sub, h_old_sub, nx, lens[mpi_rank]/nx, d_red, c);
    propagateOneTimeStep1(h_old_sub, h_tmp_sub, nx, lens[mpi_rank]/nx, d_red, c);
    propagateOneTimeStep1(h_tmp_sub, h_old_sub, nx, lens[mpi_rank]/nx, d_red, c);
    propagateOneTimeStep1(h_old_sub, h_tmp_sub, nx, lens[mpi_rank]/nx, d_red, c);
    propagateOneTimeStep1(h_tmp_sub, h_old_sub, nx, lens[mpi_rank]/nx, d_red, c);
    
    propagateOneTimeStep2(h_old_sub, h_new_sub, nx, lens_ret[mpi_rank]/nx, d_red, c, up);
    //propagateOneTimeStep1(h_old_sub, h_new_sub, nx, lens[mpi_rank]/nx, d_red, c);

    /*
    for(int i=1; i<lens[mpi_rank]/nx-1; i++){ 
      for(int j=1; j<nx-1; j++){ 
	int idx = i*nx + j;
	float term1 = c * h_old_sub[idx];
	float term2 = d_red * (h_old_sub[idx-nx] + h_old_sub[idx+nx] + h_old_sub[idx-1] + h_old_sub[idx+1]);
	h_new_sub[idx-nx] = term1 + term2;
      }
    }
    */
    MPI_Gatherv(h_new_sub, lens_ret[mpi_rank], MPI_FLOAT, h_old, lens_ret, poss_ret, MPI_FLOAT, 0, MPI_COMM_WORLD);
  }
  //MPI_Gatherv(h_new_sub, lens_ret[mpi_rank], MPI_FLOAT, h_old, lens_ret, poss_ret, MPI_FLOAT, 0, MPI_COMM_WORLD);

  if(mpi_rank==0){ 
    float average = calc_average(h_old, nx, ny, up);
    float average_distance = calc_distance_from_average(h_old, nx, ny, average,up);
    printf("\nResults\n** After %d iterations we obtain:\n", nbr_iterations);
    printf("** Average temperature: %.2f\n** Average distance from average temperature: %.2f\n",
    	   average, average_distance);
    //print_grid(h_old, nx, ny);
    
    free(h_old); free(h_new);
  }
  free(lens_ret); free(poss_ret);
  free(h_old_sub); free(h_new_sub);
  free(h_tmp_sub);
  free(lens); free(poss);
  MPI_Finalize();
  return 0;
}




