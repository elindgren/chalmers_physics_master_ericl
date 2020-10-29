#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
#include"headerfile.h"

void print_grid(float* h_old, int nx, int ny){
  for(int i=0; i<ny; i++){
    for(int j=0; j<nx; j++)
      printf("%010.2f    ", h_old[i*nx + j]);
    printf("\n");
  }
}

static inline
void propagateOneTimeStep1(float* ht, float* htp, int nx, int ny, float d, float c, int os){
  /* Propagates the grid ht into htp. 
     os is offset that specifies which row to start and end the calculation at.
     We dont want to calculate unnecessary rows which are send for padding */ 
  for(int i=1+os; i<ny-1-os; i++){
    for(int j=1; j<nx-1; j+=2){
      // this loop is unrolled for performance only works for even gridsizes
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
  /* Propagated the grid ht into htp.
     htp is here assumed to be the sub array returned by gatherv 
     (different size then in propagateOneTimeStep1)
  */
  for(int i=0; i<ny; i++){
    for(int j=1; j<nx-1; j+=2){
      // this loop is unrolled for performance only works for even gridsizes
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
  /*** MPI setup ***/
  MPI_Init(&argc,&argv);
  int nbr_mpi_proc, mpi_rank;
  MPI_Comm_size(MPI_COMM_WORLD, &nbr_mpi_proc);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);  
  int scatter_root = 0;
  int up = 25; // Amount of extra rows that need to be included for unrolling iteration loop.
               // All data sent by scatterv  need to be padded with up rows so the return rows are
               // correct after a number of timesteps without syncronization.
               // This means that we can perform up timesteps without synchronization (scatterv/gatherv).
  
  /*** setup grid, sizes and parameters for calculation ***/
  int nbr_iterations = 7;
  float d = 1/30.;
  float* h_old;
  float* h_new;
  int nx, ny;          // number of rows and columes with added padding of zeros
  int nx_red, ny_red;  // number of rows and colums we acctually want to calculate
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

  }
  /*** setup scatter and broadcast data to all processes ***/
  MPI_Bcast(&nx, 1, MPI_INT, scatter_root, MPI_COMM_WORLD);
  MPI_Bcast(&ny, 1, MPI_INT, scatter_root, MPI_COMM_WORLD);
  MPI_Bcast(&ny_red, 1, MPI_INT, scatter_root, MPI_COMM_WORLD);
  MPI_Bcast(&nx_red, 1, MPI_INT, scatter_root, MPI_COMM_WORLD);
  MPI_Bcast(&nbr_iterations, 1, MPI_INT, scatter_root, MPI_COMM_WORLD);
  MPI_Bcast(&d, 1, MPI_FLOAT, scatter_root, MPI_COMM_WORLD);
  int nbr_row = (ny_red - 1)/nbr_mpi_proc + 1;             // nbr of rows each prossec will compute
  int* lens = (int*) malloc(sizeof(int)*nbr_mpi_proc);     // lens for scatered data
  int* poss = (int*) malloc(sizeof(int)*nbr_mpi_proc);     // position of scattered data
  int* lens_ret = (int*) malloc(sizeof(int)*nbr_mpi_proc); // lens of gathered data
  int* poss_ret = (int*) malloc(sizeof(int)*nbr_mpi_proc); // postion of gathered data
  for(int ix = 0; ix<nbr_mpi_proc; ix++){
    // calculate lens and position for date to send and recieve for all processes
    lens[ix] = ix < nbr_mpi_proc-1 ? (nbr_row + 2*up)*nx : (ny_red - nbr_row*ix+2*up)*nx;
    poss[ix] = nbr_row*nx*ix;
    lens_ret[ix] = lens[ix] - 2*up*nx;
    poss_ret[ix] = poss[ix] + up*nx;
  }

  // allocate sub arrays of the grid on all processes to which data is scattered
  float* h_old_sub = (float*) malloc(sizeof(float)*lens[mpi_rank]); 
  float* h_tmp_sub = (float*) malloc(sizeof(float)*lens[mpi_rank]); 
  float* h_new_sub = (float*) malloc(sizeof(float)*lens_ret[mpi_rank]); // sub array wich are returned by gatherv
  for(int i = 0; i<lens[mpi_rank]; i++){
    h_old_sub[i] = 0;
    h_tmp_sub[i] = 0;
  }
  for(int i = 0; i<lens_ret[mpi_rank];i++){
    h_new_sub[i] = 0;
  }

  /*** perform iteration loop ***/
  float c = 1 - d;       //speedup propagation 
  float d_red = d*0.25;  //speedup propagation
  for(int ix = 0; ix<nbr_iterations;ix+=up){
    // unroll iteration to reduce syncronization when calling scatterv / gatherv
    MPI_Scatterv(h_old, lens, poss, MPI_FLOAT, h_old_sub, lens[mpi_rank], MPI_FLOAT,0,MPI_COMM_WORLD);
    for(int i = 0; i<up-2; i+=2){
      // only works for odd up
      propagateOneTimeStep1(h_old_sub, h_tmp_sub, nx, lens[mpi_rank]/nx, d_red, c, i);
      propagateOneTimeStep1(h_tmp_sub, h_old_sub, nx, lens[mpi_rank]/nx, d_red, c, i+1);
    }
    propagateOneTimeStep2(h_old_sub, h_new_sub, nx, lens_ret[mpi_rank]/nx, d_red, c, up);
    MPI_Gatherv(h_new_sub, lens_ret[mpi_rank], MPI_FLOAT, h_old, lens_ret, poss_ret, MPI_FLOAT, 0, MPI_COMM_WORLD);
  }
  
  /**** Print resulting average and clean up alocations ****/
  if(mpi_rank==0){ 
    float average = calc_average(h_old, nx, ny, up);
    float average_distance = calc_distance_from_average(h_old, nx, ny, average,up);
    printf("\nResults\n** After %d iterations we obtain:\n", nbr_iterations);
    printf("** Average temperature: %.2f\n** Average distance from average temperature: %.2f\n",
    	   average, average_distance);
    free(h_old); free(h_new);
  }
  free(lens_ret); free(poss_ret);
  free(h_old_sub); free(h_new_sub);
  free(h_tmp_sub);
  free(lens); free(poss);
  MPI_Finalize();
  return 0;
}




