#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<math.h>
#include<omp.h>
#include "header_file.h"


void calc_dist_bl(int** block1, int bl_size, int* dist_entries, int threads){
  /* Calculate the distances between all points in a block, and count them once. */
  #pragma omp parallel for reduction(+:dist_entries[:3465])
  for(int i=0; i<bl_size; i++){
    for(int j=i+1; j<bl_size; j++){
      int dist_sq = calc_dist_sq(block1[i], block1[j]);
      int dist = (int)(1e-1 * sqrt(dist_sq)); // Truncating
      dist_entries[dist] += 1;
    }
  }  
}

void calc_dist_bl_pair(int** block1, int** block2, int bl_size1, int bl_size2, int *dist_entries){
  /* Calculate the pairwise distances between points in different blocks. */
  #pragma omp parallel for reduction(+:dist_entries[:3465])
  for(int i=0; i<bl_size1; i++){
    for(int j=0; j<bl_size2; j++){
      int dist_sq = calc_dist_sq(block1[i], block2[j]);
      int dist = (int) (1e-1 * sqrt(dist_sq));
      dist_entries[dist] += 1;
    }
  }  
}
int main(int argc, char** argv){
  /* 
   * This program calculates the distances between the points in the 
   * file cells. 
   *
   * To fulfill the memory requirement of this program, we settled for an 
   * approach where we split the data into blocks and compute all 
   * distances between all points in and between all blocks. 
   *
   * The program uses at most two blocks. To minimize reading from the file 
   * the program tries to fit as much data as possible into one block. 
   * */

  /* Config and setup */
  int threads = parsing_cmdl_arg(argc, argv);
  char data[] = "./cells";
  FILE* file;
  
  omp_set_num_threads(threads);

  int N_points = 3465; // Maximum distance is (10--10)^2*3 = 1200 units. 
                       // The sqrt of this is 34.64, to two decimal places.
                       // Converting this to int (3464), we thus need to fit
                       // 3465 points.
  int *dist_entries = (int*)malloc(sizeof(int)*N_points); 
  for(int i=0; i<N_points; i++){
    dist_entries[i] = 0;
  }

  // Open file
  file = fopen(data, "r");
  if(file == NULL){
    printf("The file could not be opened\n");
    exit(1);
  }
  
  /* Setup blocks */

  // We have four major allocations in our program. These are dependent on the block size. To make sure that fit into 1 GiBi, we need to calculate the maximum block size. 
  // This was calculated as:
  //    *block_entries_1: 4*3*bl_size
  //    *block_entries_2: -||-
  //    **block1: 8*bl_size
  //    **block2: -||-
  //    + much smaller allocations, like *dist_entries
  //    -----------------------
  //    = 40*bl_size. 
  // This should be equal to at most 1 GiBi = 2^30
  //    bl_size < 2^30/40 \approx 2^25
  // Setting bl_size = 2^24 is thus on the safe side, giving us room for other allocations.
  int max_n_lines = 1<<24; // Maximum number of lines that will fit in 1 GiB 
  
  // Calculate block size and number of blocks - maximum 2^32 lines
  fseek(file, 0, SEEK_END);
  int char_per_line = 24;
  unsigned int nbr_of_lines_in_file = ftell(file)/char_per_line;
  int nbr_bl = (int)((float)nbr_of_lines_in_file/max_n_lines + 1-1e-32); // Calculate number of lines, but avoid the limiting case where nbr_lines/max = exact by not adding 1, but slightly less.
  int bl_size = nbr_bl > 1 ? max_n_lines : nbr_of_lines_in_file; 

  // Allocate memory for blocks
  int *block_entries1 = (int*)malloc(sizeof(int)*(bl_size*3+1));
  int *block_entries2 = (int*)malloc(sizeof(int)*(bl_size*3+1));
  int **block1 = (int**)malloc(sizeof(int*)*bl_size);
  int **block2 =  (int**)malloc(sizeof(int*)*bl_size);


  /* Calculation */
  for(int i=0; i<bl_size; i++){
    block1[i] = block_entries1 + i*3;
    block2[i] = block_entries2 + i*3;  
  }
  for(size_t bl = 0; bl<nbr_bl; bl++){
    int bl_size1 = read_data(block1, bl*bl_size, bl_size, file);
    calc_dist_bl(block1, bl_size1, dist_entries, threads);
 
    for(size_t bl2 = bl+1; bl2<nbr_bl; bl2++ ){
      int bl_size2 = read_data(block2, bl2*bl_size, bl_size, file);
      calc_dist_bl_pair(block1, block2, bl_size1, bl_size2, dist_entries);
    }
  }

  /* Print results */
  long int sum = 0; // DEBUG
  for(int i=0; i<N_points; i+=1){
    if(dist_entries[i] > 0){
      printf("%05.2f %d\n", (float)i/100, dist_entries[i]);
      sum += dist_entries[i]; // DEBUG
    }
  }
  long int checksum = ((long int)nbr_of_lines_in_file-1)*(long int)nbr_of_lines_in_file/2;
  if(sum != checksum){
    printf("Sum of nbr of distances is invalid.\n The sum is: %lu \n The sum should be: %lu \nExiting...\n", sum, checksum);
    exit(1);
  }
  
  /* Cleanup */
  fclose(file);
  free(block1); free(block_entries1);
  free(block2); free(block_entries2);
  free(dist_entries);
  return 0;
}
