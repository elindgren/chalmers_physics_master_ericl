#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<CL/cl.h>
#include"headerfile.h"


char * kernel_string = R"(
__kernel void prop_kern(
    __global const float *ht,
    __global float *htp,
    int ny,
    float e,
    float f
)
{
  int i = get_global_id(0) + 1; // Skip padded boundary of 0's
  int j = get_global_id(1) + 1;
  int idx = i*ny+j;
  htp[idx] = e * ht[idx] +\
      f * (ht[idx-ny] + ht[idx+ny] + ht[idx+1] + ht[idx-1]);
  
}
)";


void print_grid(float** h_old, int nx, int ny){
  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++)
      printf("%010.2f    ", h_old[i][j]);
    printf("\n");
  }
}

int main(int argc, char** argv){
  int nbr_iterations = 7;
  float d = 1/30.;
  parsing_cmdl_args(argc, argv, &nbr_iterations, &d);
  /************* Obtain initial conditions ***********/
  // Open file with initial conditions
  FILE *init_file;
  char filename[] = "./diffusion";
  init_file = fopen(filename, "r");
  if(init_file==NULL){
     printf("Init file could not be opened.");
     exit(1);
  }
  
  int nx, ny;
  read_header(init_file, &nx, &ny);
  printf("** Rows = %d Cols: = %d \n\n", nx,ny);
  // Allocate one extra space all the way around for boundaries
  nx += 2;
  ny += 2;
  // Allocate data structures	
  float* h_old_ent = (float*)malloc(sizeof(float)*nx*ny);
  float* h_new_ent = (float*)malloc(sizeof(float)*nx*ny);
  float** h_old = (float**)malloc(sizeof(float*)*nx);
  float** h_new = (float**)malloc(sizeof(float*)*nx);
  // Initialize contigously
  for(int x = 0; x<nx; x++){
    h_old[x] = h_old_ent + x * ny;
    h_new[x] = h_new_ent + x * ny;
    for(int y = 0; y<ny; y++){
      h_old[x][y] = 0;
      h_new[x][y] = 0;
    } 
  }
  read_data(init_file, h_old);
  // Close init file
  fclose(init_file);

  /*********** OpenCl configuration ************/
  printf("OpenCl\n");
  // Platform
  cl_int error;
  cl_platform_id platform_ids[2]; 
  cl_uint nmb_platforms;
  if(clGetPlatformIDs(2, platform_ids, &nmb_platforms) != CL_SUCCESS){
    printf("Cannot get OpenCl platform.\n");
    exit(1);
  }
  printf("** Number of platforms: %d -- IDs: %x and %x \n", nmb_platforms, platform_ids[0], platform_ids[1]);
  cl_platform_id platform_id = platform_ids[0]; // CPU in our case
  // Device
  cl_device_id device_id;
  cl_uint nmb_devices;
  if(clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_GPU, 1, &device_id, &nmb_devices) != CL_SUCCESS){
    printf("Cannot get OpenCl devices.\n");
    exit(1);
  }
  printf("** Number of devices: %d\n", nmb_devices);
  // Context
  cl_context context;
  cl_context_properties properties[] = {
    CL_CONTEXT_PLATFORM,
    (cl_context_properties) platform_id,
    0
  };
  context = clCreateContext(properties, 1, &device_id, NULL, NULL, NULL);
  // Command queue
  cl_command_queue command_queue;
  command_queue = clCreateCommandQueue(context, device_id, 0, &error);
  if(error != CL_SUCCESS){
    printf("Cannot create context or command queue.");
    exit(1);
  }
  // Create OpenCl program
  cl_program program;
  program = clCreateProgramWithSource(context, 1, 
      (const char**)&kernel_string, NULL, &error);
  clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
  cl_kernel kernel_odd, kernel_even;
  kernel_odd = clCreateKernel(program, "prop_kern", &error);
  kernel_even = clCreateKernel(program, "prop_kern", &error);
  if(error != CL_SUCCESS){
    printf("Failed to create kernel.");
    size_t log_size = 0;
    clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG,
        0, NULL, &log_size);
    char * log =  calloc(log_size, sizeof(char));
    clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG,
       log_size, log, NULL);
    printf("%s\n", log);
    free(log);
    exit(1);
  }
  // Setup buffers
  cl_mem buffer_h_old = clCreateBuffer(context, CL_MEM_READ_WRITE,
      nx*ny*sizeof(float), NULL, &error);
  cl_mem buffer_h_new = clCreateBuffer(context, CL_MEM_READ_WRITE,
      nx*ny*sizeof(float), NULL, &error);
  if(error != CL_SUCCESS){
    printf("Cannot setup buffers\n");
    exit(1);
  }
  // Write to buffer
  clEnqueueWriteBuffer(command_queue, buffer_h_old, CL_TRUE,
      0, nx*ny*sizeof(float), h_old_ent, 0, NULL, NULL);
  clEnqueueWriteBuffer(command_queue, buffer_h_new, CL_TRUE,
      0, nx*ny*sizeof(float), h_new_ent, 0, NULL, NULL);

  int nx_red = nx-2; // Actual matrix size - without boundary padding
  int ny_red = ny-2;
  float e = 1-d;
  float f = d*0.25;
  const size_t global[] = {nx_red,ny_red}; 
  // Prepare two kernels - one for even and one for odd iterations and fill their arguments beforehand
  clSetKernelArg(kernel_odd, 0, sizeof(cl_mem), &buffer_h_new);
  clSetKernelArg(kernel_odd, 1, sizeof(cl_mem), &buffer_h_old);
  clSetKernelArg(kernel_odd, 2, sizeof(int), &ny);
  clSetKernelArg(kernel_odd, 3, sizeof(int), &e);
  clSetKernelArg(kernel_odd, 4, sizeof(float), &f);

  clSetKernelArg(kernel_even, 0, sizeof(cl_mem), &buffer_h_old);
  clSetKernelArg(kernel_even, 1, sizeof(cl_mem), &buffer_h_new);
  clSetKernelArg(kernel_even, 2, sizeof(int), &ny);
  clSetKernelArg(kernel_even, 3, sizeof(int), &e);
  clSetKernelArg(kernel_even, 4, sizeof(float), &f);

  /*********** Perform calculation *************/
  printf("\nConfiguration complete, beginning calculations\n");
  for(int ix = 0; ix <nbr_iterations-1; ix+=2){
    clEnqueueNDRangeKernel(command_queue, kernel_even, 2, NULL,
    			   (const size_t*)&global, NULL, 0, NULL, NULL);
    clEnqueueNDRangeKernel(command_queue, kernel_odd, 2, NULL,
    			   (const size_t*)&global, NULL, 0, NULL, NULL);
  }
  /* Handle even or odd number of iterations */
  if(nbr_iterations%2 == 0){
    clEnqueueReadBuffer(command_queue, buffer_h_old, CL_TRUE,
      0, nx*ny*sizeof(float), h_old_ent, 0, NULL, NULL);
  }else{
    clEnqueueNDRangeKernel(command_queue, kernel_even, 2, NULL,
    			   (const size_t*)&global, NULL, 0, NULL, NULL);
    clEnqueueReadBuffer(command_queue, buffer_h_new, CL_TRUE,
      0, nx*ny*sizeof(float), h_old_ent, 0, NULL, NULL);
  }
  
  /*********** Calculate final averages ***********/ 
  float average = calc_average(h_old_ent, nx, ny);
  float average_distance = calc_distance_from_average(h_old_ent, nx, ny, average);
  printf("\nResults\n** After %d iterations we obtain:\n", nbr_iterations);
  printf("** Average temperature: %.2f\n** Average distance from average temperature: %.2f\n", average, average_distance);
  
  /* Cleanup */
  free(h_old_ent); free(h_new_ent);
  free(h_old); free(h_new);
  // Release OpenCl structures
  clReleaseContext(context);
  clReleaseCommandQueue(command_queue);
  clReleaseProgram(program);
  clReleaseKernel(kernel_odd);
  clReleaseKernel(kernel_even);
  clReleaseMemObject(buffer_h_old);
  clReleaseMemObject(buffer_h_new);
  return 0;
}
