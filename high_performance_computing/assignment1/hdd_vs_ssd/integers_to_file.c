#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void benchmark_drive(char* name, FILE* file, size_t size, int data_size, int benchmark_iters){
  // Initialization
  struct timespec bench_start_time;
  struct timespec bench_stop_time;
  double  byte_write_time = 0;
  double byte_read_time = 0;
  double array_write_time = 0;
  double array_read_time = 0;
  int tmp_int = 0;
  int *data = (int*)malloc(sizeof(int)*size);
  int *read_data = (int*)malloc(sizeof(int)*size);
  for(int i=0; i<size; i++)
    data[i] = i;
    read_data[0] = 0;
  // Benchmark
  for(int i=0; i<benchmark_iters; i++){
    //******** Write bytes
    timespec_get(&bench_start_time, TIME_UTC);
    for(int j=0; j<size; j++){
      fwrite(&j, sizeof(int), 1, file);  
    }
    fflush(file);
    fseek(file, 0, SEEK_SET);
    timespec_get(&bench_stop_time, TIME_UTC);
    byte_write_time += difftime(bench_stop_time.tv_sec, bench_start_time.tv_sec)*1000 + (bench_stop_time.tv_nsec - bench_start_time.tv_nsec)/1000000;
    //******** Read bytes
    timespec_get(&bench_start_time, TIME_UTC);
    for(int j=0; j<size; j++){
      fread(&tmp_int, sizeof(int), 1, file);  
    }
    timespec_get(&bench_stop_time, TIME_UTC);
    byte_read_time += difftime(bench_stop_time.tv_sec, bench_start_time.tv_sec)*1000 + (bench_stop_time.tv_nsec - bench_start_time.tv_nsec)/1000000;

    //******** Write array
    timespec_get(&bench_start_time, TIME_UTC);
    fwrite(data, sizeof(int), size, file);
    fseek(file, 0, SEEK_SET);
    fflush(file);
    timespec_get(&bench_stop_time, TIME_UTC);
    array_write_time += difftime(bench_stop_time.tv_sec, bench_start_time.tv_sec)*1000 + (bench_stop_time.tv_nsec - bench_start_time.tv_nsec)/1000000;
    //******** Read array
    timespec_get(&bench_start_time, TIME_UTC);
    fread(read_data, sizeof(int), size, file);
    timespec_get(&bench_stop_time, TIME_UTC);
    array_read_time += difftime(bench_stop_time.tv_sec, bench_start_time.tv_sec)*1000 + (bench_stop_time.tv_nsec - bench_start_time.tv_nsec)/1000000;
    //printf("elem: %d ", read_data[rand()%size]); // print random element 
  }
  byte_write_time /= benchmark_iters;
  byte_read_time /= benchmark_iters;
  array_write_time /= benchmark_iters;
  array_read_time /= benchmark_iters;
  
  printf("\n\%s::\n", name);
  printf("\tByte write: %fmus/iter - %fMB/s\n", byte_write_time, data_size/(byte_write_time/1000));
  printf("\tByte read: %fmus/iter - %fMB/s\n", byte_read_time, data_size/(byte_read_time/1000));
  printf("\tByte write: %fmus/iter - %fMB/s\n", array_write_time, data_size/(array_write_time/1000));
  printf("\tByte read: %fmus/iter - %fMB/s\n", array_read_time, data_size/(array_read_time/1000));
  
  free(data); 
  free(read_data);
}


int
main(
    int argc,
    char* argv[]
){
  int benchmark_iters = 10;
  size_t size = 1048576;
  int data_size = size*4/1000000; // Data size in MB

  // Write to HDD
  FILE *hdd = fopen("./tmpHDD", "w+");
  benchmark_drive("HDD", hdd, size, data_size, benchmark_iters);
  fclose(hdd);
  // Write to SSD
  FILE *ssd = fopen("/run/mount/scratch/hpcuser272/tmpSSD", "w+");
  benchmark_drive("SSD", ssd, size, data_size, benchmark_iters);
  fclose(ssd);

}
