#include<stdio.h>
#include<unistd.h>
#include<stdlib.h>


void parsing_cmdl_args(int argc, char** argv, int *nbr_iteration, float *d, int *local){
  // Use getopt to read number of iterations and the diffusion constant 
  char* nvalue = NULL;
  char* dvalue = NULL;
  char* lvalue = NULL;
  int c;
  opterr = 0;
  int optind = 0;
  while((c = getopt(argc, argv, "d:n:l:")) != -1)
    switch (c){
      case 'd':
        dvalue = optarg;
        break;
      case 'n':
        nvalue = optarg;
        break;
      case 'l':
        lvalue = optarg;
        break;
      default:
        printf("Invalid flag supplied to program. Exiting... \n");
        exit(1);
    }
  // Parse cvalue if not null
  printf("Parameters\n");
  if(nvalue != NULL && dvalue != NULL){
    int ntmp = atoi(nvalue);
    float dtmp = atof(dvalue);
    if(ntmp <= 0 || dtmp <=0){
      printf("Invalid parameters. Exiting...\n");      
      exit(1);
    } else {
      *nbr_iteration = ntmp;
      *d = dtmp;
    }
  } else {
    printf("** Using default values for n and d\n");
  }

  if(lvalue != NULL){
    int ltmp = atoi(lvalue);
    if(ltmp <= 0){
      printf("Invalid local size. Exiting...\n");
      exit(1);
    } else{
      *local = ltmp;
    }
  }
  printf("** Number of iterations: %d\n", *nbr_iteration);
  printf("** Diffusion constant: %f\n", *d);
  printf("** OpenCl Local size: (%d, %d) \n", *local, *local);
}

