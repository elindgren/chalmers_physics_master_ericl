#include<stdio.h>
#include<unistd.h>
#include<stdlib.h>

int parsing_cmdl_arg(int argc, char** argv){
  int threads = 1; // Default use 1

  // Use getopt to read number of threads from the command line
  char* cvalue = NULL;
  int c;
  opterr = 0;
  while((c = getopt(argc, argv, "t:")) != -1)
    switch (c){
      case 't':
        cvalue = optarg;
        break;
      default:
        printf("Invalid flag supplied to program. Exiting... \n");
        exit(1);
    }
  // Parse cvalue if not null
  if(cvalue != NULL){
    threads = atoi(cvalue);
    if(threads <= 0){
      printf("Invalid number of threads supplied. Exiting...\n");      
      exit(1);
    }  
  }

  return threads;
}

