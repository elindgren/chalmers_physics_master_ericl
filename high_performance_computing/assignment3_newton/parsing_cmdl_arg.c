#include<stdio.h>
#include<unistd.h>
#include<stdlib.h>
#include "argstruct.h"

struct Args parsing_cmdl_arg(int argc, char** argv){
  int threads = 1; // Default use 1
  int lines = 1000; // Default use 1
  int d = 2;
  // Use getopt to read number of threads and lines from the command line
  char* tvalue = NULL;
  char* lvalue = NULL;
  char* dvalue = NULL;
  int c;
  opterr = 0;
  int optind = 0;
  while(optind < argc){
    if((c = getopt(argc, argv, "t:l:")) != -1)
      switch (c){
        case 't':
          tvalue = optarg;
          break;
        case 'l':
          lvalue = optarg;
          break;
        default:
         printf("Invalid flag supplied to program. Exiting... \n");
         exit(1);
    }else{
      if(optind = argc-1){
        dvalue = argv[optind];
      }
      optind++;
    }

  }
  // Parse cvalue if not null
  if(tvalue != NULL && lvalue != NULL && dvalue != NULL){
    threads = atoi(tvalue);
    lines = atoi(lvalue);
    d = atoi(dvalue);
    if(threads <= 0 || lines <=0 || d<=0 || d>10){
      printf("Invalid number of threads supplied. Exiting...\n");      
      exit(1);
    }  
  } else {
    printf("Using default values\n");
  }

  printf("Number of threads: %d\n", threads);
  printf("Number of lines: %d\n", lines);
  printf("Degree of polynomial: %d\n", d);
  struct Args a = {threads, lines, d};
  return a;
}

