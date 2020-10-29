#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "headerfile.h"
#include "argstruct.h"
#include <threads.h>
#define PI 3.14159265358979323846

typedef struct {
  int val; 
  char pad[60]; //cacheline - sizeof(int)
} int_padded;

typedef struct{
  char** attractors;
  char** convergences;
  int ib;
  const int *d;
  const int *conv_cap;
  const double *h;
  const int *istep;
  int tx;
  const int *lines;
  const double complex* roots;
  int_padded *status;
  cnd_t *cnd;
  mtx_t *mtx;
}thrd_info_t;

typedef struct{
  char** attractors;
  char** convergences;
  const char** grayscale;
  const int *lines;
  const int *nthrds;
  const int *d;
  const int *conv_cap;
  FILE* attr_file;
  FILE* conv_file;
  mtx_t *mtx;
  cnd_t *cnd;
  int_padded *status;
}write_info_t;

/*
 * Thread for performing Newton algorithm calculations.
 * Each threads performs calculations on one row at a time,
 * for all lines that has been allocated to that thread.
 */
int main_thrd(void* args){
  thrd_info_t* thrd_info = (thrd_info_t*) args;
  char** attractors = thrd_info->attractors;
  char** convergences = thrd_info->convergences;
  int ib = thrd_info->ib;
  const int lines = *thrd_info->lines;
  const int d = *thrd_info->d;
  const int conv_cap = *thrd_info->conv_cap;
  const int istep = *thrd_info->istep;
  const double h = *thrd_info->h;
  const double complex* roots = thrd_info->roots;
  int_padded *status = thrd_info->status;
  int tx = thrd_info->tx;
  mtx_t *mtx = thrd_info->mtx;
  cnd_t *cnd = thrd_info->cnd;
  double complex x0;

  for(int ix = ib; ix<lines; ix+=istep){
    char* attractor = (char*) malloc(sizeof(char)* lines);
    char* convergence = (char*) malloc(sizeof(char)* lines);
    for(int jx = 0; jx <lines; jx++){
      // Perform calculation for every pixel on this row
      x0 = -2 + h*jx + I*(h*ix-2);
      newton_iteration(&attractor[jx], &convergence[jx], x0, d, conv_cap, roots);
    } 
    mtx_lock(mtx);
    attractors[ix] = attractor;
    convergences[ix] = convergence;
    status[tx].val = ix+istep;
    mtx_unlock(mtx);
    cnd_signal(cnd);
  }
}

/* 
 * Thread for writing to file. Utilizes spurious waking as in the example
 * code on the course page (notice the while(1) loop; it relies on being
 * awoken from cnd_wait by the system). 
 */
int main_thrd_write(void* args){
    write_info_t* thrd_info = (write_info_t*) args;
    char** attractors = thrd_info->attractors;
    char** convergences = thrd_info->convergences;
    const char** grayscale = thrd_info->grayscale;
    const int lines = *thrd_info->lines;
    const int nthrds = *thrd_info->nthrds;
    const int d = *thrd_info->d;
    const int conv_cap = *thrd_info->conv_cap;
    FILE* attr_file = thrd_info->attr_file;
    FILE* conv_file = thrd_info->conv_file;
    int_padded *status = thrd_info->status;
    mtx_t *mtx = thrd_info->mtx;
    cnd_t *cnd = thrd_info->cnd;
    
    for(int ix = 0, ibnd; ix<lines;){
      mtx_lock(mtx);
      while(1){
        ibnd = lines;
        for(int tx = 0; tx<nthrds; tx++){
          if(status[tx].val < ibnd)
            ibnd = status[tx].val; 
        } 
        if(ibnd<=ix)
          cnd_wait(cnd,mtx);
        else{
          mtx_unlock(mtx);
          break;
        }
      }
      // Write the rows that are finished to file
      for(;ix<ibnd;ix++){
        write_row_to_file(attractors[ix], convergences[ix], attr_file, conv_file, lines, d, conv_cap, grayscale);
        // Cleanup
        free(attractors[ix]); free(convergences[ix]);
      }
    }
}

int main(int argc, char* argv){
  /* Fetch parameters */
  struct Args args = parsing_cmdl_arg(argc, argv);
  int threads = args.threads;
  int lines = args.lines;
  int d = args.d;
  /* General parameters */
  int conv_cap = 100;
  double delta_h = 4.0/(lines - 1);
  double complex roots[d];
  for(int ix = 0; ix < d; ix++){
    roots[ix] = cos(2*PI*ix/d) + I * sin(2*PI*ix/d);
  }

  /* Alocate structures */
  // The output from each thread will be a row - save each row in a global array
  char** attractors = (char**) malloc(sizeof(char*)*lines);
  char** convergences = (char**) malloc(sizeof(char*)*lines);

  /* Initialize grayscale values */
  char* grayscale_ent = (char*) malloc(sizeof(char)*(conv_cap+1)*12);
  char** grayscale = (char**) malloc(sizeof(char*)*(conv_cap+1));

  for(int li=0; li<conv_cap+1; li++){
    grayscale[li] =  grayscale_ent + li*12;
    sprintf(grayscale[li], "%03d %03d %03d ", li,li,li);
  }
  
  /* Initialize files */
  FILE* attr_file = initialize_attractor_file(lines, d);
  FILE* conv_file = initialize_convergence_file(lines, d, conv_cap);

  /* Threading setup */
  thrd_t thrds[threads];
  thrd_info_t thrd_info[threads];
  
  thrd_t thrd_write;
  write_info_t write_info;

  mtx_t mtx; 
  mtx_init(&mtx, mtx_plain);

  cnd_t cnd;
  cnd_init(&cnd);

  int_padded status[threads];
  
  /* Setup and start all threads */
  for(int tx = 0; tx<threads; tx++){
    thrd_info[tx].attractors = attractors;
    thrd_info[tx].convergences = convergences;
    thrd_info[tx].ib = tx;
    thrd_info[tx].d = &d;
    thrd_info[tx].conv_cap = &conv_cap;
    thrd_info[tx].h = &delta_h;
    thrd_info[tx].istep = &threads;
    thrd_info[tx].tx = tx;
    thrd_info[tx].lines = &lines;
    thrd_info[tx].roots = roots;
    thrd_info[tx].status = status;
    thrd_info[tx].cnd = &cnd;
    thrd_info[tx].mtx = &mtx;
    status[tx].val = 0;
    int r = thrd_create(thrds + tx, main_thrd, (void*) (thrd_info + tx));
    if( r != thrd_success ){
      fprintf(stderr,"Failed to create thread\n");
      exit(1);
    }
    thrd_detach(thrds[tx]);
  }
  
  /* Start writing thread */
  {
    write_info.attractors = attractors;
    write_info.convergences = convergences;
    write_info.grayscale = (const char**)grayscale;
    write_info.d = &d;
    write_info.conv_cap = &conv_cap;
    write_info.lines = &lines;
    write_info.nthrds = &threads;
    write_info.attr_file = attr_file;
    write_info.conv_file = conv_file;
    write_info.status = status;
    write_info.cnd = &cnd;
    write_info.mtx = &mtx;
    int r = thrd_create(&thrd_write, main_thrd_write, (void*) &write_info);
    if( r != thrd_success ){
      fprintf(stderr,"Failed to create write thread\n");
      exit(1);
    }
  }
  {
    int r;
    thrd_join(thrd_write, &r);
  }

  /* Cleanup */
  mtx_destroy(&mtx); cnd_destroy(&cnd);
  free(attractors); free(convergences);
  free(grayscale); free(grayscale_ent);
  return 0;
}
