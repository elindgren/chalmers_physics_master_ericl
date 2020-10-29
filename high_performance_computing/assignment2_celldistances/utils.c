#include<stdio.h>
#include<stdlib.h>

int calc_dist_sq(int* a, int* b){
  // Computes the distance squared between 3D points a and b.
  // Note that it is assumed that a and b are vectors of length 3.
  int d0 = a[0] - b[0]; 
  int d1 = a[1] - b[1]; 
  int d2 = a[2] - b[2];
  return d0*d0 + d1*d1 + d2*d2;
}

void parse_string_point(int* coordinates, char* point){
/*
  input 
  - coordinates is the place where to store the cordinates of the point
  - point is the string with coordinates read from the data file

*/
  int jx;
  for(int ix = 0, jx = 0; jx <3; ix += 8, jx++){
    int sign  = -1*(((int) point[0 + ix]) - 44); // + 1 for + and -1 for -
    int d1000 = ((int) point[2 + ix] - 48) * 1000;
    int d100  = ((int) point[4 + ix] - 48) * 100;
    int d10   = ((int) point[5 + ix] - 48) * 10;
    int d1    = ((int) point[6 + ix] - 48) * 1;

    coordinates[jx] = sign*( d1000 + d100 + d10 + d1);
  }
}

int read_data(int** matrix, int start_line, int nbr_lines, FILE* file){
/*
  indata
  - file is the file where the data should be read from
  - matrix where the read points are stored
  - start_line is where the program starts reading the file
  - nbr of lines we want the function to read

  Returns the number of lines that were read.
*/
  int char_per_line = 24;
  
  /* Temporary data structures */
  char *data_entries = (char*) malloc(sizeof(char)*char_per_line*nbr_lines);
  char **data = (char**) malloc(sizeof(char*)*nbr_lines);

  fseek(file, start_line*char_per_line, SEEK_SET); // Move file cursor

  for(int ix = 0; ix < nbr_lines; ix++){ // initialize contigously
    data[ix] = data_entries + ix*char_per_line;
  }
  
  fread(data_entries, sizeof(char), char_per_line*nbr_lines, file);

  for(size_t row = 0; row<nbr_lines; row++){
    parse_string_point(matrix[row], data[row]);
  }

  free(data_entries);
  free(data);
  int last_line = ftell(file)/char_per_line;
  int nbr_of_lines_read = last_line - start_line;
  return nbr_of_lines_read;
}

