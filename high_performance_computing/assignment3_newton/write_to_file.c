#include <stdio.h>
#include <complex.h>
#include <stdlib.h>
#include <math.h>
const char colors[][6] = {
		"1 0 0 ",
		"0 1 0 ",
		"0 0 1 ",
		"1 1 0 ",
		"1 0 1 ",
		"0 1 1 ",
		"2 1 0 ",
		"2 0 1 ",
		"0 2 1 ",
		"1 2 0 ",
		"1 0 2 ",
		"0 1 2 "
};

int nbrDigits(int n){
  return (int)(log10(n)+1);
}

FILE* initialize_attractor_file(int lines, int d){
  FILE* attr_file;
  char attr_file_name[40];

  sprintf(attr_file_name, "newton_attractors_x%d.ppm", d);

  // Open files 
  attr_file = fopen(attr_file_name, "w");
  
  /* Write header */
  // For the attractors, we need d different colors since that is the number of roots from 
  // the fundamental theorem of algebra.
  int colordepth = 2;
  int digitsInLines = nbrDigits(lines);
  int digitsInColordepth = nbrDigits(colordepth);
  char header[2*digitsInLines+6+digitsInColordepth];
  sprintf(header, "P3\n%d %d\n%d\n", lines, lines, colordepth);
  fwrite(header, sizeof(char), sizeof(header), attr_file);
  

  return attr_file;
}       

FILE* initialize_convergence_file(int lines, int d, int convergence_cap){
  FILE* conv_file; 
  char conv_file_name[40];

  sprintf(conv_file_name, "newton_convergence_x%d.ppm", d);

  // Open files 
  conv_file = fopen(conv_file_name, "w");

  // Write headers
  int digitsInLines = nbrDigits(lines);
  int digitsInGreyDepth = nbrDigits(convergence_cap);
  char header[2*digitsInLines+6+digitsInGreyDepth];
  sprintf(header, "P3\n%d %d\n%d\n", lines, lines, convergence_cap);
  fwrite(header, sizeof(char), sizeof(header), conv_file);
  

  return conv_file;
}       

void write_row_to_file(char* attractor, char* convergence, FILE* attr_file, FILE* conv_file, int lines, int d, int convergence_cap, const char** grayscale){
  // We may assume both the number of colors, and the number of convergence iterations
  // as small, i.e. in the range < 100. Hence they will fit in a standard 1 byte char. 
  // When writing to file, we need twice the size of a line to make size for whitespace.
 
  int a_size = lines*6+1;
  int c_size = lines*12+1;
  char *attr_row = (char*)malloc(sizeof(char)*a_size); 
  char *conv_row = (char*)malloc(sizeof(char)*c_size); 
  for(int i=0; i<lines; i++){ 
    char c = convergence[i];
    for(int ix=0; ix<6; ix++){ 
      attr_row[i*6 + ix] = colors[(int)attractor[i]][ix];  
      conv_row[i*12 + 2*ix] = grayscale[c][2*ix];
      conv_row[i*12 + 2*ix+1] = grayscale[c][2*ix+1];
    }
  }
  
  attr_row[a_size-1] = 10; //\n
  conv_row[c_size-1] = 10;
  fwrite(attr_row, sizeof(char), a_size, attr_file); 
  fwrite(conv_row, sizeof(char), c_size, conv_file); 
  
  free(attr_row); free(conv_row); 
}

void close_files(FILE* attr_file, FILE* conv_file){
  // Cleanup
  fclose(attr_file); 
  fclose(conv_file);
}
