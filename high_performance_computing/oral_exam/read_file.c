#include<stdio.h>
#include<stdlib.h>
#include<string.h>

void parse_line(char* line, float** h_old){
  char delim[2] = " ";
  char* index1 = strtok(line, delim);
  char* index2 = strtok(NULL, delim);
  char* value = strtok(NULL, delim);

  int x = atoi(index1)+1; // Shift with one to compensate for border
  int y = atoi(index2)+1;
  float h = atof(value);

  h_old[x][y] = h;
}

/* 
 * Reads the header information for the initialization
 * of the grid using strtok. 
 */
void read_header(FILE* file, int *nx, int *ny){
  fseek(file, 0, SEEK_END);
  int size = ftell(file);
  fseek(file,0,SEEK_SET);
  
  char* data =(char*) malloc(sizeof(char)*size);
  fread(data, sizeof(char), size, file);

  char* line = strtok(data, "\n"); // Split data on newline
  
  char* cnx;
  char* cny;
  const char delim[2]=" ";

  cnx = strtok(line, delim); // read first two lines
  cny = strtok(NULL, delim);


  *nx = atoi(cnx);
  *ny = atoi(cny);

  free(data);
}

void read_data(FILE* file, float** h_old)
{
  /* Read data again - a bit unnecessary but oh well */
  fseek(file, 0, SEEK_END);
  int size = ftell(file);
  fseek(file,0,SEEK_SET);
  
  char* data =(char*) malloc(sizeof(char)*size);
  fread(data, sizeof(char), size, file);
  
  /* Skip the first two lines */
  char* line = strtok(data, "\n"); // first line
  line = strtok(NULL, "\n");

  /* Read the rest of the file into h_old*/
  do{
    parse_line(line, h_old);
    line = strtok(NULL, "\n");
  }while(line != NULL);
  

  free(data);
}
