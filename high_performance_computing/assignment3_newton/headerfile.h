#include<complex.h>
struct Args parsing_cmdl_arg(int argc, char* argv);
void newton_iteration(char* atractor,  char* convergence, double complex x0, int d, int cap, const double complex* roots);
FILE* initialize_attractor_file(int lines, int d);
FILE* initialize_convergence_file(int lines, int d, int convergence_cap);
void write_row_to_file(char* attractor, char* convergences, FILE* attr_file, FILE* conv_file, int lines, int d, int convergence_cap, const char** grayscale);
void close_files(FILE* attr_file, FILE* conv_file);

