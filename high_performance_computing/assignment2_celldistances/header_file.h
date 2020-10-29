#include <x86intrin.h>

int read_data(int** matrix, int start_line, int nbr_lines, FILE* file);
int calc_dist_sq(int* a, int* b);
int parsing_cmdl_arg(int argc, char** argv);









/* Legacy code */

/*
 * Slow implementation of vectorized version of calculating distances.
static inline
int calc_dist_sq1(int a[4], int b[4]){
  // Vectorized version of calc_dist_sq below
  int d[4] = {0, 0, 0, 0};
  __m128i ar = _mm_loadu_si128((__m128i*)(a));
  __m128i br = _mm_loadu_si128((__m128i*)(b));
  // perform subtraction
  __m128i cr = _mm_sub_epi32(ar, br);
  // Multiply, i.e. square
  __m128i dr = _mm_mul_epi32(cr, cr);

  _mm_store_si128((__m128i*)(d), dr); 
  
  return d[0] + d[1] + d[2];
}
*/
