// this include is needed, remember to add the gsl library when you compile as in E1.
#include <gsl/gsl_rng.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main() {
	// declarations and initialization
	double u;
	const gsl_rng_type *T;
	gsl_rng *q;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	q = gsl_rng_alloc(T);
	gsl_rng_set(q,time(NULL));

	// generate uniform random number
	u = gsl_rng_uniform(q);

	printf("%e\n",u);

	// free memory
	gsl_rng_free(q);
}
