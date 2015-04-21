#ifndef TSTTOOLS_H
#define TSTTOOLS_H

#include "common.h"
#include <time.h>

/**
 * init_prg() - initialize the pseudorandom number generator
 *
 * Takes the current time as the seed and prints it to standard output, to
 * make the experiment reproducible.
 */
static inline void init_prg(void)
{
	time_t seed = time(NULL);
	printf("Seed: %d\n", (int) seed);
	srand(seed);
}

#define MAX_RAND_DIM 1024

static inline uint random_dim(void)
{
	uint d = (uint)(rand() % (MAX_RAND_DIM + 1));
	return (d != 0 ? d : 1);
}

static inline float random_float(void)
{
	return (float)rand();
}

static inline float *random_vector(uint n)
{
	float *v = create_vector(n);
	for (uint k = 1; k <= n; k++) {
		V_IDX(v, k) = random_float();
	}
	return v;
}

static inline float *random_matrix(uint n, uint m)
{
	float *M = create_matrix(n, m);
	for (uint i = 1; i <= n; i++) {
		for (uint j = 1; j <= m; j++) {
			M_IDX(M, n, i, j) = random_float();
		}
	}
	return M;
}

#endif /* TSTTOOLS_H */
