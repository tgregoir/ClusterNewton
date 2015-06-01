/*
 *    This file is part of CNewt.
 *
 *    CNewt is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    CNewt is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with CNewt.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef TSTTOOLS_H
#define TSTTOOLS_H

#ifdef __cplusplus
extern "C" {
#endif

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

#ifdef __cplusplus
}
#endif

#endif /* TSTTOOLS_H */
