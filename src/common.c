#include "linalg.h"

#include "common.h"

/**
 * create_vector() - memory allocation for vectors
 * @n:        Dimension.
 *
 * Return: a pointer to a non-initialized vector.
 */
float *create_vector(uint n)
{
	assert(n);
	float *M = (float *)malloc(sizeof(float) * n);
	assert(M);
	return M;
}

/**
 * create_matrix() - memory allocation for matrices
 * @n:        Row dimension.
 * @m:        Column dimension.
 *
 * Return: a pointer to a non-initialized matrix.
 */
float *create_matrix(uint n, uint m)
{
	assert(n);
	assert(m);
	float *M = (float *)malloc(sizeof(float) * n * m);
	assert(M);
	return M;
}
