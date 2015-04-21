#ifndef COMMON_H
#define COMMON_H

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

typedef unsigned int uint;

/**
 * V_IDX() - indexing function for vectors
 * @v:         Vector, represented by a float *.
 * @i:         Index, with 1 <= i <= d, where d is the dimension of v.
 *
 * Bounds are NOT checked.
 *
 * Return: v(i).
 */
#define V_IDX(v, i) (v[i - 1])

/**
 * M_IDX() - indexing function for matrices
 * @A:        Matrix, represented by a float *.
 * @n:        Row dimension of A.
 * @i:        Index, with 1 <= i <= n.
 * @j:        Index, with 1 <= j <= m, where m is the column dimension of A.
 *
 * A is assumed to be written in column-major order. Bounds are NOT checked.
 *
 * Return: A(i,j).
 */
#define M_IDX(A, n, i, j) (A[(j - 1) * (n) + i - 1])

/**
 * M_COL() - get a pointer to a column in a matrix
 * @M:         Matrix, a float *.
 * @n:         Row dimension of M.
 * @j:         Column index.
 *
 * M is assumed to be written in column-major order. Bounds are NOT checked.
 *
 * Returns: a pointer to A(., j).
 */
#define M_COL(M, n, j) (&M_IDX(M, n, 1, j))

float *create_vector(uint n);
float *create_matrix(uint n, uint m);

#endif /* COMMON_H */
