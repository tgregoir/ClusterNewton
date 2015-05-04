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
#define V_IDX(v, i) ((v)[(i) - 1])

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
#define M_IDX(A, n, i, j) ((A)[(j - 1) * (n) + (i) - 1])

/**
 * M_COL() - get a pointer to a column in a matrix
 * @A:         Matrix, a float *.
 * @n:         Row dimension of A.
 * @j:         Column index.
 *
 * M is assumed to be written in column-major order. Bounds are NOT checked.
 *
 * Returns: a pointer to A(., j).
 */
#define M_COL(A, n, j) (&(M_IDX(A, n, 1, j)))

float *create_vector(uint);
float *create_matrix(uint, uint);

void print_vector_(uint, float *, const char *);
void print_matrix_(uint, uint, float *, const char *);

/**
 * print_vector() - prints a vector
 * @n:                Its size.
 * @v:                The vector.
 */
#define print_vector(n, v) do { print_vector_(n, v, #v); } while (0)

/**
 * print_matrix() - prints a matrix
 * @m:                Its row dimension.
 * @n:                Its column dimension.
 * @A:                The matrix.
 */
#define print_matrix(m, n, A) do { print_matrix_(m, n, A, #A); } while (0)

void m_add(uint, uint, float *, float *);
void m_scale_cols(uint, uint, float *, float *);
void m_scale_rows_inv(uint, uint, float *, float *);

#endif /* COMMON_H */
