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

void print_vector_(uint n, float *v, const char *str)
{
	printf("%s = [ ", str);
	for (uint i = 1; i <= n; i++) {
		printf("%f; ", V_IDX(v, i));
	}
	printf("]\n");
}

void print_matrix_(uint m, uint n, float *A, const char *str)
{
	printf("%s = ...\n [ ", str);
	for (uint i = 1; i <= m; i++) {
		for (uint j = 1; j <= n; j++) {
			printf("%f ", M_IDX(A, m, i, j));
		}
		if (i < m) {
			printf(";\n  ");
		}
	}
	printf("]\n");
}

/**
 * m_add() - add a matrix to another one
 * @m:         Row dimension.
 * @l:         Column dimension.
 * @ldA:       Leading dimension of matrix A.
 * @A:         Target.
 * @ldB:       Leading dimension of matrix B.
 * @B:         Matrix to be added to A.
 *
 * Performs A <- A + B.
 */
void m_add(uint m, uint l, uint ldA, float *A, uint ldB, float *B)
{
	for (uint j = 1; j <= l; j++) {
		for (uint i = 1; i <= m; i++) {
			M_IDX(A, ldA, i, j) += M_IDX(B, ldB, i, j);
		}
	}
}

/**
 * m_scale_cols() - right multiplication with a diagonal matrix
 * @m:                Row dimension of the matrix A.
 * @n:                Column dimension of the matrix A.
 * @A:                A matrix
 * @sv:               A length-n vector.
 *
 * Performs A <- A diag(sv).
 */
void m_scale_cols(uint m, uint n, float *A, float *sv)
{
	for (uint j = 1; j <= n; j++) {
		float s = V_IDX(sv, j);
		for (uint i = 1; i <= m; i++) {
			M_IDX(A, m, i, j) *= s;
		}
	}
}

/**
 * m_scale_rows_inv() - left product with the inverse of a diagonal matrix
 * @m:                Row dimension of the matrix A.
 * @n:                Column dimension of the matrix A.
 * @A:                A matrix
 * @sv:               A length-m vector.
 *
 * Performs A <- diag(sv)^(-1) A.
 */
void m_scale_rows_inv(uint m, uint n, float *A, float *sv)
{
	for (uint j = 1; j <= n; j++) {
		for (uint i = 1; i <= m; i++) {
			M_IDX(A, m, i, j) /= V_IDX(sv, i);
		}
	}
}
