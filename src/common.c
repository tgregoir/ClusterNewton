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
#include "common.h"

#include <lapacke.h>

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
		printf("%e; ", V_IDX(v, i));
	}
	printf("];\n");
}

void print_matrix_(uint m, uint n, float *A, const char *str)
{
	printf("%s = ...\n [ ", str);
	for (uint i = 1; i <= m; i++) {
		for (uint j = 1; j <= n; j++) {
			printf("%e ", M_IDX(A, m, i, j));
		}
		if (i < m) {
			printf(";\n  ");
		}
	}
	printf("];\n");
}

/**
 * m_copy() - replaces a matrix with a copy of another one
 * @m:         Row dimension.
 * @l:         Column dimension.
 * @ldA:       Leading dimension of matrix A.
 * @A:         Target.
 * @ldB:       Leading dimension of matrix B.
 * @B:         Matrix to be added to A.
 *
 * Performs A <- B.
 */
void m_copy(uint m, uint n, uint ldA, float *A, uint ldB, float *B)
{
	for (uint j = 1; j <= n; j++) {
		for (uint i = 1; i <= m; i++) {
			M_IDX(A, ldA, i, j) = M_IDX(B, ldB, i, j);
		}
	}
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
 * m_add() - subtracts a matrix from another one
 * @m:         Row dimension.
 * @l:         Column dimension.
 * @ldA:       Leading dimension of matrix A.
 * @A:         Target.
 * @ldB:       Leading dimension of matrix B.
 * @B:         Matrix to be subtracted to A.
 *
 * Performs A <- A - B.
 */
void m_sub(uint m, uint l, uint ldA, float *A, uint ldB, float *B)
{
	for (uint j = 1; j <= l; j++) {
		for (uint i = 1; i <= m; i++) {
			M_IDX(A, ldA, i, j) -= M_IDX(B, ldB, i, j);
		}
	}
}

/**
 * m_scale() - mutliplies a matrix by a constant
 * @m:               Row dimension of A.
 * @n:               Column dimension of A.
 * @ldA:             Leading dimension of A.
 * @k:               Constant.
 *
 * Sets A <- kA.
 */
void m_scale(uint m, uint n, uint ldA, float *A, float k)
{
	for (uint j = 1; j <= n; j++) {
		for (uint i = 1; i <= m; i++) {
			M_IDX(A, ldA, i, j) *= k;
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

/**
 * m_replicate() - puts several copies of a column vector into a matrix
 * @n:               Dimension of the vector/row dimension of the matrix.
 * @v:               The vector.
 * @m:               Column dimension of the matrix.
 * @A:               Where to put the result.
 *
 * Puts m copies of v into A.
 */
void m_replicate(uint n, float *v, uint m, float *A)
{
	for (uint j = 1; j <= m; j++) {
		for (uint i = 1; i <= n; i++) {
			M_IDX(A, n, i, j) = V_IDX(v, i);
		}
	}
}

/**
 * m_transpose() - matrix transposition
 * @n:               Row dimension of B.
 * @m:               Column dimension of B.
 * @A:               Where to store the result (m-by-n).
 * @B:               Input matrix (n-by-m).
 *
 * Sets A <- B'. This is a naive algorithm with poor cache performance.
 */
void m_transpose(uint n, uint m, float *A, float *B)
{
	for (uint j = 1; j <= m; j++) {
		for (uint i = 1; i <= n; i++) {
			M_IDX(A, m, j, i) = M_IDX(B, n, i, j);
		}
	}
}

/**
 * v_norm() - 2-norm of a vector
 * @n:          Row dimension of v.
 * @v:          Input vector.
 *
 * Returns ||v||_2.
 */
float v_norm(uint n, float *v)
{
	return LAPACKE_slange(LAPACK_COL_MAJOR, 'F', n, 1, v, n);
}
