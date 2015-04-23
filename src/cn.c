#include "cn.h"

#include "common.h"
#include <cblas.h>
#include <lapacke.h>

/**
 * random_pts_in_box - samples random points in a box
 * @m:      Dimension of the space.
 * @l:      Number of points to sample.
 * @xh:     An m-dimensional vector.
 * @v:      An m-dimensional vector.
 * @X:      The (m + l)-by-l matrix where the result is written.
 *
 * Samples l points uniformly at random in the box
 *   { x : abs((x(i) - xh(i)) / (xh(i) * v(i))) < 1 }.
 * x has dimension m.
 *
 * An l-by-l identity matrix is used to pad X. See cluster_newton() for
 * the rationale behind this decision.
 */
void random_pts_in_box(uint m, uint l, float *xh, float *v, float *X)
{
	assert (m != 0);
	assert (l != 0);

	for (uint i = 1; i <= l; i++) {
		for (uint k = 1; k <= m; k++) {
			float r = 2. * (float)rand() / (float)RAND_MAX - 1.;
			M_IDX(X, m + l, k, i) = V_IDX(xh, k) *
			                        (1. + V_IDX(v, k) * r);
		}

		for (uint k = 1; k <= l; k++) {
			M_IDX(X, m + l, m + k, i) = (k == i ? 1. : 0.);
		}
	}
}

/**
 * perturbate() - creates a perturbated target vector
 * @l:             Number of perturbed target vectors.
 * @n:             Dimension of the target vector.
 * @ys:            Target vector.
 * @eta:           Relative magnitude of the random perturbation.
 * @Ys:            The n-by-l matrix where the result is written.
 *
 * The resulting matrix ys satisfies:
 *   max (ys(i,j) - ys(i) / ys(i)) <= eta,
 * where the max is taken over i = 1, 2, ..., n.
 */
void perturbate(uint l, uint n, float *ys, float eta, float *Ys)
{
	for (uint j = 1; j <= l; j++) {
		for (uint i = 1; i <= n; i++) {
			float r = eta * (2. * (float)rand()
					 / (float)RAND_MAX - 1.);
			M_IDX(Ys, n, i, j) = V_IDX(ys, i) * (1. + r);
		}
	}
}

/**
 * multi_eval() - evaluates a function at multiple points
 * @m:              Number of parameters of the function.
 * @f:              The function to evaluate.
 * @l:              Number of points.
 * @X:              Coordinates of the points, one point per column.
 * @y:              Vector in which to store the result.
 *
 * This function assumes COLUMN-MAJOR ORDER.
 */
void multi_eval(uint m, uint n, void (*f)(float *, float *),
		uint l, float *X, float *Y)
{
	for (uint j = 1; j <= l; j++) {
		f(M_COL(X, m + l, j), M_COL(Y, n, j));
	}
}

/**
 * least_squares() - solve an overdetermined linear system
 * @m:                 Row dimension of A.
 * @n:                 Column dimension of A.
 * @A:                 An m-by-n matrix.
 * @l:                 Column dimension of B.
 * @B:                 An n-by-l matrix.
 * @X:                 An n-by-m matrix in which to store the result.
 *
 * Solves XA = B in the sense of least squares.
 */
void least_squares(uint m, uint n, float *A, uint l, float *B, float *X)
{
	float *C = create_matrix(m, m);
	float *D = X;

	//print_matrix(m, n, A);
	//print_matrix(l, n, B);

	//printf("m=%u n=%u\n", m, n);
	/* C = AA' */
	cblas_ssyrk(CblasColMajor, CblasLower, CblasNoTrans,
	            m, n, 1.0f, A, m, 0.0f, C, m);
	//print_matrix(m, m, C);
	/* D = BA' */
	cblas_sgemm(CblasColMajor, CblasNoTrans, CblasTrans,
	            l, m, n, 1.0f, B, l, A, m, 0.0f, D, l);
	//print_matrix(l, m, D);
	/* Solve XC=D for X. This overrides D, hence X. */
	int *ipiv = (int*)malloc(sizeof(int) * m);
	assert(ipiv);
	LAPACKE_ssysv(LAPACK_ROW_MAJOR, 'u', m, l, C, m, ipiv, D, l);

	free(ipiv);
	free(C);
}

/**
 * minimum_norm() - solve an underdetermined linear system
 * @m:                Number of equations.
 * @n:                Number of unknowns.
 * @A:                LHS, an m-by-n matrix.
 * @l:                Column dimension of the RHS.
 * @B:                RHS, an m-by-l matrix.
 * @X:                An n-by-l matrix in which the result is stored.
 *
 * Finds the solution of AX = B of minimum Frobenius norm. B is modified.
 */
void minimum_norm(uint m, uint n, float *A, uint l, float *B, float *X)
{
	float *C = create_matrix(m, m);

	/* C = AA' */
	cblas_ssyrk(CblasColMajor, CblasUpper, CblasNoTrans,
	            m, n, 1.0f, A, m, 0.0f, C, m);
	//print_matrix(m, m, C);
	/* Solve CX = B */
	int *ipiv = (int*)malloc(sizeof(int) * m);
	assert(ipiv);
	LAPACKE_ssysv(LAPACK_COL_MAJOR, 'u', m, l, C, m, ipiv, B, m);
	//print_matrix(m, l, B);

	/* Multiply the result by A' on the left */
	cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, n, l, m, 1.0f, A,
	            m, B, m, 0.0f, X, n);

	// clean up
	free(ipiv);
	free(C);
}

void cluster_newton(uint m, uint n, void (*f)(float *, float *), float *ys,
                    float *xh, float *v,
                    uint l, float eta, uint K)
{
	// safety checks
	assert(m > n);
	assert(n > 0);
	assert(l > 0);

	/* 1.1 */ float *X = create_matrix(m + l, l);
	random_pts_in_box(m, l, xh, v, X);

	/* 1.2 */ float *Ys = create_matrix(n, l);
	perturbate(l, n, ys, eta, Ys);

	float *Y = create_matrix(n, l);

	/* A and Y0 are stored in the same matrix. This comes handy when one
	 * has to solve the overdetermined linear system in 2.2. */
	float *A_Y0 = create_matrix(n, m + l);
	float *A = A_Y0;
	float *Y0 = M_COL(A_Y0, n, m + 1);
	float *S = create_matrix(m, l);

	for (uint k = 0; k <= K; k++) {
		/* 2.1 */ multi_eval(m, n, f, l, X, Y);

		/* 2.2 */ least_squares(m + l, l, X, n, Y, A_Y0);

		/* 2.3 */
		// Y0 <-- Ys - AX - Y0
		cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
			    n, l, m, -1.0f, A, n, X, m + l, -1.0f, Y0, n);
		m_add(n, l, Y0, Ys);

		m_scale_cols(n, m, A, xh);
		minimum_norm(n, m, A, l, Y0, S);
		m_scale_rows_inv(n, m, S, xh);

		/* 2.4 */
		for (uint j = 1; j <= l; j++) {
			while (0) { // FIXME
				for (uint i = 1; i <= m; i++) {
					M_IDX(S, m, i, j) /= 2.0f;
				}
			}
		}
		m_add(m, l, X, S);
	}

	// cleaning up
	free(S);
	free(A_Y0);
	free(Y);
	free(Ys);
	free(X);
}
