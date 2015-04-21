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
 * @X:      The m-by-l matrix where the result is written.
 *
 * Samples l points uniformly at random in the box
 *   { x : abs((x(i) - xh(i)) / (xh(i) * v(i))) < 1 }.
 * x has dimension m.
 */
void random_pts_in_box(uint m, uint l, float *xh, float *v, float *X)
{
	assert (m != 0);
	assert (l != 0);

	for (uint i = 1; i <= l; i++) {
		for (uint k = 1; k <= m; k++) {
			float r = 2. * (float)rand() / (float)RAND_MAX - 1.;
			M_IDX(X, m, k, i) = V_IDX(xh, k) *
			                    (1. + V_IDX(v, k) * r);
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
		f(M_COL(X, m, j), M_COL(Y, n, j));
	}
}

/**
 * least_squares() - solve an overdetermined linear system
 * @m:                 Row dimension of the matrix.
 * @n:                 Column dimension of the matrix.
 * @A:                 LHS.
 * @p:                 Column dimension of the RHS.
 * @X:                 Matrix in which to store the result.
 *
 *
 */
void least_squares(uint m, uint n, float *A, uint p, float *B, float *X)
{
	float *C = create_matrix(n, n);
	float *D = X;

	/* C = A'A */
	cblas_ssyrk(CblasColMajor, CblasUpper, CblasTrans,
	            n, m, 1.0f, A, m, 0.0f, C, n);
	/* D = A'B */
	cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans,
	            n, p, m, 1.0f, A, m, B, m, 0.0f, D, n);
	/* Solve CX = D for X. This overrides D, hence X. */
	int *ipiv = (int*)malloc(sizeof(int) * n);
	assert(ipiv);
	LAPACKE_ssysv(LAPACK_COL_MAJOR, 'u', n, p, C, n, ipiv, D, n);

	free(ipiv);
	free(C);
}

void fit_hyperplane(uint m, uint n, uint l, float *X, float *Y,
		    float *A, float *Y0)
{
}

void cluster_newton(uint m, uint n, void (*f)(float *, float *), float *ys,
                    float *xh, float *v,
                    uint l, float eta, uint K)
{
	// safety checks
	assert(m > n);
	assert(n > 0);
	assert(l > 0);

	// 1.1.
	float *X = create_matrix(m, l);
	random_pts_in_box(m, l, xh, v, X);

	// 1.2.
	float *Ys = create_matrix(n, l);
	perturbate(l, n, ys, eta, Ys);

	float *Y = create_matrix(n, l);
	for (uint k = 0; k <= K; k++) {
		// 2.1.
		multi_eval(m, n, f, l, X, Y);

		// 2.2.
		//fit_hyperplane(m, n, l, X, Y, A, Y0);
	}

	// cleaning up
	free(Y);
	free(Ys);
	free(X);
}
