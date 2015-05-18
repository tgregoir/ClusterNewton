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
#include "cn.h"

#include "common.h"
#include <cblas.h>
#include <lapacke.h>
#include <math.h>

/**
 * random_pts_in_box - samples random points in a box
 * @m:      Dimension of the space.
 * @l:      Number of points to sample.
 * @xh:     An m-dimensional vector.
 * @v:      An m-dimensional vector.
 * @X:      The (m + 1)-by-l matrix where the result is written.
 *
 * Samples l points uniformly at random in the box
 *   { x : abs((x(i) - xh(i)) / (xh(i) * v(i))) < 1 }.
 * x has dimension m.
 *
 * A 1-by-l block of ones is used to pad X. See cluster_newton() for
 * the rationale behind this decision.
 */
void random_pts_in_box(uint m, uint l, float *xh, float *v, float *X)
{
	for (uint j = 1; j <= l; j++) {
		for (uint i = 1; i <= m; i++) {
			float r = 2. * (float)rand() / (float)RAND_MAX - 1.;
			M_IDX(X, m + 1, i, j) = V_IDX(xh, i) *
			                        (1. + V_IDX(v, i) * r);
		}
		M_IDX(X, m + 1, m + 1, j) = 1.;
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
		f(M_COL(X, m + 1, j), M_COL(Y, n, j));
	}
}

/**
 * pinv_ls() - solve an overdetermined linear system
 * @m:                 Row dimension of A.
 * @n:                 Column dimension of A, n <= m.
 * @A:                 An m-by-n matrix.
 * @l:                 Column dimension of B.
 * @B:                 An l-by-n matrix.
 * @X:                 An l-by-m matrix in which to store the result.
 *
 * Computes the least-squares solution of an overdetermined linear system,
 * by computing a Moore-Penrose pseudoinverse. Solves XA = B for X.
 */
void pinv_ls(uint m, uint n, float *A, uint l, float *B, float *X)
{
	assert (n <= m);

	/* create a copy of A */
	float *cA = create_matrix(m, n);
	for (uint i = 1; i <= m; i++) {
		for (uint j = 1; j <= n; j++) {
			M_IDX(cA, m, i, j) = M_IDX(A, m, i, j);
		}
	}

	float *invA = create_matrix(m, m);

	/* set up an identity function */
	for (uint i = 1; i <= m; i++) {
		for (uint j = 1; j <= m; j++) {
			M_IDX(invA, m, i, j) = (i == j ? 1.0f : 0.0f);
		}
	}

	/* compute the pseudoinverse */
	float *S = create_vector(m);
	int rank;
	LAPACKE_sgelss(LAPACK_COL_MAJOR, m, n, m, cA, m, invA, m, S,
	               -1.0f, &rank);

	/* multiply the RHS by the pseudoinverse */
	cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
	            l, m, n, 1.0f, B, l, invA, m, 0.0f, X, l);

	free(S);
	free(invA);
	free(cA);
}

/**
 * normal_ls() - solve an overdetermined linear system
 * @m:                 Row dimension of A.
 * @n:                 Column dimension of A.
 * @A:                 An m-by-n matrix.
 * @l:                 Column dimension of B.
 * @B:                 An l-by-n matrix.
 * @X:                 An l-by-m matrix in which to store the result.
 *
 * Computes the least-squares solution of an overdetermined linear system,
 * by solving the normal equations. Solves XA = B for X.
 */
void normal_ls(uint m, uint n, float *A, uint l, float *B, float *X)
{
	float *C = create_matrix(m, m);
	float *D = create_matrix(m, l);

	/*
	 * Dimensions:
	 *   - A is m by n.
	 *   - B is l by n.
	 *   - C is m by m.
	 *   - D is m by l.
	 */

	/* C = AA' */
	cblas_ssyrk(CblasColMajor, CblasLower, CblasNoTrans,
	            m, n, 1.0f, A, m, 0.0f, C, m);
	/* D = AB' */
	cblas_sgemm(CblasColMajor, CblasNoTrans, CblasTrans,
	            m, l, n, 1.0f, A, m, B, l, 0.0f, D, m);
	/* solve CX=D for X
	 * this overrides D */
	int *ipiv = (int*)malloc(sizeof(int) * m);
	assert(ipiv);
	LAPACKE_ssysv(LAPACK_COL_MAJOR, 'l', m, l, C, m, ipiv, D, m);

	/* transpose the result */
	m_transpose(m, l, X, D);

	free(ipiv);
	free(D);
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

	/*
	 * Dimensions
	 *   - A is m by n.
	 *   - B is m by l.
	 *   - C is m by m.
	 */

	/* C = AA' */
	cblas_ssyrk(CblasColMajor, CblasUpper, CblasNoTrans,
	            m, n, 1.0f, A, m, 0.0f, C, m);
	/* solve CX = B */
	int *ipiv = (int*)malloc(sizeof(int) * m);
	assert(ipiv);
	LAPACKE_ssysv(LAPACK_COL_MAJOR, 'u', m, l, C, m, ipiv, B, m);

	/* multiply the result by A' on the left */
	cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans,
	            n, l, m, 1.0f, A, m, B, m, 0.0f, X, n);

	/* clean up */
	free(ipiv);
	free(C);
}

/**
 * cluster_newton() - the cluster Newton method to solve inverse problems
 * @m:      Dimension of the parameter space.
 * @n:      Dimension of the result space.
 * @f:      A function that maps vectors of size m to vectors of size n.
 * @ys:     Target vector, of dimension n.
 * @xh:     Center of the initial box. Vector of size m.
 * @v:      Relative size of the initial box. Vector of size m.
 * @l:      Number of candidates to generate.
 * @eta:    Target accuracy.
 * @K:      Number of iterations.
 * @Xf:     Where to store the result. Matrix of size l by m.
 * @r:      Where to store the residuals. Vector of size l. Can also be
 *          set to NULL, if the user does not need to compute them.
 */
void cluster_newton(uint m, uint n, void (*f)(float *, float *), float *ys,
                    float *xh, float *v,
                    uint l, float eta, uint K, float *Xf, float * r)
{
	/* safety checks */
	assert(m > n);
	assert(n > 0);
	assert(l > 0);

	/* 1.1 */ float *X = create_matrix(m + 1, l);
	random_pts_in_box(m, l, xh, v, X);

	/* 1.2 */ float *Ys = create_matrix(n, l);
	perturbate(l, n, ys, eta, Ys);

	float *Y = create_matrix(n, l);

	/* A and y0 are stored in the same matrix
	 * handy when solving the overdetermined linear system in 2.2 */
	float *A_y0 = create_matrix(n, m + 1);
	float *A = A_y0;
	float *y0 = M_COL(A_y0, n, m + 1);
	float *Y0 = create_matrix(n, l);
	float *S = create_matrix(m, l);

	for (uint k = 0; k <= K; k++) {
		/* 2.1 */ multi_eval(m, n, f, l, X, Y);

		/* 2.2 */ normal_ls(m + 1, l, X, n, Y, A_y0);
		/* 2.2 */ //pinv_ls(l, m + 1, X, n, Y, A_y0);
		m_replicate(n, y0, l, Y0);

		/* 2.3 */
		/* Y0 <-- Ys - AX - Y0 */
		cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
			    n, l, m, -1.0f, A, n, X, m + 1, -1.0f, Y0, n);
		m_add(n, l, n, Y0, n, Ys);

		m_scale_cols(n, m, A, xh);
		minimum_norm(n, m, A, l, Y0, S);
		m_scale_rows_inv(m, l, S, xh);

		/* 2.4 */
		for (uint j = 1; j <= l; j++) {
			/* FIXME */
			while (0) {
				for (uint i = 1; i <= m; i++) {
					M_IDX(S, m, i, j) /= 2.0f;
				}
			}
		}
		m_add(m, l, m + 1, X, m, S);
	}

	/* copy the result */
	m_copy(m, l, m, Xf, m + 1, X);

	/* compute the residuals */
	if (r) {
		for (uint j = 1; j <= l; j++) {
			V_IDX(r, j) = 0.0f;
			for (uint i = 1; i <= n; i++) {
				float rel = (M_IDX(Y, n, i, j) - V_IDX(ys, i))
				            / V_IDX(ys, i);
				V_IDX(r, j) += pow(fabs(rel), 2.0f);
			}
			V_IDX(r, j) = sqrt(V_IDX(r, j));
		}
	}

	/* cleaning up */
	free(S);
	free(Y0);
	free(A_y0);
	free(Y);
	free(Ys);
	free(X);
}
