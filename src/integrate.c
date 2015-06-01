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
 * rk4() - Fourth-order Runge-Kutta method
 * @n:         A positive integer.
 * @f:         f : R x R^n -> R^n.
 * @t0:        Initial time.
 * @y:         Vector of size n. Input: y(t0). Output: y(t1).
 * @t1:        Final time.
 * @N:         Number of steps.
 *
 * Computes y(t1), where y' = f(t,y) and y(t0) = y0.
 */
void rk4(uint n, void (*f)(float, float *, float *),
         float t0, float *y, float t1, uint N)
{
	assert(t0 < t1);

	/* allocate memory */
	float *k1 = create_vector(n);
	float *k2 = create_vector(n);
	float *k3 = create_vector(n);
	float *k4 = create_vector(n);
	float *z = create_vector(n);

	float h = (t1 - t0) / N;
	float t = t0;

	for (uint i = 1; i <= N; i++) {
		/* k1 = f(t,y) */
		f(t, y, k1);

		/* k2 = f(t + h/2, y + h/2 k1) */
		for (uint j = 1; j <= n; j++) {
			V_IDX(z, j) = V_IDX(y, j) + 0.5f * h * V_IDX(k1, j);
		}
		f(t + 0.5f * h, z, k2);

		/* k2 = f(t + h/2, y + h/2 k2) */
		for (uint j = 1; j <= n; j++) {
			V_IDX(z, j) = V_IDX(y, j) + 0.5f * h * V_IDX(k2, j);
		}
		f(t + 0.5f * h, z, k3);

		/* k2 = f(t + h, y + h k3) */
		for (uint j = 1; j <= n; j++) {
			V_IDX(z, j) = V_IDX(y, j) + h * V_IDX(k3, j);
		}
		f(t + h, z, k4);

		/* y <- y + h/6 (k1 + 2k2 + 2k3 + k4) */
		for (uint j = 1; j <= n; j++) {
			V_IDX(y, j) += (V_IDX(k1, j) + 2.0f * V_IDX(k2, j)
			                + 2.0f * V_IDX(k3, j)
			                + V_IDX(k4, j)) * h / 6.0f;
		}
		t += h;
	}

	/* clean up */
	free(z);
	free(k4);
	free(k3);
	free(k2);
	free(k1);
}

/**
 * bdf1() - Backwards Euler method
 * @n:          Dimension of the problem.
 * @f:          f : R x R^n -> R^n.
 * @df:         The Jacobian of f.
 * @t0:         Initial time.
 * @y:          Vector of size n. Input: y(t0). Output: y(t1).
 * @t1:         Final time.
 * @N:          Number of steps.
 * @tol:        Tolerance internally used in Newton's method.
 *
 * Computes y(t1), where y' = f(t,y) and y(t0) = y0.
 */
void bdf1(uint n, void (*f)(float, float *, float *),
	  void (*df)(float, float *, float *),
          float t0, float *y, float t1, uint N, float tol)
{
	assert(t0 < t1);

	/* allocate memory */
	float *x = create_vector(n);
	float *z = create_vector(n);
	float *D = create_matrix(n, n);
	int *ipiv = (int*)malloc(sizeof(int) * n);
	assert(ipiv);

	float h = (t1 - t0) / (float)N;
	float t = t0;

	for (uint i = 1; i <= N; i++) {
		/* we're going to solve y_{i+1} = y_i + hf(t_{i+1}, y_{i+1})
		 * for y_{i+1} */
		t += h;

		/* F(t,x) = x - y - hf(t,x)
		 * J(t,x) = 1 - hDf(t,x)
		 *
		 * Newton iteration:
		 *   x0 = y
		 *   x_{i+1} = x_i + J(t,x_i)^(-1)F(t,x_i)
		 */

		/* x <- y */
		m_copy(n, 1, n, y, n, x);

		do {
			/* D = 1 - hJ(t,x) */
			df(t, x, D);
			m_scale(n, n, n, D, -h);
			for (uint i = 1; i <= n; i++) {
				M_IDX(D, n, i, i) += 1.0f;
			}

			/* z = x - y - hf(t,y) */
			f(t, y, z);
			for (uint i = 1; i <= n; i++) {
				V_IDX(z, i) = V_IDX(x, i) - V_IDX(y, i)
				              - h * V_IDX(z, i);
			}

			/* Replace z with D^(-1)z */
			LAPACKE_sgesv(LAPACK_COL_MAJOR, n, 1, D, n, ipiv,
			              z, n);

			/* x += z */
			m_add(n, 1, n, z, n, x);
		} while (v_norm(n, z) > tol);

		/* y <- x */
		m_copy(n, 1, n, x, n, y);
	}

	free(D);
	free(z);
	free(x);
}
