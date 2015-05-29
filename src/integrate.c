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
