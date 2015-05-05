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

#include <math.h>

void f(float *in, float *out)
{
	float x1 = V_IDX(in, 1);
	float x2 = V_IDX(in, 2);
	V_IDX(out, 1) = (x1 * x1 + x2 * x2);
	V_IDX(out, 1) += sin(10000. * x1) * sin(10000. * x2) / 100.;
}

int main(void)
{
	srand(1429874166);

	uint m = 2;
	uint n = 1;
	uint l = 100;
	float ys[] = { 100. };
	float xh[] = { 2.5, 2.5 };
	float v[] = { 1., 1. };
	float eta = 0.1;
	uint K = 10;

	float *res = create_matrix(m, l);
	printf("m=%u, n=%u, l=%u\n", m, n, l);
	cluster_newton(m, n, f, ys, xh, v, l, eta, K, res);
	//print_matrix(m, l, res);
	for (uint j = 1; j <= l; j++) {
		float x = M_IDX(res, m, 1, j);
		float y = M_IDX(res, m, 2, j);
		printf("%f\n", sqrt(x * x + y * y));
	}

	free(res);

	return 0;
}
