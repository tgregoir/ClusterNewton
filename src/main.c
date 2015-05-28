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
#include "integrate.h"

#include <math.h>
#include <unistd.h>

void f(float *in, float *out)
{
	float x1 = V_IDX(in, 1);
	float x2 = V_IDX(in, 2);
	V_IDX(out, 1) = (x1 * x1 + x2 * x2);
	V_IDX(out, 1) += sin(10000. * x1) * sin(10000. * x2) / 100.;
	//usleep(370000);
}

int main(void)
{
	srand(1429874166);
	//srand(time(NULL));

	uint m = 2;
	uint n = 1;
	uint l = 10;
	float ys[1] = { 100.0f };
	float xh[2] = { 2.5f, 2.5f };
	float v[2] = { 1.0f, 1.0f };
	float eta = 0.01f;
	uint K = 10;

	float *X = create_matrix(m, l);
	float *r = create_vector(l);
	//printf("m=%u, n=%u, l=%u, K=%u\n", m, n, l, K);
	cluster_newton(m, n, f, ys, xh, v, l, eta, K, X, r);

	print_vector(l, r);
	//print_matrix(m, l, X);
	//printf("plot(X(1,:),X(2,:), '.'), axis equal, xlim([-15 15]), ylim([-15 15])\n");
	//printf("t = linspace(0, 2*pi);\n");
	//printf("hold on, plot(10. * cos(t), 10. * sin(t))\n");

	free(r);
	free(X);

	return 0;
}
