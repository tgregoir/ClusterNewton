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

/*
void f_influenza_kinetics(float *in, float *out)
{
	float x1 = in[0], x2 = in[1], x3 = in[2];
	float x4 = in[3], x5 = in[4], x6 = in[5];
	float x7 = in[6];

	uint Nsteps = 180;
	float deltaT = 180.0f / Nsteps;

	float *u1 = create_vector(Nsteps);
	float *u2 = create_vector(Nsteps);
	float *u3 = create_vector(Nsteps);
	float *u4 = create_vector(Nsteps);

	V_IDX(u1, 1) = x5;
	V_IDX(u2, 1) = 0.0f;
	V_IDX(u3, 1) = 0.0f;
	V_IDX(u4, 1) = x7;

	free(u4);
	free(u3);
	free(u2);
	free(u1);
}
*/

void f(float *in, float *out)
{
	float x1 = V_IDX(in, 1);
	float x2 = V_IDX(in, 2);
	V_IDX(out, 1) = (x1 * x1 + x2 * x2);
	V_IDX(out, 1) += sin(10000. * x1) * sin(10000. * x2) / 100.;
	//usleep(370000);
	//V_IDX(out, 1) = x1 + x2;
}

/* Y = y y'
 * Y' = y' -y */
void f_cos(float t, float *y, float *F)
{
	F[0] = y[1];
	F[1] = -y[0];
}

int main(void)
{
	float y0[2] = { 1.0f, 0.0f };
	rk4(1, 1, f_cos, 0.0f, y0, 3.14157f, 20);
	printf("RK4: %f\n", y0[0]);

	//srand(324635343);
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
