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

/* For the sake of clarity, we'll ignore x[0]. */
static float x[8];

/** F_influenza() - Forward problem for Influenza Kinetics model
 * @t:                Time (unused here).
 * @u:                Vector of size 4.
 * @d:                Output, vector of size 4.
 *
 * The Influenza Kinetics model (Baccam et al.) is given by the differential
 * system: u' = F_influenza(t, u).
 */
void F_influenza(float t, float *u, float *d)
{
	float u1 = V_IDX(u, 1);
	float u2 = V_IDX(u, 2);
	float u3 = V_IDX(u, 3);
	float u4 = V_IDX(u, 4);

	V_IDX(d, 1) = -x[1] * u1 * u4;
	V_IDX(d, 2) = x[1] * u1 * u4 - u2 / x[2];
	V_IDX(d, 3) = u2 / x[2] - u3 / x[3];
	V_IDX(d, 4) = x[4] * u3 / x[5] - x[6] * u4;
}

/** dF_influenza() - Jacobian of F_influenza()
 * @t:                Time (unused here).
 * @u:                Vector of size 4.
 * @J:                Output, 4-by-4 matrix.
 */
void dF_influenza(float t, float *u, float *J)
{
	float u1 = V_IDX(u, 1);
	float u4 = V_IDX(u, 4);

	M_IDX(J, 4, 1, 1) = -x[1] * u4;
	M_IDX(J, 4, 2, 1) = x[1] * u4;
	M_IDX(J, 4, 3, 1) = 0.0f;
	M_IDX(J, 4, 4, 1) = 0.0f;

	M_IDX(J, 4, 1, 2) = 0.0f;
	M_IDX(J, 4, 2, 2) = -1.0f / x[2];
	M_IDX(J, 4, 3, 2) = 1.0f / x[2];
	M_IDX(J, 4, 4, 2) = 0.0f;

	M_IDX(J, 4, 1, 3) = 0.0f;
	M_IDX(J, 4, 2, 3) = 0.0f;
	M_IDX(J, 4, 3, 3) = -1.0f / x[3];
	M_IDX(J, 4, 4, 3) = x[4] / x[5];

	M_IDX(J, 4, 1, 4) = -x[1] * u1;
	M_IDX(J, 4, 2, 4) = x[1] * u1;
	M_IDX(J, 4, 3, 4) = 0.0f;
	M_IDX(J, 4, 4, 4) = -x[6];
}

/* Times at which experimental data have been gathered. */
static const float tf[22] = {
	4.5f,   12.0f,  20.0f,  27.0f,  35.0f,  43.0f,  51.0f,  58.0f,
	66.0f,  73.0f,  81.0f,  89.0f,  96.0f,  105.0f, 112.0f, 120.0f,
	127.0f, 135.0f, 143.0f, 151.0f, 159.0f, 166.0f
};

/* Experimental data */
static const float target[22] = {
	1.21f, 1.40f, 2.83f, 4.13f, 4.42f, 5.04f, 4.62f, 4.44f,
	5.13f, 4.33f, 4.02f, 3.13f, 2.90f, 3.0f,  3.02f, 3.12f,
	1.0f,  2.10f, 1.12f, 0.79f, 0.17f, 0.19f
};

/* Number of steps in the RK4 scheme */
static const uint N[22] = {
	200, 200, 200, 200, 200, 200, 200, 200,
	200, 200, 200, 200, 200, 200, 200, 200,
	200, 200, 200, 200, 200, 200
};

void fwd_influenza(float *X, float *Y)
{
	/* set up the parameters so that F_influenza() can access them */
	for (uint i = 1; i <= 7; i++) {
		x[i] = V_IDX(X, i);
	}

	/* for each possible final time, simulate the system */
	float u[4];
	for (uint i = 1; i <= 22; i++) {
		V_IDX(u, 1) = x[5];
		V_IDX(u, 2) = 0.0f;
		V_IDX(u, 3) = 0.0f;
		V_IDX(u, 4) = x[7];
		//rk4(4, F_influenza, 0.0f, u, V_IDX(tf, i), V_IDX(N, i));
		bdf1(4, F_influenza, dF_influenza,
		     0.0f, u, V_IDX(tf, i), V_IDX(N, i), 0.001);
		V_IDX(Y, i) = V_IDX(u, 4);
	}
}

void influenza(void)
{
	float X[7] = { 0.3f, 1.2f, 0.7f, 3.3f, 0.4f, 0.7f, 1.1f };
	float Y[22] = { 0.0f };
	fwd_influenza(X, Y);
	print_vector(22, Y);
}
