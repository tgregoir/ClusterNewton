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
static float x[14];

/** F_HIV() - Forward problem for HIV Kinetics model
 * @u:                Vector of size 4.
 * @d:                Output, vector of size 4.
 *
 * The HIV Kinetics model (Miao et al.) is given by the differential
 * system: u' = F_HIV(t, u).
 */
static void F_HIV(float t, float *u, float *d)
{
	float u1 = V_IDX(u, 1);
	float u2 = V_IDX(u, 2);
	float u3 = V_IDX(u, 3);
	float u4 = V_IDX(u, 4);

	V_IDX(d, 1) = (x[1] - x[5] * u2 - x[6] * u3 - x[7] * u4) * u1;
	V_IDX(d, 2) = (x[2] + x[5] * u1 - x[8] * u3) * u2
	               + x[7] * u4 * u1 / 4.0f;
	V_IDX(d, 3) = (x[3] + x[6] * u1 - x[9] * u2) * u3
	               + x[7] * u4 * u1 / 4.0f;
	V_IDX(d, 4) = (x[4] + x[7] * u1 / 2.0f) * u4
	               + (x[8] + x[9]) * u3 * u2;
}

/* Times at which experimental data have been gathered. */
static const float tf[5] = {
	70.0f, 95.0f, 115.0f, 140.0f, 160.0f
};

/* Experimental data */
static const float target[5] = {
	0.0f
};

/* Number of steps in the RK4 scheme */
static const uint N[22] = {
	20, 20, 20, 20, 20, 20, 20, 20,
	20, 20, 20, 20, 20, 20, 20, 20,
	20, 20, 20, 20, 20, 20
};

void fwd_HIV(float *X, float *Y)
{
	/* set up the parameters so that F_HIV() can access them */
	for (uint i = 1; i <= 13; i++) {
		x[i] = V_IDX(X, i);
	}

	/* for each possible final time, simulate the system */
	float u[4];
	for (uint i = 1; i <= 5; i++) {
		V_IDX(u, 1) = x[10];
		V_IDX(u, 2) = x[11];
		V_IDX(u, 3) = x[12];
		V_IDX(u, 4) = x[13];
		rk4(4, F_HIV, 0.0f, u, V_IDX(tf, i), V_IDX(N, i));
		/* FIXME: is u4 the quantity of interest here? */
		V_IDX(Y, i) = V_IDX(u, 4);
	}
}

void hiv(void)
{
	float X[13] = {
		0.3f, 1.2f, 0.7f, 3.3f, 0.4f, 0.7f, 1.1f,
		0.5f, 3.3f, 0.2f, 1.4f, 0.1f, 3.7f
	};
	float Y[5] = { 0.0f };
	fwd_HIV(X, Y);
	print_vector(5, Y);
}
