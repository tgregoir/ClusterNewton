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

static float x[8];

/** F_influenza() - Forward problem for Influenza Kinetics model
 * @u:                Vector of size 4.
 * @d:                Output, vector of size 4.
 *
 * The Influenza Kinetics model (Baccam et al.) is given by the differential
 * system: u' = F_influenza(t, u).
 */
static void F_influenza(float t, float *u, float *d)
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

/* Times at which experimental data have been gathered. */
const float tf[22] = {
	4.5f,   12.0f,  20.0f,  27.0f,  35.0f,  43.0f,  51.0f,  58.0f,
	66.0f,  73.0f,  81.0f,  89.0f,  96.0f,  105.0f, 112.0f, 120.0f,
	127.0f, 135.0f, 143.0f, 151.0f, 159.0f, 166.0f
};

/* Experimental data */
const float target[22] = {
	1.21f, 1.40f, 2.83f, 4.13f, 4.42f, 5.04f, 4.62f, 4.44f,
	5.13f, 4.33f, 4.02f, 3.13f, 2.90f, 3.0f,  3.02f, 3.12f,
	1.0f,  2.10f, 1.12f, 0.79f, 0.17f, 0.19f
};

const uint N[22] = { 20 };

void fwd_influenza(float *X, float *Y)
{
	/* set up the parameters so that F_influenza() can access them */
	for (uint i = 1; i <= 7; i++) {
		x[i] = V_IDX(X, i);
	}

	/* for each possible final time, simulate the system */
	for (uint i = 1; i <= 22; i++) {
		float u[4] = { x[5], 0.0f, 0.0f, x[7] };
		rk4(4, 4, F_influenza, 0.0f, u, V_IDX(tf, i), V_IDX(N, i));
		V_IDX(Y, i) = V_IDX(u, 4);
	}
}

void influenza(void)
{
}
