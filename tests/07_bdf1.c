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
#include "tsttools.h"

#include <math.h>

static const float a = 4.0f;
static const float c = 1.0f;

void f_volt(float t, float *x, float *y)
{
	float x1 = V_IDX(x, 1);
	float x2 = V_IDX(x, 2);
	V_IDX(y, 1) = a * (x1 - x1 * x2);
	V_IDX(y, 2) = -c * (x2 - x1 * x2);
}

void df_volt(float t, float *x, float *J)
{
	float x1 = V_IDX(x, 1);
	float x2 = V_IDX(x, 2);

	M_IDX(J, 2, 1, 1) = a * (1.0f - x2);
	M_IDX(J, 2, 2, 1) = c * x2;
	M_IDX(J, 2, 1, 2) = -a * x1;
	M_IDX(J, 2, 2, 2) = -c * (1.0f - x1);
}

int main(void)
{
	float y0[2] = { 2.0f, 1.0f };
	bdf1(2, f_volt, df_volt, 0.0f, y0, 10.0f, 1000, 0.0001f);
	printf("BDF1: %f, %f\n", y0[0], y0[1]);

	//assert(fabs(y0[0] + 1.0f) < 0.01f);
	//assert(fabs(y0[1]) < 0.01f);
	return 0;
}

