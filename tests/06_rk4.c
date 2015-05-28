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

void f_cos(float t, float *y, float *F)
{
	F[0] = y[1];
	F[1] = -y[0];
}

int main(void)
{
	float y0[2] = { 1.0f, 0.0f };
	rk4(1, 1, f_cos, 0.0f, y0, 3.14157f, 20);
	printf("RK4: %f, %f\n", y0[0], y0[1]);

	assert(fabs(y0[0] - 1.0f) < 0.01f);
	assert(fabs(y0[1]) < 0.01f);
	return 0;
}
