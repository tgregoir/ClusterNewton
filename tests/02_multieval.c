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
#include "tsttools.h"

void f(float *in, float *out)
{
	out[0] = in[0] * in[1];
}

int main(void)
{
	init_prg();

	uint l = random_dim();
	float *X = random_matrix(2 + l, l);
	float *Y = random_vector(l);

	multi_eval(2, 1, f, l, X, Y);

	for (uint i = 1; i <= l; i++) {
		float x[2] = { M_IDX(X, 2 + l, 1, i), M_IDX(X, 2 + l, 2, i) };
		float y;
		f(x, &y);
		assert(y == V_IDX(Y, i));
	}

	free(Y);
	free(X);

	return 0;
}
