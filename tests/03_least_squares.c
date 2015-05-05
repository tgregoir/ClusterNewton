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
#include <math.h>

int main(void)
{
	uint m = 1;
	uint n = 2;
	uint l = 1;
	float *A = create_matrix(m, n);
	float *B = create_matrix(l, n);

	M_IDX(A, m, 1, 1) = 2;
	M_IDX(A, m, 1, 2) = 3;
	M_IDX(B, l, 1, 1) = 6;
	M_IDX(B, l, 1, 2) = 6;

	float *X = create_matrix(l, m);
	normal_ls(m, n, A, l, B, X);
	assert(fabs(M_IDX(X, l, 1, 1) - 30. / 13. < 0.00001f));

	free(X);
	free(B);
	free(A);

	return 0;
}
