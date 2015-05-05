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

int main(void)
{
	uint m = 1;
	uint n = 2;
	uint l = 1;
	float A[] = {
		1., 1.
	};
	float B[] = {
		2.
	};

	float *X = create_matrix(n, l);
	minimum_norm(m, n, A, l, B, X);
	print_matrix(n, l, X);

	free(X);

	return 0;
}

