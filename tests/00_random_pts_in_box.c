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
	init_prg();

	uint m = random_dim();
	uint l = random_dim();
	//printf("m=%u l=%u\n",m,l);

	float *xh = random_vector(m);
	float *v = random_vector(m);
	float *X = create_matrix(m + 1, l);

	random_pts_in_box(m, l, xh, v, X);

	for (uint k = 1; k <= l; k++) {
		for (uint i = 1; i <= m; i++) {
			assert (abs(M_IDX(X, m + 1, i, k) - V_IDX(xh, i))
			        <= abs(V_IDX(xh, i)) * abs(V_IDX(v, i)));
		}
	}

	free(X);
	free(v);
	free(xh);

	return 0;
}
