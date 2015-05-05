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

	uint l = random_dim();
	uint n = random_dim();
	float *ys = random_vector(n);
	float eta = 0.001;

	float *Ys = random_matrix(n, l);

	perturbate(l, n, ys, eta, Ys);

	for (uint j = 1; j <= l; j++) {
		for (uint i = 1; i <= n; i++) {
			assert(abs((M_IDX(Ys, n, i, j) - V_IDX(ys, i))
				   / V_IDX(ys, i)) <= eta);
		}
	}

	free(Ys);
	free(ys);

	return 0;
}
