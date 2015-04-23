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
	float *X = create_matrix(m + l, l);

	random_pts_in_box(m, l, xh, v, X);

	for (uint k = 1; k <= l; k++) {
		for (uint i = 1; i <= m; i++) {
			assert (abs(M_IDX(X, m + l, i, k) - V_IDX(xh, i))
			        <= abs(V_IDX(xh, i)) * abs(V_IDX(v, i)));
		}
	}
	for (uint k = 1; k <= l; k++) {
		for (uint j = 1; j <= l; j++) {
			if (j == k)
				assert(M_IDX(X, m + l, m + k, j) == 1.);
			else
				assert(M_IDX(X, m + l, m + k, j) == 0.);
		}
	}

	free(X);
	free(v);
	free(xh);

	return 0;
}
