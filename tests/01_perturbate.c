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
