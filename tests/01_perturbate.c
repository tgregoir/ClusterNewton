#include "cn.h"

/**
 * perturbate_test() - a simple test for perturbate()
 *
 * Creates a vector, perturbates it and checks that it satisfies the bound
 * given in the documentation of perturbate().
 */
int main(void)
{
	uint l = 10;
	float yv = 5.;
	float eta = 0.1;

	float *ys = create_vector(l);

	for (uint i = 1; i <= l; i++) {
		V_IDX(ys, i) = (float)rand();
	}

	perturbate(l, yv, eta, ys);

	for (uint i = 1; i <= l; i++) {
		assert(abs((V_IDX(ys, i) - yv) / yv) <= eta);
	}

	free(ys);

	return 0;
}
