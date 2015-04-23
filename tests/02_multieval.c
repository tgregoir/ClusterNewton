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
