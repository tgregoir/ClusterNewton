#include "cn.h"

#include <math.h>

void f(float *in, float *out)
{
	float x1 = V_IDX(in, 1);
	float x2 = V_IDX(in, 2);
	V_IDX(out, 1) = (x1 * x1 + x2 * x2);
	V_IDX(out, 1) += sin(10000. * x1) * sin(10000. * x2) / 100.;
}

int main(void)
{
	srand(1429874166);

	uint m = 2;
	uint n = 1;
	uint l = 3;
	float ys[] = { 100. };
	float xh[] = { 2.5, 2.5 };
	float v[] = { 1., 1. };
	float eta = 0.1;
	uint K = 5;

	float *res = create_matrix(m, l);
	printf("m=%u, n=%u, l=%u\n", m, n, l);
	cluster_newton(m, n, f, ys, xh, v, l, eta, K, res);
	//print_matrix(m, l, res);
	for (uint j = 1; j <= l; j++) {
		float x = M_IDX(res, m, 1, j);
		float y = M_IDX(res, m, 2, j);
		printf("%f\n", sqrt(x * x + y * y));
	}

	free(res);

	return 0;
}
