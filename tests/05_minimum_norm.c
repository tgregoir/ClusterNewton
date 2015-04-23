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

