#include "cn.h"
#include "tsttools.h"
#include <math.h>

int main(void)
{
	uint m = 2;
	uint n = 3;
	uint l = 1;
	float A[] = {
		1.,  -1.,
		1., 1.,
		2., 1.
	};
	float B[] = {
		2., 4., 8.
	};

	float *X = create_matrix(l, m);
	least_squares(m, n, A, l, B, X);
	print_matrix(l, m, X);
	//assert(fabs(M_IDX(X, l, 1, 1) - 30. / 13. < 0.00001f));

	free(X);

	return 0;
}

