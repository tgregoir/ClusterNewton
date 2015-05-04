#include "cn.h"
#include "tsttools.h"
#include <math.h>

int main(void)
{
	uint m = 1;
	uint n = 2;
	uint l = 1;
	float *A = create_matrix(m, n);
	float *B = create_matrix(l, n);

	M_IDX(A, m, 1, 1) = 2;
	M_IDX(A, m, 1, 2) = 3;
	M_IDX(B, l, 1, 1) = 6;
	M_IDX(B, l, 1, 2) = 6;

	float *X = create_matrix(l, m);
	normal_ls(m, n, A, l, B, X);
	assert(fabs(M_IDX(X, l, 1, 1) - 30. / 13. < 0.00001f));

	free(X);
	free(B);
	free(A);

	return 0;
}
