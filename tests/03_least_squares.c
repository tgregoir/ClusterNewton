#include "cn.h"
#include "tsttools.h"
#include <math.h>

int main(void)
{
	uint m = 2;
	uint n = 1;
	uint l = 1;
	float *A = create_matrix(m, n);
	float *B = create_matrix(m, l);

	M_IDX(A, m, 1, 1) = 2;
	M_IDX(A, m, 2, 1) = 3;
	M_IDX(B, m, 1, 1) = 6;
	M_IDX(B, m, 2, 1) = 6;

	float *X = create_matrix(n, l);
	least_squares(m, n, A, l, B, X);
	assert(fabs(M_IDX(X, n, 1, 1) - 30. / 13. < 0.00001f));

	free(X);
	free(B);
	free(A);

	return 0;
}
