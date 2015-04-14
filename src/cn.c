#include "cn.h"

#include "common.h"
#include "linalg.h"

/**
 * choose_init_pts() - samples random points in a box
 * @d:      Dimension of the space.
 * @l:      Number of points to sample.
 * @xh:     A d-dimensional vector.
 * @v:      A d-dimensional vector.
 * @X:      The d-by-l matrix where the result is written.
 *
 * Samples l points uniformly at random in the box
 *   { x : abs((x(i) - xh(i)) / (xh(i) * v(i))) < 1 }.
 * x has dimension d.
 */
void choose_init_pts(uint d, uint l, float *xh, float *v, float *X)
{
	for (uint i = 1; i <= l; i++) {
		for (uint k = 1; k <= d; k++) {
			float r = 2. * (float)rand() / (float)RAND_MAX - 1.;
			M_IDX(X, d, k, i) = V_IDX(xh, k) *
			                    (1. + V_IDX(v, k) * r);
		}
	}
}

/**
 * perturbate() - creates a pertubated target vector
 * @l:             Size of the vector.
 * @yv:            Target value.
 * @eta:           Relative magnitude of the random perturbation.
 * @ys:            The l-dimensional vector where the result is written.
 *
 * Sets ys(i) = yv for all i, then introduces a random, uniform perturbation
 * such that
 *   abs((ys(j) - yv) / yv) < eta.
 */
void perturbate(uint l, float yv, float eta, float *ys)
{
	for (uint i = 1; i <= l; i++) {
		float r = eta * (2. * (float)rand() / (float)RAND_MAX - 1.);
		V_IDX(ys, i) = yv * (1. + r);
	}
}

/**
 * multi_eval() - evaluates a function at multiple points
 * @m:              Number of parameters of the function.
 * @f:              The function to evaluate.
 * @l:              Number of points.
 * @X:              Coordinates of the points, one point per column.
 * @y:              Vector in which to store the result.
 *
 * This function assumes COLUMN-MAJOR ORDER.
 */
void multi_eval(uint m, float (*f)(float *), uint l, float *X, float* y)
{
	for (uint k = 1; k <= l; k++) {
		V_IDX(y, k) = f(&M_IDX(X, m, 1, k));		
	}
}

void cluster_newton(uint m, uint l, float *xh, float *v, float eta, float yv,
                    uint K1, float (*f)(float *))
{
	// safety checks
	assert(m > 0);
	assert(l > 0);

	// 1.1.
	float *X = create_matrix(m, l);
	choose_init_pts(m, l, xh, v, X);

	// 1.2.
	float *ys = create_vector(l);
	perturbate(l, yv, eta, ys);

	float *y = create_vector(l);
	for (uint k = 0; k <= K1; k++) {
		// 2.1.
		multi_eval(m, f, l, X, y);
	}

	// cleaning up
	free(y);
	free(ys);
	free(X);
}
