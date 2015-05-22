/*
 *    This file is part of CNewt.
 *
 *    CNewt is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    CNewt is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with CNewt.  If not, see <http://www.gnu.org/licenses/>.
 */
extern "C" {
#include "common.h"
}

void multi_eval_sequential(uint m, uint n, void (*f)(float *, float *),
                           uint l, float *X, float *Y)
{
	for (uint j = 1; j <= l; j++) {
		f(M_COL(X, m + 1, j), M_COL(Y, n, j));
	}
}

__global__ void eval_fct_kernel(const float *X, float *Y,
                                uint m, uint n, uint l)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;

	if (i <= l) {
		const float *in = &X[(m + 1) * i];
		float *out = &Y[n * i];
		float x1 = V_IDX(in, 1);
		float x2 = V_IDX(in, 2);
		V_IDX(out, 1) = (x1 * x1 + x2 * x2);
		V_IDX(out, 1) += sin(10000.f * x1) * sin(10000.f * x2) / 100.f;
	}
}

void multi_eval_gpu(uint m, uint n, uint l, float *X, float *Y)
{
	// Load X to device memory
	uint sizeX = (m + 1) * l * sizeof(float);
	float *devX = NULL;
	cudaMalloc(&devX, sizeX);	
	cudaMemcpy(devX, X, sizeX, cudaMemcpyHostToDevice);

	// Allocate Y in device memory
	uint sizeY = n * l * sizeof(float);
	float *devY = NULL;
	cudaMalloc(&devY, sizeY);

	// Invoke kernel
	int threadsPerBlock = 256;
	int blocksPerGrid = (l + threadsPerBlock - 1) / threadsPerBlock;
	eval_fct_kernel<<<blocksPerGrid, threadsPerBlock>>>(devX, devY, m, n, l);

	// Read Y from device memory
	cudaMemcpy(Y, devY, sizeY, cudaMemcpyDeviceToHost);

	// Free device memory
	cudaFree(devX);
	cudaFree(devY);
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
extern "C" void multi_eval(uint m, uint n, void (*f)(float *, float *),
                           uint l, float *X, float *Y)
{
	//multi_eval_sequential(m, n, f, l, X, Y);
	multi_eval_gpu(m, n, l, X, Y);
}
