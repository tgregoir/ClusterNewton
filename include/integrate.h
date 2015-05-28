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

#ifndef INTEGRATE_H
#define INTEGRATE_H

#include "common.h"

/**
 * rk4() - Fourth-order Runge-Kutta method
 * @n:         A positive integer.
 * @m:         A positive integer.
 * @f:         f : R x R^n -> R^m.
 * @t0:        Initial time.
 * @y:         Vector of size m. Input: y(t0). Output: y(t1).
 * @t1:        Final time.
 * @N:         Number of steps.
 *
 * Computes y(t1), where y' = f(t,y) and y(t0) = y0.
 */
void rk4(uint, uint, void (*)(float, float *, float *), float, float *,
         float, uint);

#endif /* INTEGRATE_H */
