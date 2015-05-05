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
#ifndef CN_H
#define CN_H

#include "common.h"

void random_pts_in_box(uint, uint, float *, float *, float *);
void perturbate(uint, uint, float *, float, float *);
void multi_eval(uint, uint, void (*)(float *, float *),
                uint, float *, float*);
void pinv_ls(uint, uint, float *, uint, float *, float *);
void normal_ls(uint, uint, float *, uint, float *, float *);
void minimum_norm(uint, uint, float *, uint, float *, float *);

void cluster_newton(uint, uint, void (*)(float *, float *), float *,
                    float *, float *, uint, float, uint, float *, float *);

#endif /* CN_H */
