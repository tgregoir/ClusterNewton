#ifndef CN_H
#define CN_H

#include "common.h"

void choose_init_pts(uint, uint, float *, float *, float *);
void perturbate(uint, float, float, float *);
void multi_eval(uint, float (*)(float *), uint, float *, float*);

void cluster_newton(uint, uint, float *, float *, float, float, uint,
                    float (*)(float *));

#endif /* CN_H */
