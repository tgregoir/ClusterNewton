#ifndef CN_H
#define CN_H

#include "common.h"

void random_pts_in_box(uint, uint, float *, float *, float *);
void perturbate(uint, uint, float *, float, float *);
void multi_eval(uint, uint, void (*)(float *, float *),
                uint, float *, float*);
void least_squares(uint, uint, float *, uint, float *, float *);

void cluster_newton(uint, uint, void (*)(float *, float *), float *,
                    float *, float *, uint, float, uint);

#endif /* CN_H */
