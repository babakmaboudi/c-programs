#ifndef DEVICEFUNCTIONS_H
#define DEVICEFUNCTIONS_H

#include <cmath>
#include "particle.h"
#include "cell.h"

#define TILTE_SIZE 2*80
#define	density 0.0005
#define	mass 0.1
#define	cutoff 0.01
#define	min_r cutoff/100
#define	dt 0.0005

__global__ void find_cell(Particle*,Cell*,int,int);
__global__ void find_particles(Particle*,Cell*,int*,int,int,int);
__global__ void compute_force(Particle*,Cell*,int*,int,int,int,int*);
__global__ void move(Particle*,int);

__device__ void apply_force(float*,float*,float,float,float,float);

#endif
