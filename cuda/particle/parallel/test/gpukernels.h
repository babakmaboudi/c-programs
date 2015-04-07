#ifndef GPUKERNELS_H
#define GPUKERNELS_H

#include "particle.h"

__global__ void kernel_function(Particle* particles, int num_particles)
{
	int idx = threadIdx.x + blockIdx.x*blockDim.x;
	if(idx < num_particles)
		particles[idx].randomize();
}



#endif
