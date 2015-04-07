#include <stdio.h>
#include <stdlib.h>
#include "particle.h"
#include "gpukernels.h"

//__global__ void kernel_function(Particle* particles, int num_particles)
//{
//	int idx = threadIdx.x + blockIdx.x*blockDim.x;
//	if(idx < num_particles)
//		particles[idx].randomize();
//}

int main()
{
	int num_particles = 2*512;
	Particle* particle_array = new Particle[ num_particles ];
	Particle* device_parray = NULL;
	cudaMalloc(&device_parray, num_particles*sizeof(Particle));
	
	dim3 grid_size;
	grid_size.x = 2;

	dim3 block_size;
	block_size.x = 512;

	kernel_function<<<grid_size,block_size>>>(device_parray,num_particles);

	cudaMemcpy(particle_array, device_parray, num_particles*sizeof(Particle), cudaMemcpyDeviceToHost);

	for(int i = 0 ; i < num_particles ; i++)
		particle_array[i].print();
	return 0;
}
