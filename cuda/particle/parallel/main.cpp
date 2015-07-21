#include <iostream>
#include <cstdlib>
#include <ctime>
#include "particle.h"
#include "cell.h"
#include "hostfunctions.h"
#include "devicefunctions.h"

using namespace std;

int main()
{
	srand48( 0 );
	int num_particles = 1000;
	int num_cells = 5;
	int MAX_ITER = 1000;

	Particle* particles;
	particles = initiate(num_particles);

	Cell* domain;
	int* neighbour_cells;
	domain = decompose_domain(num_particles,num_cells,&neighbour_cells);
	
	int average = num_particles/(num_cells*num_cells);
	
	Particle* device_particles;
	Cell* device_domain;
	int* device_neighbours;
	int* device_tracker;

	cudaMalloc(&device_particles, num_particles*sizeof(Particle));
	cudaMemcpy(device_particles, particles, num_particles*sizeof(Particle), cudaMemcpyHostToDevice);

	cudaMalloc(&device_domain, num_cells*num_cells*sizeof(Cell));
	cudaMemcpy(device_domain, domain, num_cells*num_cells*sizeof(Cell), cudaMemcpyHostToDevice);

	cudaMalloc(&device_neighbours, num_cells*num_cells*8*sizeof(int));
	cudaMemcpy(device_neighbours, neighbour_cells, num_cells*num_cells*8*sizeof(int), cudaMemcpyHostToDevice);

	cudaMalloc(&device_tracker, num_cells*num_cells*average*2*sizeof(int));

	dim3 grid_size;
	grid_size.x = num_cells*num_cells;

	dim3 block_size;
	block_size.x = num_particles*2/grid_size.x;

	clock_t begin = clock();
	for(int i=0 ; i<MAX_ITER ; i++)
	{
		find_cell<<<grid_size,block_size>>>(device_particles,device_domain,num_particles,num_cells*num_cells);
		find_particles<<<grid_size,block_size>>>(device_particles,device_domain,device_tracker,num_particles,num_cells*num_cells,average*2);
		compute_force<<<grid_size,block_size>>>(device_particles,device_domain,device_tracker,num_particles,num_cells*num_cells,average*2,device_neighbours);
		move<<<grid_size,block_size>>>(device_particles,num_particles);
	}
	clock_t end = clock();
	double elapsed_time = double(end - begin) / CLOCKS_PER_SEC;
	cout << elapsed_time << endl;

	return 0;
}
