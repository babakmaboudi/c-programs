#include <iostream>
#include <cstdlib>
#include "particle.h"
#include "hostfunctions.h"
#include "cell.h"

using namespace std;

int main()
{
	int num_particles = 1000;
	int num_cells = 5;

	Particle* particles;
	particles = initiate(num_particles);

	Cell* domain;
	int* neighbour_cells;
	domain = decompose_domain(num_particles,num_cells,&neighbour_cells);

	for(int i=0 ; i<num_cells*num_cells ; i++)
	{
		cout << i << " : ";
		for(int j=0 ; j<8 ; j++)
			cout << neighbour_cells[i*8+j] << " ";
		cout << endl;
	}

//	Particle* device_particles;
//	Cell* device_domain;
//	int** device_neighbours;
//
//	cudaMalloc(&device_particles, num_particles*sizeof(Particle));
//	cudaMemcpy(device_particles, particles, num_particles*sizeof(Particle), cudaMemcpyHostToDevice);
//
//	cudaMalloc(&device_domain, num_cells*num_cells*sizeof(Cell));
//	cudaMemcpy(device_domain, domain, num_cells*num_cells*sizeof(Cell), cudaMemcpyHostToDevice);
//
//	cudaMalloc(&device_neighbours, num_cells*num_cells*sizeof(int*));
//	for(int i=0; i<num_cells*num_cells ; i++)
//	{
//		cudaMalloc(&device_neighbours[i], 8*sizeof(int));
//		cudaMemcpy(device_neighbours[i], neighbour_cells[i], 8*sizeof(int), cudaMemcpyHostToDevice);
//	}
//
//	Particle* test = new Particle[num_particles];
//	cudaMemcpy(&test, device_particles, num_particles*sizeof(Particle), cudaMemcpyHostToDevice);
//
//	for(int i=0 ; i<num_particles ; i++)
//	{
//		cout << particles[i].pos_x << endl;
//	}

	return 0;
}
