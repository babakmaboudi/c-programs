#ifndef PARTICLE_H
#define PARTICLE_H

#include <math.h>
#include <iostream>
#include <vector>

using namespace std;

class Particle
{
	private:
		float x;
		float y;
		float z;
		vector<float> coords;
	public:
		__host__ __device__ void randomize();
		void print();
};

#endif
