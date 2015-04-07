#include "particle.h"

__host__ __device__ void Particle::randomize()
{
	x = 1;
	y = 1;
	z = 1;
}

void Particle::print()
{
	cout << x << " " << y << " " << z << endl;
}
