#ifndef MODEL_H
#define MODEL_H

#include <iostream>
#include "../libraries/math/matrix.h"

class Model
{
	private:
		int NDIM;			// number of dimensions
		double dt;			// time step
		double cut_off;			// cut off radius
		Matrix region;			// size of region
		Matrix num_mols_dir;		// number of particles in each direction
		int num_mols;
		double density;
		double temperature;
		double vel_mag;			// initial velocity magnitude

		vector<Mol> mols;
	public:
		Model();	
		void init_coords();
};

#endif
