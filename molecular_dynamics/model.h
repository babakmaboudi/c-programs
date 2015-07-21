#ifndef MODEL_H
#define MODEL_H

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include "vec.h"
#include "mol.h"

using namespace std;

class Model
{
	private:
		int NDIM;			// number of dimensions
		double dt;			// time step
		double cut_off;			// cut off radius
		Vec region;			// size of region
		Vec num_mols_dir;		// number of particles in each direction
		int num_mols;
		double density;
		double temperature;
		double vel_mag;			// initial velocity magnitude

		vector<Mol> mols;
	public:
		Model();	
		void init_coords();
		void init_vels();
		void initialize();

		void boundary_conditions();
		void compute_forces();
		void leap_frog(int);

		void single_step();
		void simulate();
		void evaluate_parameters();

		void wrap(Vec*);

		void export_file(vector< vector<Mol> >*);
		void print();
};

#endif
