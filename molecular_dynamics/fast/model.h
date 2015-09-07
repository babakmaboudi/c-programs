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
		// physical parameters
		int NDIM;			// number of dimensions
		real dt;			// time step
		real cut_off;			// cut off radius
		Vec region;			// size of region
		Vec num_mols_dir;		// number of particles in each direction
		int num_mols;
		real density;
		real temperature;
		real vel_mag;			// initial velocity magnitude
		
		// particle variables
		vector<Mol> mols;

		// parameter estimation parameters
		real potential;
		real virSum;			// no clue what this is!!! probably f.r


		// cell variables
		vector< vector<int> > cell_list;	// list of particles in each cell
		Vec_int num_cells_dir;
		int num_cells;
		Vec cell_width;				// physical width of cell in each direction
		vector<Vec_int> neigs;			// vector pointing to neighbors
	public:
		Model();	
		void init_coords();
		void init_vels();
		void initiate_cells();
		void initialize();

		void boundary_conditions();
		void compute_forces();
		void compute_forces_cell();
		void leap_frog(int);

		void single_step();
		void simulate();
		void correct_temperature();
		void evaluate_parameters();

		void wrap(Vec*);
		void wrap_cell(Vec_int*,Vec*);

		void export_file(vector< vector<Mol> >*);
		void print();
};

#endif
