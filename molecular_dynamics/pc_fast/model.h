#ifndef MODEL_H
#define MODEL_H

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include "vec.h"
#include "matrix.h"
#include "site.h"
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
		Vec body_inert;
		
		// particle variables
		vector<Mol> mols;

		// site variables
		int sites_mol;			// number of interaction sites in each molecule
		vector< Site > ref_sites;
		vector< vector<Site> > all_sites;	// site of all molecules

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

		// initialization
		void init_coords();
		void init_vels();
		void init_accels();
		void initiate_cells();
		void initialize();

		// force calculation
		void boundary_conditions();
		void compute_forces();
		void compute_forces_cell();

		// moment calculation
		void compute_quat_accels();
		void compute_ang_vel(int,Vec*);
		void compute_torqs();
		void compute_rotation_mat(Matrix*,Quat*,int);
		void compute_site_coords();

		// time stepping
		void leap_frog(int);
		void predictor_step();
		void corrector_step();
		void pc(Mol*,real*,real*,int);
		void single_step();
		void simulate();

		// parameter estimation
		void correct_temperature();
		void evaluate_parameters();

		// boundary calculation
		void wrap(Vec*);
		void wrap_cell(Vec_int*,Vec*);

		// print
		void export_file(vector< vector<Mol> >*);
		void print();
};

#endif
