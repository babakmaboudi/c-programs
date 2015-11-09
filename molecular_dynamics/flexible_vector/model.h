#ifndef MODEL_H
#define MODEL_H

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <cassert>
#include "vec.h"
#include "site.h"
#include "mol.h"
#include "../../libraries/math/matrix.h"

using namespace std;

class Model
{
	private:
		// physical parameters
		int NDIM;			// number of dimensions
		double dt;			// time step
		double cut_off;			// cut off radius
		Vec region;			// size of region
		Vec num_mols_dir;		// number of particles in each direction
		int num_mols;
		double density;
		double temperature;
		double vel_mag;			// initial velocity magnitude
		Vec mol_inert;
		double bCon;			// electric charge density constant
		
		// particle variables
		vector<Mol> mols;
		int num_chains;
		int chain_length;
		double bond_lim;
		Vec_int initUchain;		// no clue what this is
		Vec_int initUcell;

		// site variables
		int num_sites;			// number of interaction sites in each molecule
		vector< Site > ref_sites;
		vector< vector<Site> > all_sites;	// site of all molecules

		// parameter estimation parameters
		double potential;
		double virSum;			// no clue what this is!!! probably f.r


		// cell variables
		vector< vector<int> > cell_list;	// list of particles in each cell
		Vec_int num_cells_dir;
		int num_cells;
		Vec cell_width;				// physical width of cell in each direction
		vector<Vec_int> neigs;			// vector pointing to neighbors

		int vtk_count;
	public:
		Model();

		// initial condition & initialization
		void init_coords();
		void init_coords_chain();
		void assign_to_chain();
		void init_vels();
		void init_accels();
		void initiate_cells();
		void define_mol();
		void initialize();

		// force calculation
		void boundary_conditions();
		void compute_forces();
		void compute_bonded_forces();

		// time stepping
		void leap_frog(int);
		void predictor_step();
		void corrector_step();
		void pc(Mol*,double*,double*,int);
		void single_step();
		void simulate();

		// parameter estimation
		void evaluate_parameters();

		// boundary calculation
		void wrap(Vec*);
		void wrap_cell(Vec_int*,Vec*);
		void compute_shift(Vec*,Vec*);

		// print
		void export_file(vector< vector<Mol> >*);
		void export_vtk();
		void print();
};

#endif
