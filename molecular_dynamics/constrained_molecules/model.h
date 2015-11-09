#ifndef MODEL_H
#define MODEL_H

#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include "vec.h"
#include "matrix.h"
#include "mol.h"
#include "cons.h"

class Model
{
	private:
		// physical parameters
		int NDIM;		// number of dimensions
		real dt;		// time step
		real cut_off;		// cut off radius
		Vec region;		// size of region
		int num_mols;
		real density;
		real temperature;
		real vel_mag;		// initial velocity magnitude

		// particle varaibles
		vector<Mol> mols;
		int num_chains;
		int chain_length;
		real bond_length;
		real bond_angle;
		real bond_lim;
		Vec_int initUchain;
		Vec_int initUcell;

		// state parameters
		real potential;
		real virSum;

		// constraint parameters
		int num_constraints;
		vector<Cons> cons;
		Matrix m_mat;
		real* cons_vec;

		int vtk_count;
	public:
		Model();

		void init_coords_chain();
		void assign_to_chain();
		void init_vels();
		void init_accels();
		void build_constraint_matrix();
		void initialize();

		void compute_forces();
		void compute_chain_torsion_forces();
		void compute_chain_angel_forces();
		void compute_constraints();

		void predictor_step();
		void corrector_step();
		void pc(Mol*,real*,real*,int);
		void single_step();
		void simulate();

		void evaluate_parameters();

		void boundary_conditions();
		void wrap(Vec* v);

		void lin_eq(Matrix,real*,int);

		void export_vtk();
};

using namespace std;

#endif
