#ifndef MODEL_H
#define MODEL_H

#include <iostream>
#include <cmath>
#include <cassert>
#include "/Users/maboudi/c-programs/libraries/math/matrix.h"

using namespace std;

class Model
{
	private:
		double dt;
		double cut_off;
		Matrix region;
		Matrix num_mols_dir;

		int num_mols;
		double density;
		double temperature;
		double vel_mag;

		Matrix z;
		int num_chains;
		int chain_length;
		double bond_lim;
		Matrix initUchain;
		Matrix initUcell;
		double mol_dist;

		int vtk_count;

		double potential;

		Matrix q;
		Matrix p;
		Matrix f;
		
		Matrix snaps;

		Matrix A;
	public:
		Model();

		void initial_condition(Matrix*,Matrix*);

		void leap_frog(Matrix,Matrix,int);

		void deriv(Matrix*,int);
		void compute_forces();
		void compute_bonded_forces();

		void boundary_conditions(Matrix*);
		void wrap(Matrix*);

		void simulate();
		void generate_snapshots();
		void generate_reduced_basis();
		
		void evaluate_parameters();

		void export_vtk(Matrix*);
};

#endif
