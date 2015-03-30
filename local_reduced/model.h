#ifndef MODEL_H
#define MODEL_H

#include <iostream>
#include "parameters.h"
#include "tools.h"
#include "grid.h"
#include "../libraries/math/matrix.h"

using namespace std;

class Model
{
	private:
		Parameters param;
		Matrix Phi;

		Matrix cluster_centers;
		int nearest;
		vector<Matrix> local_phi;
	public:
		void sv2oe(vector<double>,vector<double>,double,vector<double>*);

		int find_closest_center(Matrix);
		void compute_significant_subspace(Matrix*,Matrix*);

		void set_parameters(double,Parameters*);
		void initiate_from_file(char*);

		double compute_hamiltonian(Matrix,Matrix,Parameters);

		void deriv(Matrix*,double,Matrix*,Parameters*);
		void deriv_symp(Matrix*,double,Matrix*,Parameters*,int);
		void deriv_symp_reduced(Matrix*,double,Matrix*,Parameters*,int);

		void integ4_symp(void (Model::*func)(Matrix*,double,Matrix*,Parameters*,int),double*,Matrix,Matrix,Parameters*,double,Matrix*,Matrix*,Matrix*,int);
		void integ4_symp_reduced(void (Model::*func)(Matrix*,double,Matrix*,Parameters*,int),double*,Matrix,Matrix,Parameters*,double,Matrix*,Matrix*,Matrix*,int);

		void build_reduced_model(int,double);
		void test_reduced_model();

		void single_sat();
};

#endif
