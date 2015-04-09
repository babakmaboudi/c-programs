#ifndef WAVE_H
#define WAVE_H

#include <iostream>
#include <cmath>
#include "parameters.h"
#include "../libraries/math/matrix.h"

class Wave
{
	private:
		Parameters param;

		Matrix Dxx;
		Matrix Phi;
	public:
		Wave();

		Matrix initial_condition(int);
		void compute_significant_subspace(Matrix*,Matrix*,double);
		
	
		void deriv_symp(Matrix*,double,Matrix*,Parameters*,int);
		
		void integ2_symp(void (Wave::*func)(Matrix*,double,Matrix*,Parameters*,int),double*,Matrix,Matrix,Parameters*,double,Matrix*,Matrix*,Matrix*,int);

		void build_reduced_model(int,double);
		void test_reduced_model();
		void solver(int);
};

#endif
