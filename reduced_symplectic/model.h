#ifndef MODEL_H
#define MODEL_H

#include <iostream>
#include <cassert>
#include "parameters.h"
#include "tools.h"
#include "grid.h"
#include "../libraries/math/matrix.h"

using namespace std;

class Model
{
	private:
		Parameters param;

		Matrix phi1;
		Matrix phi2;
		Matrix L;
	public:
		void sv2oe(double*,double*,double*);
		void initiate_from_file(char*);

		void deriv(Matrix*,double,Matrix*,Parameters*);
		void euler(void (Model::*func)(Matrix*,double,Matrix*,Parameters*),double*,Matrix,Parameters*,double,Matrix*,Matrix*,int);
		
		void deriv_symp(Matrix*,double,Matrix*,Parameters*,int);
		void deriv_symp_test(Matrix*,double,Matrix*,Parameters*,int);
		void deriv_symp_reduced(Matrix*,double,Matrix*,Parameters*,int);
		void euler_symp(void (Model::*func)(Matrix*,double,Matrix*,Parameters*,int),double*,Matrix,Parameters*,double,Matrix*,Matrix*,int);
		void integ2_symp(void (Model::*func)(Matrix*,double,Matrix*,Parameters*,int),double*,Matrix,Parameters*,double,Matrix*,Matrix*,int);
		void integ4_symp(void (Model::*func)(Matrix*,double,Matrix*,Parameters*,int),double*,Matrix,Parameters*,double,Matrix*,Matrix*,int);
		void integ4_symp_reduced(void (Model::*func)(Matrix*,double,Matrix*,Parameters*,int),double*,Matrix,Matrix,Parameters*,double,Matrix*,Matrix*,Matrix*,int);

		void build_reduced_model(int,double);
		void test_reduced_model();

		void single_sat();
};

#endif
