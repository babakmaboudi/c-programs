#ifndef MODEL_H
#define MODEL_H

#include <iostream>
#include "parameters.h"
#include "tools.h"
#include "../libraries/math/matrix.h"

using namespace std;

class Model
{
	private:
		Parameters param;
	public:
		void sv2oe(double*,double*,double*);
		void initiate_from_file(char*);

		void deriv(Matrix*,double,Matrix*,Parameters*);
		void euler(void (Model::*func)(Matrix*,double,Matrix*,Parameters*),double*,Matrix,Parameters*,double,Matrix*,Matrix*,int);
		
		void deriv_symp(Matrix*,double,Matrix*,Parameters*,int);
		void euler_symp(void (Model::*func)(Matrix*,double,Matrix*,Parameters*,int),double*,Matrix,Parameters*,double,Matrix*,Matrix*,int);

		void single_sat();
};

#endif
