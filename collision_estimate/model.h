#ifndef MODEL_H
#define MODEL_H

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <cmath>
#include <vector>
#include "parameters.h"
#include "tools.h"
#include "basis.h"
#include "../libraries/math/matrix.h"

class Model
{
	private:
		Parameters param;
	public:
		void sv2oe(double*,double*,double*);
		void initiate_from_file(char*);

		void dynamic_3d_param(Matrix*,double,Matrix*,Parameters*);
		void explicit_rk6(void (Model::*func)(Matrix*,double,Matrix*,Parameters*),double*,Matrix,Parameters*,double,Matrix*,Matrix*,int);

		void single_sat();

		void SC_gaussian(int,int);
};

#endif
