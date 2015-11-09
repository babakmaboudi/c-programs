#ifndef MODEL_H
#define MODEL_H

#include <iostream>
#include <cmath>
#include "tools.h"
#include "../libraries/math/matrix.h"

using namespace std;

class Model
{
	private:
		double L;
		int N;
		double dx;
		double x0;

		double T;
		double dt;
		int MAX_ITER;
		double v;

		Matrix Dxx;
		Matrix bd;
	public:
		void initiate();

		void initial_condition(Matrix*, Matrix*);
		void sine_Gordon();

		void deriv_symp(Matrix*,double,Matrix*,int);
		void nonlin(Matrix*,Matrix*);

		void integ2_symp(void (Model::*func)(Matrix*,double,Matrix*,int),double*,Matrix,Matrix,double,Matrix*,Matrix*,Matrix*);
		void integ4_symp(void (Model::*func)(Matrix*,double,Matrix*,int),double*,Matrix,Matrix,double,Matrix*,Matrix*,Matrix*);

		void DEI();
};

#endif
