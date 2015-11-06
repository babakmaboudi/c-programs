#ifndef MODEL_H
#define MODEL_H

#include <iostream>
#include <cmath>
#include <ctime>
#include "../libraries/math/matrix.h"

class Model
{
	private:
		int N;
		double L;
		double T;
		double dt;
		double dx;
		double c;
		int K;

		Matrix q;
		Matrix p;

		Matrix Dxx;

		Matrix J2n;
		Matrix J2k;
		Matrix A;
		Matrix A_inv;
		Matrix Ar;
		Matrix Aq;
		Matrix Ap;

		double TOL;
	public:
		Model();
		void initial_condition(Matrix*,Matrix*);
		void build_stiffness();

		void leap_frog(void (Model::*func)(Matrix*,int),Matrix*,Matrix*,Matrix,Matrix,int);
		void deriv(Matrix*,int);
		void deriv_reduced(Matrix*,int);

		void build_reduced_basis();

		void symplectic_qr(Matrix*,Matrix*);
		void symplectify(Matrix*,Matrix*,Matrix*);
		void normalize(Matrix*,Matrix*);
		double error_function(Matrix*,Matrix*, int ind);
		double hamiltonian(Matrix*,Matrix*, int ind);

		void simulate();
		void simulate_reduced();
};

#endif
