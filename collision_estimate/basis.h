#ifndef BASIS_H
#define BASIS_H

#include <iostream>
#include <cmath>
#include "polynomial.h"
#include <vector>
#include <Accelerate/Accelerate.h>
#include "../libraries/math/matrix.h"
typedef __CLPK_integer integer;



using namespace std;

class Basis
{
	private:
		int max_degree;
		Polynomial* pol;
		vector< vector<double> > recurrence;
	public:
		Basis();

		void initialize_hermite(int);

		Polynomial* get_element(int);

		void get_quadrature_points(Matrix*,Matrix*);

		static vector< vector<int> > compute_mixed_expansions(int,int);
		static unsigned factorial(unsigned);

		void print_basis();
};

#endif
