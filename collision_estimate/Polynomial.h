#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <iostream>
#include <cmath>
#include <cassert>
#include <Accelerate/Accelerate.h>
#include "../libraries/math/matrix.h"
typedef __CLPK_integer integer;

using namespace std;

class Polynomial
{
	private:
		int degree;
		Matrix coef;
	public:
		Polynomial();
//		~Polynomial();
		void initiate(int);
		void initiate_with_reference(int,double*);
	
		void set_coef(int,double);
		void set_coef_zeros();

		Matrix get_polynomial();
		double get_coef(int);
		int get_degree();

		Polynomial scalar(double);
		Polynomial operator+(Polynomial);
		Polynomial operator-(Polynomial);
		Polynomial operator*(Polynomial);
		void operator=(Polynomial);

		double evaluate(double);
		Matrix extract_roots();

		void print();
};

#endif
