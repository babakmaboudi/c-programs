#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <iostream>
#include <cmath>
#include <cassert>
#include <Accelerate/Accelerate.h>
#include <vector>
typedef __CLPK_integer integer;

using namespace std;

class Polynomial
{
	private:
		int degree;
		double* coef;
	public:
		Polynomial();
//		~Polynomial();
		void initiate(int);
		void initiate_with_reference(int,double*);
	
		void set_coef(int,double);
		void set_coef_zeros();
		
		double* get_polynomial();
		double get_coef(int);
		int get_degree();

		Polynomial scalar(double);
		Polynomial operator+(Polynomial);
		Polynomial operator-(Polynomial);
		Polynomial operator*(Polynomial);
		void operator=(Polynomial);

		double evaluate(double);
		double* extract_roots();

		void print();
};

#endif
