#ifndef MATRIX_COMPLEX_H
#define MATRIX_COMPLEX_H

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <vector>
#include <cassert>
#include <fstream>
#include <cmath>
#include <Accelerate/Accelerate.h>
#include "../libraries/math/matrix.h"

typedef __CLPK_integer integer;

using namespace std;

class Matrix_Complex
{
	private:
		Matrix real;
		Matrix imag;
	public:
		Matrix_Complex();
		Matrix_Complex(int);
		Matrix_Complex(int,int);

		void initiate_array(double*,double*,int,int);
		void initiate_vector(vector<double>,vector<double>,int,int);
		void initiate_vector_vector(vector< vector<double> >,vector< vector<double> >);
		void initiate_matrix(Matrix,Matrix);

		void zeros(int,int);
		void zeros(int);
		void ones(int,int);
		void ones(int);
		void rand(int,int);
		void rand(int);
		void clear();

		int get_num_rows();
		int get_num_cols();
		void get_size(int*,int*);
		Matrix get_real();
		Matrix get_imag();
		int length();

		double& at(int,int,char);	// 'r' for real part, 'i' for imaginary part
		double& at(int,char);

		void operator=(Matrix_Complex);
		Matrix_Complex operator+(Matrix_Complex);
		void operator+=(Matrix_Complex);
		Matrix_Complex operator-(Matrix_Complex);
		void operator-=(Matrix_Complex);

		Matrix_Complex conj();
		Matrix_Complex tr();

		void svd(Matrix_Complex*,Matrix*,Matrix_Complex*);

		friend ostream& operator<<(ostream&,Matrix_Complex&);
		void save();
};

#endif
