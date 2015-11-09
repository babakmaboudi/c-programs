#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <vector>
#include <cassert>
#include "vec.h"

using namespace std;

class Matrix
{
	private:
		int nrows;
		int ncols;
		vector<real> mat;
	public:
		Matrix();
		Matrix(int,int);

		void set(int, int);
		
		real& at(int,int);
		real& at(int);

		Matrix tr();
		real& operator()(int,int);
		void operator=(Matrix);
		Vec operator*(Vec);

		void operator()();
		void print();
};

#endif
