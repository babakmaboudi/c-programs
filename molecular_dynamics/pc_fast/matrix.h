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
		
		real& at(int,int);

		Matrix tr();
		void operator=(Matrix);
		Vec operator*(Vec);

		void print();
};

#endif
