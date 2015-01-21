#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <vector>
#include <cassert>
#include <fstream>
#include <Accelerate/Accelerate.h>

typedef __CLPK_integer integer;


using namespace std;

class Matrix
{
	private:
		vector<double> mat;
		int nrows;
		int ncols;
	public:
		Matrix();
		~Matrix();

		void initiate_array(double*,int,int);
		void initiate_vector(vector<double>,int,int);
		void initiate_vector_vector(vector< vector<double> >);

		void zeros(int,int);
		void ones(int,int);
		void eye(int);
		void clear();

		int get_num_rows();
		int get_num_cols();
		void get_size(int*,int*);
		double get_element(int,int);
		vector<double> get_matrix();
		Matrix get_row(int);
		Matrix get_col(int);

		double& at(int,int);
		void add_row(Matrix);
		void add_col(Matrix);

		void operator=(Matrix);
		Matrix operator+(Matrix);
		Matrix operator-(Matrix);
		Matrix operator*(Matrix);
		Matrix scalar(double);
		Matrix tr();

		void svd(vector< vector<double> >*,vector<double>*,vector< vector<double> >*);
		Matrix inv();

		friend ostream& operator<<(ostream&,Matrix&);
		void print();
		void save(char*);
};

#endif
