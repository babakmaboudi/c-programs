#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <vector>
#include <cassert>
#include <fstream>
#include <cmath>
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
		Matrix(int);
		Matrix(int,char);	// 'c' for column vector 'r' for row vector
		Matrix(int,int);
		~Matrix();

		void initiate_array(double*,int,int);
		void initiate_vector(vector<double>,int,int);
		void initiate_vector_vector(vector< vector<double> >);

		void zeros(int,int);
		void ones(int,int);
		void eye(int);
		void jay(int);
		void rand(int,int);
		void rand(int);
		void clear();

		int get_num_rows();
		int get_num_cols();
		void get_size(int*,int*);
		int length();
		double get_element(int,int);
		vector<double> get_matrix();
		vector<double>* get_matrix_ptr();
		Matrix get_row(int);
		Matrix get_col(int);
		Matrix get_submat(int,int,int,int);

		void set_row(int,Matrix);
		void set_col(int,Matrix);
		void set_submat(int,int,Matrix);

		double& at(int,int);
		double& at(int);
		void add_row(Matrix);
		void add_col(Matrix);
		void append(Matrix,char);	// 'r' for row-wise append and 'c' for column-wise append

		void operator=(Matrix);
		Matrix operator+(Matrix);
		void operator+=(Matrix);
		Matrix operator-(Matrix);
		void operator-=(Matrix);
		Matrix operator*(Matrix);
		Matrix operator+(double);
		void operator+=(double);
		Matrix operator-(double);
		void operator-=(double);
		Matrix operator*(double);
		void operator*=(double);
		Matrix operator/(double);
		void operator/=(double);

		Matrix scalar(double);
		Matrix tr();

		double trace();
		double norm2();
		double frobenius();
		void linspace(double,double,int);

		void svd(vector< vector<double> >*,vector<double>*,vector< vector<double> >*);
		void svd(Matrix*,Matrix*,Matrix*);
		Matrix inv();
		void lu(Matrix*,int*);
		double det();
		void eig(Matrix*);
		void eig(Matrix*,Matrix*,Matrix*);
		void eig_sym(Matrix*);
		void eig_sym(Matrix*,Matrix*);

		friend ostream& operator<<(ostream&,Matrix&);
		void print();
		void save(char*);
		void open(int,int,char*);
};

#endif
