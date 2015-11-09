#include "matrix.h"

Matrix::Matrix()
{
	nrows = 0;
	ncols = 0;
}

Matrix::Matrix(int num_rows, int num_cols)
{
	nrows = num_rows;
	ncols = num_cols;
	mat = vector<real>(nrows*ncols,0.0);
}

void Matrix::set(int num_rows, int num_cols)
{
	nrows = num_rows;
	ncols = num_cols;
	mat = vector<real>(nrows*ncols,0.0);
}

real& Matrix::at(int i, int j)
{
	assert(i<nrows && j<ncols);
	return mat[i*ncols + j];
}

real& Matrix::at(int i)
{
	assert(i<mat.size());
	return mat[i];
}

real& Matrix::operator()(int i, int j)
{
	assert(i<nrows && j<ncols);
	return mat[i*ncols + j];
}

Matrix Matrix::tr()
{
	Matrix res;
	for(int i=0 ; i<nrows ; i++)
		for(int j=0 ; j<ncols ; j++)
			res.at(i,j) = at(j,i);
	return res;
}

void Matrix::operator=(Matrix mat)
{
	assert(nrows == mat.nrows && ncols == mat.ncols);
	for(int i=0 ; i<nrows ; i++)
		for(int j=0 ; j<ncols ; j++)
			at(i,j) = mat.at(i,j);
}

Vec Matrix::operator*(Vec v)
{
	Vec res;
	res.x() = v.x() * at(0,0) + v.y() * at(0,1) + v.z() * at(0,2);
	res.y() = v.x() * at(1,0) + v.y() * at(1,1) + v.z() * at(1,2);
	res.z() = v.x() * at(2,0) + v.y() * at(2,1) + v.z() * at(2,2);
	return res;
}

void Matrix::operator()()
{
	for(int i=0 ; i<nrows ; i++)
	{
		for(int j=0 ; j<ncols ; j++)
		{
			cout.width(10);
			cout << right << (*this)(i,j);
		}
		cout << endl;
	}
	cout << "-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-" << endl;	
}

void Matrix::print()
{
	for(int i=0 ; i<nrows ; i++)
	{
		for(int j=0 ; j<ncols ; j++)
		{
			cout.width(10);
			cout << right << mat[i*nrows + j];
		}
		cout << endl;
	}
	cout << "-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-" << endl;
}
