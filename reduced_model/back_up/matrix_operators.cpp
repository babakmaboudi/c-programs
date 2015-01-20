#include "matrix_operators.h"

/* 
This routine comoputes the svd decomposition
of an m x n matrix A. such that A = U*S*V'.
The first argument is the input matrix to be
decomposed. The outputmatrices U_mat , S_vec
and V_mat are passed by reference as a 
vector<vector<double>>. each vector of U_mat
contain one right singular vector and each 
vector of V_mat contain one right singular 
vector. this means: 
	A = U_mat' * diag(S) * V_mat.
*/ 

void svd(vector< vector<double> > mat,vector< vector<double> >* U_mat, vector<double>* S_vec, vector< vector<double> >* V_mat)
{
	int nrows = mat.size();
	int ncols = mat[0].size();

	double* A = new double[nrows*ncols];

	// matrix must be stored column wise
	for(int i = 0 ; i < ncols ; i++)
	{
		for(int j = 0 ; j < nrows ; j++)
			A[ i*nrows + j ] = (mat[j])[i];
	}

	char JOBU = 'S';
	char JOBVT = 'S';
	__CLPK_integer M = nrows;
	__CLPK_integer N = ncols;
	__CLPK_integer LDA = nrows;
	int min;
	if(nrows > ncols)
		min = ncols;
	else
		min = nrows;
	double* S = new double[min];
	double* U = new double[nrows * min];
	__CLPK_integer LDU = nrows;

	double* VT = new double[min * ncols];
	__CLPK_integer LDVT = min;

	double work_size;
	__CLPK_integer LWORK = -1;
	__CLPK_integer INFO;

	dgesvd_(&JOBU,&JOBVT,&M,&N,A,&LDA,S,U,&LDU,VT,&LDVT,&work_size,&LWORK,&INFO);
	LWORK = static_cast<int>(work_size);
	double* WORK = new double[LWORK];
	
	dgesvd_(&JOBU,&JOBVT,&M,&N,A,&LDA,S,U,&LDU,VT,&LDVT,WORK,&LWORK,&INFO);

	assert(INFO == 0);

	delete[] A;
	delete[] WORK;

	vector<double> temp;
	for(int i = 0 ; i < min ; i++)
	{
		temp.clear();
		for(int j = 0 ; j < nrows ; j++)
			temp.push_back( U[i*min + j] );
		U_mat->push_back(temp);
	}

	for(int i = 0 ; i < min ; i++)
		S_vec->push_back( S[i] );

/*
	for(int i = 0 ; i < ncols ; i++)
	{
		temp.clear();
		for(int j = 0 ; j < min ; j++)
			temp.push_back( VT[i*min + j] );
		V_mat->push_back(temp);
	}
*/
	for(int i = 0 ; i < min ; i++)
	{
		temp.clear();
		for(int j = 0 ; j < ncols ; j++)
			temp.push_back( VT[j*min + i] );
		V_mat->push_back( temp );
	}

	delete[] S;
	delete[] U;
	delete[] VT;
}
