#include "matrix.h"

Matrix::Matrix()
{
	nrows = 0;
	ncols = 0;
}

Matrix::~Matrix()
{
	mat.clear();
}

void Matrix::initiate_array(double* M, int m, int n)
{
	mat.clear();
	nrows = m;
	ncols = n;
	for(int i = 0 ; i < nrows ; i++)
	{
		for(int j = 0 ; j < ncols ; j++)
			mat.push_back( M[i*ncols + j] );
	}
}


void Matrix::initiate_vector(vector<double> M, int m, int n)
{
	mat = M;
	nrows = m;
	ncols = n;
}

void Matrix::initiate_vector_vector(vector< vector<double> > M)
{
	nrows = M.size();
	ncols = M[0].size();
	for(int i = 0 ; i < nrows ; i++)
	{
		for(int j = 0 ; j < ncols ; j++)
			mat.push_back( (M[i])[j] );
	}
}

void Matrix::zeros(int m, int n)
{
	mat.clear();
	nrows = m;
	ncols = n;
	for(int i = 0 ; i < m*n ; i++)
		mat.push_back(0);
}

void Matrix::ones(int m, int n)
{
	mat.clear();
	nrows = m;
	ncols = n;
	for(int i = 0 ; i < m*n ; i++)
		mat.push_back(1);
}

void Matrix::eye(int m)
{
	mat.clear();
	nrows = m;
	ncols = m;
	for(int i = 0 ; i < m ; i++)
	{
		for(int j = 0 ; j < m ; j++)
		{
			if(i == j)
				mat.push_back(1);
			else
				mat.push_back(0);
		}
	}
}
void Matrix::clear()
{
	nrows = 0;
	ncols = 0;
	mat.clear();
}

int Matrix::get_num_rows()
{
	return nrows;
}

int Matrix::get_num_cols()
{
	return ncols;
}

void Matrix::get_size(int* m, int* n)
{
	(*m) = nrows;
	(*n) = ncols;
}

int Matrix::length()
{
	if(nrows > ncols)
		return nrows;
	else
		return ncols;
}

double Matrix::get_element(int m, int n)
{
	if( ( m<nrows ) && ( n<ncols ) )
		return mat[m*ncols + n];
	else
	{
		cout << "Index out of bound!" << endl;
		return 0;
	}
}

vector<double> Matrix::get_matrix()
{
	return mat;
}

vector<double>* Matrix::get_matrix_ptr()
{
	return &mat;
}

Matrix Matrix::get_row(int m)
{
	Matrix result;
	if( m<nrows )
	{
		result.nrows = 1;
		result.ncols = ncols;
		for(int i = 0 ; i < ncols ; i++)
			(result.mat).push_back( mat[ m*ncols + i ] );
		return result;
	}
	else
	{
		cout << "Index out of bound!" << endl;
		return result;
	}
}

Matrix Matrix::get_col(int n)
{
	Matrix result;
	if( n<ncols )
	{
		result.nrows = nrows;
		result.ncols = 1;
		for(int i = 0 ; i < nrows ; i++)
			(result.mat).push_back( mat[ i*ncols + n ] );
		return result;
	}
	else
	{
		cout << "Index out of bound!" << endl;
		return result;
	}	
}

void Matrix::set_row(int pos, Matrix row)
{
	if( (pos < nrows) && (row.ncols == ncols) )
	{
		for(int i = 0 ; i < ncols ; i++)
			this->at(pos,i) = row.at(0,i);
	}
	else
		cout << "Index out of bound!" << endl;
}

void Matrix::set_col(int pos, Matrix col)
{
	if( (pos < ncols) && (col.nrows == nrows) )
	{
		for(int i = 0 ; i < nrows ; i++)
			this->at(i,pos) = col.at(i,0);
	}
	else
		cout << "Index out of bound!" << endl;
}

double& Matrix::at( int m, int n )
{
	return mat[m*ncols + n];
}

double& Matrix::at(int pos)
{
	if( (nrows != 1) && (ncols != 1) )
	{
		cout << "Matrix must be a vector!" << endl;
		return mat[0];
	}
	else
		return mat[pos];
}

void Matrix::add_row(Matrix M)
{
	if( (M.nrows == 1) && (M.ncols == ncols) )
	{
		nrows++;
		for(int i = 0 ; i < ncols ; i++)
			mat.push_back( M.at(0,i) );
	}
	else if((nrows == 0) && (ncols == 0))
	{
		if( M.nrows == 1 )
		{
			nrows = 1;
			ncols = M.ncols;
			for(int i = 0 ; i < M.ncols ; i++)
				mat.push_back( M.at(0,i) );
		}
		else
			cout << "Dimension inconsistant!" << endl;
	}
	else
		cout << "Dimension inconsistant!" << endl;
}

void Matrix::add_col(Matrix M)
{
	if((nrows == 0) && (ncols == 0))
	{
		if( M.ncols == 1 )
		{
			ncols = 1;
			nrows = M.nrows;
			for(int i = 0 ; i < nrows ; i++)
				mat.push_back(M.at(i,0));
		}
		else
			cout << "Dimension inconsistant!" << endl;
	}
	else
	{
		(*this) = this->tr();
		this->add_row( M.tr() );
		(*this) = this->tr();
	}
}

void Matrix::operator=(Matrix M)
{
	M.get_size(&nrows,&ncols);
	mat = M.mat;
}

Matrix Matrix::operator+(Matrix M)
{
	Matrix result;
	if( (nrows == M.nrows) && (ncols == M.ncols) )
	{
		result.nrows = nrows;
		result.ncols = ncols;
		for(int i = 0 ; i < mat.size() ; i++)
			(result.mat).push_back( mat[i] + M.mat[i]);
	}
	else
	{
		cout << "Inconsistant matrix size!" << endl;
	}
	return result;
}

Matrix Matrix::operator-(Matrix M)
{
	Matrix result;
	if( (nrows == M.nrows) && (ncols == M.ncols) )
	{
		result.nrows = nrows;
		result.ncols = ncols;
		for(int i = 0 ; i < mat.size() ; i++)
			(result.mat).push_back( mat[i] - M.mat[i]);
	}
	else
	{
		cout << "Inconsistant matrix size!" << endl;
	}
	return result;
}

Matrix Matrix::operator*(Matrix M)
{
	Matrix result;
	if(ncols == M.nrows)
	{
		result.nrows = nrows;
		result.ncols = M.ncols;
		double temp;
		for(int i = 0 ; i < nrows ; i++)
		{
			for(int j = 0 ; j < M.ncols ; j++)
			{
				temp = 0;
				for(int k = 0 ; k < ncols ; k++)
				{
					temp += mat[ i*ncols+k ]*M.mat[k*M.ncols + j];
				}
				(result.mat).push_back(temp);
			}
		}
	}
	else
	{
		cout << "Inconsistant matrix size!" << endl;
	}
	return result;
}

Matrix Matrix::scalar(double c)
{
	Matrix result;
	if( (nrows) && (ncols) )
	{
		result.nrows = nrows;
		result.ncols = ncols;
		for(int i = 0 ; i < mat.size() ; i++)
			(result.mat).push_back(c * mat[i]);
	}
	else
	{
		cout << "No matrix found!" << endl;
	}
	return result;
}

Matrix Matrix::tr()
{
	Matrix result;
	result.nrows = ncols;
	result.ncols = nrows;
	for(int i = 0 ; i < ncols ; i++)
	{
		for(int j = 0 ; j < nrows ; j++)
			(result.mat).push_back( mat[ j*ncols + i ] );
	}
	return result;
}

double Matrix::norm2()
{
	double n = 0;
	for(int i = 0 ; i < mat.size() ; i++)
		n += mat[i]*mat[i];
	return sqrt(n);
}

void Matrix::svd(vector< vector<double> >* U_mat, vector<double>* S_vec, vector< vector<double> >* V_mat)
{
	double* A = new double[nrows*ncols];

	// matrix must be stored column wise
	for(int i = 0 ; i < ncols ; i++)
	{
		for(int j = 0 ; j < nrows ; j++)
			A[ i*nrows + j ] = mat[j*ncols + i];
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

void Matrix::svd(Matrix* U_mat, Matrix* S_vec, Matrix* V_mat)
{
	double* A = new double[nrows*ncols];

	// matrix must be stored column wise
	for(int i = 0 ; i < ncols ; i++)
	{
		for(int j = 0 ; j < nrows ; j++)
			A[ i*nrows + j ] = mat[j*ncols + i];
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

	U_mat->clear();
	U_mat->zeros(nrows,min);
	for(int i = 0 ; i < min ; i++)
	{
		for(int j = 0 ; j < nrows ; j++)
			 U_mat->at(j,i) = U[i*nrows + j];
	}

	S_vec->clear();
	S_vec->zeros(min,1);
	for(int i = 0 ; i < min ; i++)
		S_vec->at(i) = S[i];

	V_mat->clear();
	V_mat->zeros(min,ncols);
	for(int i = 0 ; i < min ; i++)
	{
		for(int j = 0 ; j < ncols ; j++)
			V_mat->at(i,j) = VT[j*min + i];
	}
	(*V_mat) = V_mat->tr();

	delete[] S;
	delete[] U;
	delete[] VT;
}


Matrix Matrix::inv()
{
	Matrix result;
	if( nrows == ncols )
	{
		double* A = new double[ nrows*nrows ];
		// matrix must be stored column wise
		for(int i = 0 ; i < ncols ; i++)
		{
			for(int j = 0 ; j < nrows ; j++)
				A[ i*nrows + j ] = mat[j*ncols + i];
		}

		__CLPK_integer N = nrows;
		__CLPK_integer LDA = nrows;
		__CLPK_integer* IPIV = new __CLPK_integer[ nrows+1 ];
		double work_size;
		__CLPK_integer LWORK = -1;
		__CLPK_integer INFO;

		dgetrf_(&N,&N,A,&LDA,IPIV,&INFO);
		assert(INFO == 0);

		dgetri_(&N,A,&LDA,IPIV,&work_size,&LWORK,&INFO);

		LWORK = static_cast<int>(work_size);
		double* WORK = new double[LWORK];

		dgetri_(&N,A,&LDA,IPIV,WORK,&LWORK,&INFO);

		assert(INFO == 0);

		delete[] WORK;

		result.zeros(nrows,nrows);
		for(int i = 0 ; i < nrows ; i++)
		{
			for(int j = 0 ; j < ncols ; j++)
				result.at(i,j) = A[ j*nrows + i ];
		}

		delete[] A;
	}
	else
	{
		cout << "Matrix must be a square!" << endl;
	}
	return result;
}

ostream& operator<<(ostream& os,Matrix& M)
{
	cout << "---------------------" << endl;
	for(int i = 0 ; i < M.nrows ; i++)
	{
		for(int j = 0 ; j < M.ncols ; j++)
		{
			cout.width(17);
			cout << std::setprecision(10) << M.mat[i*M.ncols + j];
		}
		cout << endl;
	}
	cout << "---------------------";
	return os;
}

void Matrix::print()
{
	if(mat.size() == 0)
		cout << "No matrix found!" << endl;
	else
	{
		for(int i = 0 ; i < nrows ; i++)
		{
			for(int j = 0 ; j < ncols ; j++)
				cout << mat[i*ncols + j] << " ";
			cout << endl;
		}
		cout << "---------------------------" << endl;
	}
}

void Matrix::save(char* path)
{
	ofstream file;
	file.open(path);

	for(int i = 0 ; i < nrows ; i++)
	{
		for(int j = 0 ; j < ncols ; j++)
			file << std::setprecision(16) << mat[i*ncols + j] << " ";
		file << endl;
	}
	file.close();	
}
