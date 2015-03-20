#include "matrix_complex.h"

Matrix_Complex::Matrix_Complex()
{

}

Matrix_Complex::Matrix_Complex(int length)
{
	real.zeros(length,1);
	imag.zeros(length,1);
}

Matrix_Complex::Matrix_Complex(int m, int n)
{
	real.zeros(m,n);
	imag.zeros(m,n);
}

void Matrix_Complex::initiate_array(double* r, double* i, int m, int n)
{
	real.initiate_array(r,m,n);
	imag.initiate_array(i,m,n);
}

void Matrix_Complex::initiate_vector(vector<double> r, vector<double> i, int m, int n)
{
	real.initiate_vector(r,m,n);
	imag.initiate_vector(i,m,n);
}

void Matrix_Complex::initiate_vector_vector(vector< vector<double> > r, vector< vector<double> > i)
{
	real.initiate_vector_vector(r);
	imag.initiate_vector_vector(i);
}

void Matrix_Complex::initiate_matrix(Matrix r, Matrix i)
{
	int r_nrows, r_ncols, i_nrows, i_ncols;
	r.get_size(&r_nrows,&r_ncols);
	i.get_size(&i_nrows,&i_ncols);
	assert( (r_nrows == i_nrows) && (r_ncols == i_ncols) );
	real = r;
	imag = i;
}

void Matrix_Complex::zeros(int m, int n)
{
	real.zeros(m,n);
	imag.zeros(m,n);
}

void Matrix_Complex::zeros(int length)
{
	real.zeros(length,1);
	imag.zeros(length,1);
}

void Matrix_Complex::ones(int m, int n)
{
	real.ones(m,n);
	imag.ones(m,n);
}

void Matrix_Complex::ones(int length)
{
	real.ones(length,1);
	imag.ones(length,1);
}

void Matrix_Complex::rand(int m, int n)
{
	real.rand(m,n);
	imag.rand(m,n);
}

void Matrix_Complex::rand(int m)
{
	real.rand(m);
	imag.rand(m);
}

void Matrix_Complex::clear()
{
	real.clear();
	imag.clear();
}

int Matrix_Complex::get_num_rows()
{
	return real.get_num_rows();
}

int Matrix_Complex::get_num_cols()
{
	return real.get_num_cols();
}

void Matrix_Complex::get_size(int* m, int* n)
{
	(*m) = real.get_num_rows();
	(*n) = real.get_num_cols();
}

Matrix Matrix_Complex::get_real()
{
	return real;
}

Matrix Matrix_Complex::get_imag()
{
	return imag;
}

int Matrix_Complex::length()
{
	return real.length();
}

double& Matrix_Complex::at(int m, int n, char ind)
{
	if(ind == 'r')
		return real.at(m,n);
	else
		return imag.at(m,n);
}

double& Matrix_Complex::at(int m, char ind)
{
	if(ind == 'r')
		return real.at(m);
	else
		return imag.at(m);
}

void Matrix_Complex::operator=(Matrix_Complex mat)
{
	real = mat.real;
	imag = mat.imag;
}

Matrix_Complex Matrix_Complex::operator+(Matrix_Complex mat)
{
	Matrix_Complex result;
	result.real = real + mat.real;
	result.imag = imag + mat.imag;
	return result;
}

void Matrix_Complex::operator+=(Matrix_Complex mat)
{
	real += mat.real;
	imag += mat.imag;
}

Matrix_Complex Matrix_Complex::operator-(Matrix_Complex mat)
{
	Matrix_Complex result;
	result.real = real - mat.real;
	result.imag = imag - mat.imag;
	return result;
}

void Matrix_Complex::operator-=(Matrix_Complex mat)
{
	real -= mat.real;
	imag -= mat.imag;
}

Matrix_Complex Matrix_Complex::conj()
{
	Matrix_Complex result;
	result.real = real;
	result.imag = imag * (-1.0);
	return result;
}

Matrix_Complex Matrix_Complex::tr()
{
	Matrix_Complex result;
	result.real = real.tr();
	result.imag = imag.tr() * (-1.0);
	return result;
}

void Matrix_Complex::svd(Matrix_Complex* U_mat, Matrix* S_vec, Matrix_Complex* V_mat)
{
	int nrows,ncols;
	int flag = 0;
	if(real.get_num_rows() < real.get_num_cols())
	{
		(*this) = this->tr();
		flag = 1;
	}

	real.get_size(&nrows,&ncols);
	__CLPK_doublecomplex* A = new __CLPK_doublecomplex[ nrows*ncols ];

	for(int i = 0 ; i < ncols ; i++)
	{
		for(int j = 0 ; j < nrows ; j++)
			A[ i*nrows + j ].r = real.at(j,i);
	}
	for(int i = 0 ; i < ncols ; i++)
	{
		for(int j = 0 ; j < nrows ; j++)
			A[ i*nrows + j ].i = imag.at(j,i);
	}

	int min;
	if( nrows < ncols )
		min = nrows;
	else
		min = ncols;

	__CLPK_integer M = nrows;
	__CLPK_integer N = ncols;
	char JOBZ='S';
	__CLPK_integer LDA = real.length();
	__CLPK_doublereal* S = new __CLPK_doublereal[min];
	__CLPK_doublecomplex* U = new __CLPK_doublecomplex[ M*min ];
	__CLPK_integer LDU = M;
	__CLPK_doublecomplex* VT = new __CLPK_doublecomplex[ min*N ];
	__CLPK_integer LDVT = N;

	__CLPK_doublecomplex worksize;
	__CLPK_integer LWORK = -1;
	double rwork_size;
	__CLPK_integer* IWORK = new __CLPK_integer[ 8*min ];
	__CLPK_integer INFO;

	zgesdd_(&JOBZ,&M,&N,A,&LDA,S,U,&LDU,VT,&LDVT,&worksize,&LWORK,&rwork_size,IWORK,&INFO);
	LWORK = static_cast<integer>( worksize.r );
	__CLPK_doublecomplex* WORK = new __CLPK_doublecomplex[ LWORK ];
	double* RWORK = new double[ LWORK ];
	zgesdd_(&JOBZ,&M,&N,A,&LDA,S,U,&LDU,VT,&LDVT,WORK,&LWORK,RWORK,IWORK,&INFO);
	assert(INFO == 0);

	U_mat->zeros(M,min);
	for(int i = 0 ; i < M ; i++)
	{
		for(int j = 0 ; j < min ; j++)
		{
			U_mat->at(i,j,'r') = U[ j*M + i ].r;
			U_mat->at(i,j,'i') = U[ j*M + i ].i;
		}
	}

	S_vec->zeros(min,1);
	for(int i = 0 ; i < min ; i++)
		S_vec->at(i) = S[i];

	V_mat->zeros(min,N);
	for(int i = 0 ; i < min ; i++)
	{
		for(int j = 0 ; j < N ; j++)
		{
			V_mat->at(i,j,'r') = VT[ i*N + j ].r;
			V_mat->at(i,j,'i') = -1.0*VT[ i*N + j ].i;
		}
	}

	if(flag == 1)
	{
		(*this) = this->tr();
		Matrix_Complex temp = (*U_mat);
		(*U_mat) = (*V_mat);
		(*V_mat) = temp;
	}

	delete[] RWORK;
	delete[] WORK;
	delete[] IWORK;
	delete[] VT;
	delete[] U;
	delete[] S;
	delete[] A;
}

ostream& operator<<(ostream& os,Matrix_Complex& M)
{
	cout << "real part:" << endl;
	cout << M.real << endl;
	cout << "imaginary part:" << endl;
	cout << M.imag << endl;
	return os;
}

void Matrix_Complex::save()
{
	char path_real[] = "out_real.txt"; 
	char path_imag[] = "out_imag.txt";
	real.save(path_real);
	imag.save(path_imag);
}
