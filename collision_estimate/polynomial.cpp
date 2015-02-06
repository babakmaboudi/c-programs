#include "polynomial.h"

Polynomial::Polynomial()
{
	degree = 0;
}

/*
Polynomial::~Polynomial()
{
	if(coef != NULL)
		delete[] coef;
	cout << "I am in the destructor with value " << degree << endl;
}
*/


void Polynomial::initiate(int N)
{
	degree = N;
	coef.zeros(1,N+1);
}

void Polynomial::initiate_with_reference(int N , double* c)
{
	degree = N;
	coef.initiate_array(c,1,N+1);
}
	
void Polynomial::set_coef(int index,double value)
{
	coef.at(index) = value;
}

void Polynomial::set_coef_zeros()
{
	coef.zeros(1,degree + 1);
}
		
Matrix Polynomial::get_polynomial()
{
	return coef;
}
		
double Polynomial::get_coef(int index)
{
	return coef.at(index);
}

int Polynomial::get_degree()
{
	return degree;
}

Polynomial Polynomial::scalar(double c)
{
	Polynomial result;
	result.initiate(degree);
	result.coef = coef.scalar(c);
	return result;
}

Polynomial Polynomial::operator+(Polynomial pol)
{
	Polynomial result;
	if( degree > pol.get_degree() )
	{
		result.initiate( degree );
		for(int i = 0 ; i < pol.get_degree() + 1 ; i++ )
			result.set_coef( i , coef.at(i)+pol.coef.at(i) );
		for(int i = pol.get_degree() + 1 ; i < degree + 1 ; i++)
			result.set_coef( i , coef.at(i) );
	}
	else
	{
		result.initiate( pol.get_degree() );
		for(int i = 0 ; i < degree + 1 ; i++ )
			result.set_coef( i , coef.at(i) + pol.coef.at(i) );
		for(int i = degree + 1 ; i < pol.get_degree() + 1 ; i++)
			result.set_coef( i , pol.coef.at(i) );
	}
	return result;
}

Polynomial Polynomial::operator-(Polynomial pol)
{
	Polynomial result;
	if( degree > pol.get_degree() )
	{
		result.initiate( degree );
		for(int i = 0 ; i < pol.get_degree() + 1 ; i++ )
			result.set_coef( i , coef.at(i)-pol.coef.at(i) );
		for(int i = pol.get_degree() + 1 ; i < degree + 1 ; i++)
			result.set_coef( i , coef.at(i) );
	}
	else
	{
		result.initiate( pol.get_degree() );
		for(int i = 0 ; i < degree + 1 ; i++ )
			result.set_coef(i,coef.at(i)-pol.coef.at(i));
		for(int i = degree + 1 ; i < pol.get_degree() + 1 ; i++)
			result.set_coef(i,-pol.coef.at(i));
	}
	return result;
}

Polynomial Polynomial::operator*(Polynomial pol)
{
	Polynomial result;
	result.initiate( degree + pol.get_degree() );
	result.set_coef_zeros();
	for(int i = 0 ; i < degree + 1 ; i++)
	{
		for( int j = 0 ; j < pol.get_degree() + 1 ; j++)
		{
			result.set_coef( i+j , result.coef.at(i+j) + coef.at(i) * pol.coef.at(j) );
		}
	}

	return result;
}

void Polynomial::operator=(Polynomial pol)
{
	coef.clear();
	initiate(pol.degree);
	for(int i = 0 ; i < degree + 1 ; i++)
		coef = pol.coef; 
}


double Polynomial::evaluate(double x)
{
	double result = 0;
	result = coef.at(degree);
	for(int i = degree - 1 ; i >= 0 ; i--)
		result = coef.at(i) + result*x;
	return result;
}

Matrix Polynomial::extract_roots()
{
	double* A = new double[degree*degree];
	double* temp = A+1;
	for(int i = 1 ; i < degree ; i++)
	{
		*(temp) = 1;
		temp += degree+1;
	}

	temp = A;
	for(int i = 0 ; i < degree ; i++ )
	{
		*(temp) = -coef.at( degree - i - 1 ) / coef.at( degree );
		temp += degree;
	}
/*
	for(int i = 0 ; i < degree ; i++)
	{
		for(int j = 0 ; j < degree ; j++ )
		{
			cout << A[ degree*j + i ] << " ";
		}
		cout << endl;
	}
*/
	char JOBVL = 'N';
	char JOBVR = 'N';
	__CLPK_integer N = degree;

	__CLPK_integer LDA = degree;
	double* WR = new double[degree];
	double* WI = new double[degree];
	double* VL;
	__CLPK_integer LDVL = degree;
	double* VR;
	__CLPK_integer LDVR = degree;
	__CLPK_integer LWORK = -1;
	__CLPK_integer INFO;
	double WORK_SIZE;

	dgeev_(&JOBVL,&JOBVR,&N,A,&LDA,WR,WI,VL,&LDVL,VR,&LDVR,&WORK_SIZE,&LWORK,&INFO);
	LWORK = static_cast<int>(WORK_SIZE);
	double* WORK = new double[LWORK];
	dgeev_(&JOBVL,&JOBVR,&N,A,&LDA,WR,WI,VL,&LDVL,VR,&LDVR,WORK,&LWORK,&INFO);

	assert( INFO == 0 );	

	Matrix result;
	result.initiate_array(WR,1,degree);

	delete[] A;
	delete[] WI;
	delete[] WORK;
	delete[] WR;

	return result;
}

void Polynomial::print()
{
	cout << coef << endl;
}

