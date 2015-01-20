#include "polynomial.h"

Polynomial::Polynomial()
{
	degree = 0;
	coef = NULL;
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
	coef = new double[N+1];
}

void Polynomial::initiate_with_reference(int N , double* c)
{
	degree = N;
	coef = c;	
}
	
void Polynomial::set_coef(int index,double value)
{
	if( coef == NULL )
		cout << "Error - set_coef ! No polynomial detected" << endl;
	else if( index > degree )
		cout << "Error - set_coef ! Index out of bound" << endl;
	else
	{
		coef[index] = value;
	}
}

void Polynomial::set_coef_zeros()
{
	for(int i = 0 ; i < degree + 1 ; i++)
		coef[i]=0;
}
		
double* Polynomial::get_polynomial()
{
	return coef;
}
		
double Polynomial::get_coef(int index)
{
	if( coef == NULL )
	{
		cout << "Error - get_coef ! No polynomial detected" << endl;
		return 0;
	}
	else if( index > degree + 1 )
	{
		cout << "Error - get_coef ! Index out of bound" << endl;
		return 0;
	}
	else
	{
		return coef[index];
	}
}

int Polynomial::get_degree()
{
	return degree;
}

Polynomial Polynomial::scalar(double c)
{
	Polynomial result;
	result.initiate( degree );
	for(int i = 0 ; i < degree + 1 ; i++)
		result.set_coef(i,c*coef[i]);
	return result;
}

Polynomial Polynomial::operator+(Polynomial pol)
{
	Polynomial result;
	if( degree > pol.get_degree() )
	{
		result.initiate( degree );
		for(int i = 0 ; i < pol.get_degree() + 1 ; i++ )
			result.set_coef(i,coef[i]+pol.get_coef(i));
		for(int i = pol.get_degree() + 1 ; i < degree + 1 ; i++)
			result.set_coef(i,coef[i]);
	}
	else
	{
		result.initiate( pol.get_degree() );
		for(int i = 0 ; i < degree + 1 ; i++ )
			result.set_coef(i,coef[i]+pol.get_coef(i));
		for(int i = degree + 1 ; i < pol.get_degree() + 1 ; i++)
			result.set_coef(i,pol.get_coef(i));
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
			result.set_coef(i,coef[i]-pol.get_coef(i));
		for(int i = pol.get_degree() + 1 ; i < degree + 1 ; i++)
			result.set_coef(i,coef[i]);
	}
	else
	{
		result.initiate( pol.get_degree() );
		for(int i = 0 ; i < degree + 1 ; i++ )
			result.set_coef(i,coef[i]-pol.get_coef(i));
		for(int i = degree + 1 ; i < pol.get_degree() + 1 ; i++)
			result.set_coef(i,-pol.get_coef(i));
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
			result.set_coef( i+j , result.get_coef(i+j) + coef[i] * pol.get_coef( j ) );
		}
	}

	return result;
}

void Polynomial::operator=(Polynomial pol)
{
	if(coef)
		delete[] coef;
	initiate(pol.degree);
	for(int i = 0 ; i < degree + 1 ; i++)
		coef[i] = pol.coef[i]; 
}


double Polynomial::evaluate(double x)
{
	double result = 0;
	if(!coef)
		cout << "Error - evaluate ! No polynomial detected" << endl;
	else
	{
		result = coef[degree];
		for(int i = degree - 1 ; i >= 0 ; i--)
			result = coef[i] + result*x;
	}
	return result;
}

double* Polynomial::extract_roots()
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
		*(temp) = -coef[ degree - i - 1 ] / coef[ degree ];
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

	delete[] A;
	delete[] WI;
	delete[] WORK;

	return WR;
}

void Polynomial::print()
{
	if( coef == NULL )
		cout << "Error - print ! No polynomial detected" << endl;
	else
	{
		for(int i = 0 ; i < degree + 1 ; i++)
			cout << coef[i] << " " ;
		cout << endl;
	}
}
