#include "basis.h"

Basis::Basis()
{
	max_degree = 0;
	pol = NULL;
}

void Basis::initialize_hermite(int N)
{
	max_degree = N;
	pol = new Polynomial[N+1];

	vector<double> temp;

	pol[0].initiate(0);
	pol[0].set_coef(0,1);

	pol[1].initiate(1);
	pol[1].set_coef(0,0);
	pol[1].set_coef(1,1);

	temp.push_back(1);
	temp.push_back(0);
	temp.push_back(0);
	recurrence.push_back(temp);

	for(int i = 2 ; i < max_degree + 1 ; i++)
	{
		pol[i].initiate(i);
		pol[i].set_coef(0, - (i-1) * pol[i-2].get_coef(0) );
		for(int j = 1 ; j < i - 1 ; j++)
		{
			pol[i].set_coef( j , pol[i-1].get_coef(j-1) - (i-1) * pol[i-2].get_coef(j) );
		}
		pol[i-1].set_coef( i-1,pol[i-2].get_coef(i-2) );
		pol[i].set_coef( i , pol[i-1].get_coef(i-1) );
		temp.clear();
		temp.push_back(1);
		temp.push_back(0);
		temp.push_back(i-1);
		recurrence.push_back(temp);
	}
	temp.clear();
}

Polynomial* Basis::get_element(int degree)
{
	return pol+degree;
}

void Basis::get_quadrature_points(Matrix* quad,Matrix* weights)
{
	double* D = new double[max_degree];
	for(int i = 0 ; i < max_degree ; i++)
		D[i] = -(recurrence[i])[1] / (recurrence[i][0]);

	double* E = new double[max_degree];
	E[0]=0;
	for(int i = 0 ; i < max_degree-1 ; i++)
		E[i] = sqrt( (recurrence[i+1])[2] / ( (recurrence[i])[0] * (recurrence[i+1])[0] ) );

	char JOBZ = 'V';
	__CLPK_integer N = max_degree;

	double* Z = new double[max_degree*max_degree];
	__CLPK_integer LDZ = N;
	double* WORK = new double[2*max_degree - 2];
	__CLPK_integer INFO;

	dstev_(&JOBZ,&N,D,E,Z,&LDZ,WORK,&INFO);

	quad->zeros(1,max_degree);
	weights->zeros(1,max_degree);

	for(int i = 0 ; i < max_degree ; i++)
	{
		quad->at(i) = D[i];
		weights->at(i) = Z[i*max_degree] * Z[i*max_degree];
	}

	delete[] D;
	delete[] E;
	delete[] Z;
}

vector< vector<int> > Basis::compute_mixed_expansions(int N, int d)
{
//	int basis_size = factorial(N+d-1) / ( factorial(N)*factorial(d-1) );

	vector< vector<int> > degrees;
	vector<int> instance;

	for(int i = 0 ; i < d ; i++)
		instance.push_back(0);
	
	degrees.push_back(instance);

	int ptr;
	int temp;
	for(int n = 1 ; n < N+1 ; n++)
	{
		instance.clear();
		for(int i = 0 ; i < d ; i++)
			instance.push_back(0);
		
		instance[0] = n;
		degrees.push_back(instance);
		ptr = 0;
		while( instance[d-1] != n )
		{
//			cout << "here 1" << endl;
			while(ptr != d-1)
			{
				instance[ptr] -= 1;
				ptr++;
				instance[ptr] += 1;
				degrees.push_back(instance);
			}

			temp = instance[ptr];
			instance[ptr] = 0;
			while( (instance[ptr] == 0) && (ptr > 0) )
				ptr--;

			if( instance[ptr] == 0 )
				instance[d-1] = n;
			else
			{
				instance[ptr] -= 1;
				ptr++;
				instance[ptr] = temp + 1;
				degrees.push_back(instance);
			}
		}
	}
	return degrees;
}

unsigned Basis::factorial(unsigned n)
{
	if (n == 1)
		return 1;
	else
		return n * factorial(n - 1);
}

void Basis::print_basis()
{
	for(int i = 0 ; i < max_degree + 1 ; i++)
		pol[i].print();
}
