#include "tools.h"

//	'b' for maximum and 's' for minimum
void arg_opt(double* opt, int* idx, Matrix* M, char order)
{
	vector<double> mat = M->get_matrix();
	(*opt) = mat[0];
	(*idx) = 0;

	int iter=1;

	if(order == 'b')
	{
		while(iter < mat.size() )
		{
			if( abs( mat[iter] ) > abs( (*opt) ) )
			{
				(*opt) = mat[iter];
				(*idx) = iter;
			}
			iter++;
		}
	}
	else
	{
		while(iter < mat.size() )
		{
			if( abs( mat[iter] ) < abs( (*opt) ) )
			{
				(*opt) = mat[iter];
				(*idx) = iter;
			}
			iter++;
		}
	}
	

}

vector<double> linspace(double min, double max, int n)
{
	vector<double> result;
	double temp;

	for(int i = 0 ; i <= n-2 ; i++)
	{
		temp = min + i*(max-min)/(floor((double)n) - 1) ;
		result.push_back( temp );
	}
	result.push_back(max);
	return result;
}

