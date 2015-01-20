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

void prog_bar(int current, int total)
{
	int bar_width = 30;
	total--;
	int progress = static_cast<double> (current) / total * bar_width;

	cout << '\xd';
	cout << '[';
	for(int i = 0 ; i < progress ; i++)
		cout << '#';
	for(int i = 0 ; i < (bar_width - progress) ; i++)
		cout << " ";
	cout << ']';
	
	int percent = static_cast<double> (current) / total * 100;

	cout << " - " << percent << '%';

	if(current == (total) )
		cout << endl;
}

