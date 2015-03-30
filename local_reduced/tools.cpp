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

vector<double> get_polar_coord(Matrix* X)
{
	vector<double> r;
	r.push_back( sqrt( pow(X->at(0,0),2) + pow(X->at(1,0),2) + pow(X->at(2,0),2) ) );
	r.push_back( asin( X->at(2,0) / r[0] ) );
	double arcsin = asin( X->at(0,0) / (r[0] * cos( r[1] )) );
	double arccos = acos( X->at(1,0) / (r[0] * cos( r[1] )) );
	if(arcsin > 0)
	{
		if(arccos < M_PI/2)
			r.push_back( arcsin );
		else
			r.push_back( arccos );
	}
	else
	{
		if(arccos < M_PI/2)
			r.push_back( arcsin );
		else
			r.push_back( M_PI - arcsin );
	}
	return r;
}

Matrix compute_mixed_expansion(int N, int d)
{
	Matrix degrees;
	Matrix instance;
	instance.zeros(1,d);

	degrees.add_row(instance);


	int ptr;
	int temp;
	for(int n = 1 ; n < N+1 ; n++)
	{
		instance.zeros(1,d);
		
		instance.at(0) = n;
		degrees.add_row(instance);

		ptr = 0;
		while( instance.at(d-1) != n )
		{
			while(ptr != d-1)
			{
				instance.at(ptr) -= 1;
				ptr++;
				instance.at(ptr) += 1;
				degrees.add_row(instance);
			}

			temp = instance.at(ptr);
			instance.at(ptr) = 0;
			while( (instance.at(ptr) == 0) && (ptr > 0) )
				ptr--;

			if( instance.at(ptr) == 0 )
				instance.at(d-1) = n;
			else
			{
				instance.at(ptr) -= 1;
				ptr++;
				instance.at(ptr) = temp + 1;
				degrees.add_row(instance);
			}
		}
	}
	return degrees;
}

Matrix kmeans( Matrix* d , int nclusters , int niter )
{
	int ndata = d->get_num_rows();
	int dimension = d->get_num_cols();

	int* perm = randperm( ndata );
	Matrix clusters;
	clusters.zeros( nclusters , dimension );

	for( int i=0 ; i<nclusters ; i++ )
		clusters.set_row( i, d->get_row(perm[i]) );

	delete[] perm;

	double tol = 1e-3;

	Matrix assign, new_clusters, count, smat;
	double n;
	for( int iter=0 ; iter<niter ; iter++ )
	{
		assign = label_vectors( d , &clusters );
		new_clusters.clear();
		new_clusters.zeros( nclusters , dimension );
		count.clear();
		count.zeros( nclusters , 1 );

		for( int ii=0 ; ii<ndata ; ii++ )
		{
			int ind = static_cast<int>( assign.at(ii,0) );
			for( int jj=0 ; jj<dimension ; jj++ )
				new_clusters.at( ind , jj ) += d->at( ii , jj );
			count.at( ind , 0 ) += 1;
		}

		for( int ii=0 ; ii<new_clusters.get_num_rows() ; ii++ )
		{
			if( count.at( ii , 0 ) > 0 )
			{
				for( int jj=0 ; jj<new_clusters.get_num_cols() ; jj++ )
					new_clusters.at(ii,jj) /= static_cast<double>( count.at(ii,0) );

			}
		}

		smat = clusters - new_clusters;
		n = smat.norm2();
		clusters = new_clusters;

		if( n < tol )
			break;
	}

	return clusters;
}

Matrix label_vectors(Matrix* data, Matrix* centers)
{
	int ndata = data->get_num_rows();
      	int dim = data->get_num_cols();
	int nc = centers->get_num_rows();
	assert( centers->get_num_cols() == dim );	

	vector<double>* data_ptr = data->get_matrix_ptr();
	vector<double>* centers_ptr = centers->get_matrix_ptr();

	Matrix assign;
	assign.zeros( ndata , 1 );

	for( int i=0 ; i<ndata ; i++ )
	{
		double min_value = 1e50;
		double min_ind = -1;

		int index1 = i*dim;
		for( int j=0 ; j<nc ; j++ )
		{
			int index2 = j*dim;
			double v = 0.0;
			for( int k=0 ; k<dim ; k++ )
			{
				v += pow( (*data_ptr)[index1 + k] - (*centers_ptr)[index2 + k] , 2 );
				if( v > min_value )
					break;
			}
			if( v < min_value )
			{
				min_value = v;
				min_ind = j;
			}
		}
		assert( min_ind != -1 );
		assign.at(i,0) = min_ind;
	}

	return assign;
}

int* randperm(int n)
{
	int* perm = new int[n];
	for(int i = 0 ; i < n ; i++)
	{
		perm[i] = i;
	}

	int temp,r;
	for(int i = 0 ; i < n ; i++)
	{
		r = rand()%(n-i) + i;
		temp = perm[r];
		perm[r] = perm[i];
		perm[i] = temp;
	}
	return perm;
}

void find(vector<int>* res, Matrix* IDX, int instance)
{
	int len = IDX->length();
	
	for(int i = 0 ; i < len ; i++)
	{
		if( IDX->at(i) == instance )
			res->push_back(i);
	}
}
