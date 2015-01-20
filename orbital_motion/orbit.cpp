#include "orbit.h"

Orbit::Orbit()
{
	a = 0;
	e = 0;
	i = 0;
	omega = 0;
	w = 0;
	nu = 0;
	X0 = new double[3];
	V0 = new double[3];
}

void Orbit::set_orbit_constants(double G_value, double M_value )
{
	G = G_value;
	M = M_value;
	mu = G*M;
}

void Orbit::set_initial_position( double value_a, double value_e, double value_i, double value_omega, double value_w, double value_nu )
{
	a = value_a;
	e = value_e;
	i = value_i;
	omega = value_omega;
	w = value_w;
	nu = value_nu;
}

void Orbit::oe2sv()
{
	double p = a*(1-e);
	X0[0] = p*( cos(omega) * cos( w + nu ) - sin(omega) * cos( i ) * sin( w + nu ) );
	X0[1] = p*( sin(omega) * cos( w + nu ) + cos(omega) * cos( i ) * sin( w + nu ) );
	X0[2] = p * sin( i ) * sin( w + nu );

	p = a*(1-pow(e,2));	
	double term = - sqrt( mu/p );
//	cout << term << endl;
	

	V0[0] = term * ( cos(omega)*( sin(w + nu) + e*sin(w) ) + sin(omega)*cos(i)*(cos(w+nu) + e*cos(w)) );
	V0[1] = term * ( sin(omega)*( sin(w + nu) + e*sin(w) ) - cos(omega)*cos(i)*(cos(w+nu) + e*cos(w)) );
	V0[2] = term * sin(i) * ( cos(w + nu) + e * cos(w) );
//	for(int i = 0 ; i < 3 ; i++)
//		cout << X0[i] << " ";
//	for(int i = 0 ; i < 3 ; i++)
//		cout << V0[i] << " ";
//	cout << endl;
}

double* Orbit::get_state_vec()
{
	return X0;
}

void Orbit::dynamic_2d(Vec* dy, double t, Vec* y)
{
	double coef = - mu / pow( pow(y->get_element(0),2) + pow(y->get_element(1),2) , 1.5 ) ;
	dy->clear();
	dy->add_element( y->get_element(2) );
	dy->add_element( y->get_element(3) );
	dy->add_element( coef * y->get_element(0) );
	dy->add_element( coef * y->get_element(1) );
}

void Orbit::dynamic_3d(Vec* dy, double t, Vec* y)
{
	double coef = - mu / pow( pow(y->get_element(0),2) + pow(y->get_element(1),2) + pow(y->get_element(2),2) , 1.5 );
	dy->clear();
	dy->add_element( y->get_element(3) );
	dy->add_element( y->get_element(4) );
	dy->add_element( y->get_element(5) );
	dy->add_element( coef * y->get_element(0) );
	dy->add_element( coef * y->get_element(1) );
	dy->add_element( coef * y->get_element(2) );
}

void Orbit::explicit_rk6(void (Orbit::*func)(Vec*,double,Vec*),double* tspan, Vec X0, double h, double* Y, double* T,int MAX_ITER)

{
	double t = tspan[0];
	Vec k1,k2,k3,k4,k5,k6,k7,current,temp;
	current = X0;

	write_matrix(&current,t,Y,T,0,4);

	for(int j = 1 ; j < MAX_ITER ; j++)
	{
		(this->*func)(&k1,t,&current);
		k1 = k1.scalar(h);

		temp = current + k1;
		(this->*func)(&k2,t+h,&temp);
		k2 = k2.scalar(h);

		temp = k1.scalar(3) + k2;
		temp = current + temp.scalar(static_cast<double>(1)/8);
		(this->*func)(&k3,t+h/2,&temp);
		k3 = k3.scalar(h);

		temp = k1.scalar(8) + k2.scalar(2) + k3.scalar(8);
		temp = current + temp.scalar(static_cast<double>(1)/27);
		(this->*func)(&k4,t+2*h/3,&temp);
		k4 = k4.scalar(h);

		temp = k1.scalar( (3*(3*sqrt(21) - 7) )) + k2.scalar( -8*(7-sqrt(21)) ) + k3.scalar( 48*(7 - sqrt(21)) ) + k4.scalar( -3*(21 - sqrt(21)) );
		temp = current + temp.scalar(static_cast<double>(1)/392);
		(this->*func)(&k5,t + (7 - sqrt(21))*h/14,&temp);
		k5 = k5.scalar(h);

		temp = k1.scalar( -5*(231 + 51*sqrt(21)) ) + k2.scalar( -40*(7 + sqrt(21)) ) + k3.scalar( -320*sqrt(21) ) + k4.scalar( 3*(21 + 121*sqrt(21)) ) + k5.scalar( 392*(6 + sqrt(21)) );
		temp = current + temp.scalar(static_cast<double>(1)/1960);
		(this->*func)(&k6,t + (7 + sqrt(21))*h/14,&temp);
		
		k6 = k6.scalar(h);
		temp = k1.scalar( 15*(22+7*sqrt(21)) ) + k2.scalar( 120 ) + k3.scalar( 40*( 7*sqrt(21) - 5 ) ) + k4.scalar( -63*(3*sqrt(21) - 2) ) + k5.scalar( - 14*(49+9*sqrt(21)) ) + k6.scalar( 70*( 7 - sqrt(21)) );
		temp = current + temp.scalar(static_cast<double>(1)/180);
		(this->*func)(&k7,t + h,&temp);
		k7 = k7.scalar(h);
		
		temp = k1.scalar(9) + k3.scalar(64) + k5.scalar(49) + k6.scalar(49) + k7.scalar(9);
		current = current + temp.scalar(static_cast<double>(1)/180);
		t = t+h;

		write_matrix(&current,t,Y,T,j,4);
	}
}

void Orbit::monte_carlo(int nsample)
{
	Vec init_con;
	init_con.add_element( X0[0] );
	init_con.add_element( X0[1] );
	init_con.add_element( V0[0] );
	init_con.add_element( V0[1] );

	double* tspan = new double[2];
	tspan[0] = 0;
	tspan[1] = 4*M_PI*sqrt(pow(a,3)/mu);

	double h = tspan[1] / 1000;

	void (Orbit::*func)(Vec*,double,Vec*) = &Orbit::dynamic_2d;

	int MAX_ITER = static_cast<int>( tspan[1] / h );
	double* Y = new double[4*MAX_ITER];
	double* T = new double[MAX_ITER];

	std::default_random_engine generator;
	std::normal_distribution<double> distribution1(X0[0],1000.0);
	std::normal_distribution<double> distribution2(X0[1],1000.0);
	std::normal_distribution<double> distribution3(V0[0],1.0);
	std::normal_distribution<double> distribution4(V0[1],1.0);

	double* results = new double[4*nsample];

	for(int i = 0 ; i < nsample ; i++)
	{
		this->progress_bar(i,nsample);
		init_con.clear();
		init_con.add_element( distribution1(generator) );
		init_con.add_element( distribution2(generator) );
		init_con.add_element( distribution3(generator) );
		init_con.add_element( distribution4(generator) );

		explicit_rk6(func, tspan, init_con, h, Y, T, MAX_ITER);
		results[i*4] = Y[(MAX_ITER-1)*4];
		results[i*4+1] = Y[(MAX_ITER-1)*4+1];
		results[i*4+2] = Y[(MAX_ITER-1)*4+2];
		results[i*4+3] = Y[(MAX_ITER-1)*4+3];
	}

	char file_name[] = "monte_carlo.txt";
	this->write_matrix_file(results,nsample,4,file_name);


	delete[] Y;
	delete[] T;
	delete[] results;
}

void Orbit::SC_gaussian(int N, int d)
{
	double* tspan = new double[2];
	tspan[0] = 0;
	tspan[1] = 4*M_PI*sqrt(pow(a,3)/mu);

	double h = tspan[1] / 500;

	void (Orbit::*func)(Vec*,double,Vec*) = &Orbit::dynamic_2d;

	int MAX_ITER = static_cast<int>( tspan[1] / h );
	double* Y = new double[d*MAX_ITER];
	double* T = new double[MAX_ITER];

	vector< vector<int> > degrees = Basis::compute_mixed_expansions(N,d);

	Basis basis;
	basis.initialize_hermite(N+1);

	vector<double> sample;
	vector<double> weights;
	
	basis.get_quadrature_points( &sample , &weights );

	vector<double> lambda;
	Polynomial* pol;
	double temp;

	for(int i = 0 ; i < N+1 ; i++)
	{
		pol = basis.get_element(i);
		temp = 0;
		for(int j = 0 ; j < sample.size() ; j++)
		{
			temp += pol->evaluate( sample[j] ) * pol->evaluate( sample[j] ) * weights[j];
		}
		lambda.push_back(temp);
	}

	vector< vector<double> > u;
	vector<double> temp1;
	vector< vector<double> > X;
	vector<double> temp2;
	vector< vector<int> > sample_index;
	vector<int> temp3;
	Vec init_con;

	for(int i = 0 ; i < sample.size() ; i++)
	{
		this->progress_bar(i,sample.size());
		for(int j = 0 ; j < sample.size() ; j++)
		{
			for(int k = 0 ; k < sample.size() ; k++)
			{
				for(int l = 0 ; l < sample.size() ; l++)
				{
					init_con.clear();
					init_con.add_element( 1000.0*sample[i] + X0[0] );
					init_con.add_element( 1000.0*sample[j] + X0[1] );
					init_con.add_element( sample[k] + V0[0] );
					init_con.add_element( sample[l] + V0[1] );
				
					explicit_rk6(func, tspan, init_con, h, Y, T, MAX_ITER);
					temp1.clear();
					temp1.push_back( Y[(MAX_ITER-1)*4] );
					temp1.push_back( Y[(MAX_ITER-1)*4+1] );
					u.push_back( temp1 );

					temp2.clear();
					temp2.push_back( sample[i] );
					temp2.push_back( sample[j] );
					temp2.push_back( sample[k] );
					temp2.push_back( sample[l] );
					X.push_back( temp2 );

					temp3.clear();
					temp3.push_back(i);
					temp3.push_back(j);
					temp3.push_back(k);
					temp3.push_back(l);
					sample_index.push_back( temp3 );
				}
			}
		}
	}

	vector<double> image_expansion_x;
	vector<double> image_expansion_y;
	double u_hat;
	double phi;
	double l;
	
	for(int k = 0 ; k < degrees.size() ; k++)
	{
		u_hat = 0;
		for(int j = 0 ; j < u.size() ; j++)
		{
			phi = this->compute_hyper_phi( &basis , degrees[k] , X[j] );
			u_hat += (u[j])[0] * phi * compute_hyper_weight( weights , sample_index[j] );
		}
		l = compute_hyper_weight( lambda , degrees[k] );
		u_hat = u_hat / l;
		image_expansion_x.push_back( u_hat );
	}

	for(int k = 0 ; k < degrees.size() ; k++)
	{
		u_hat = 0;
		for(int j = 0 ; j < u.size() ; j++)
		{
			phi = this->compute_hyper_phi( &basis , degrees[k] , X[j] );
			u_hat += (u[j])[1] * phi * compute_hyper_weight( weights , sample_index[j] );
		}
		l = compute_hyper_weight( lambda , degrees[k] );
		u_hat = u_hat / l;
		image_expansion_y.push_back( u_hat );
	}

	char file_name[] = "./results/image_expansion_x.txt";
	write_1dvectord_file(image_expansion_x,file_name);
	char file_name2[] = "./results/image_expansion_y.txt";
	write_1dvectord_file(image_expansion_y,file_name2);


	cout << "mean : " << image_expansion_x[0] << endl;

	double var = 0;
	for(int i = 1 ; i < degrees.size() ; i++)
		var += compute_hyper_lambda( lambda , degrees[i] ) * image_expansion_x[i] * image_expansion_x[i];
	cout << "var : " << var << endl;

	delete[] Y;
	delete[] T;
}

double Orbit::compute_hyper_phi(Basis* basis, vector<int> degrees , vector<double> sample)
{
	double result = 1;
	Polynomial* pol;
	for(int i = 0 ; i < degrees.size() ; i++)
	{
		pol = basis->get_element( degrees[i] );
		result *= pol->evaluate( sample[i] );
	}
	return result;
}

double Orbit::compute_hyper_lambda(vector<double> lambda, vector<int> degrees)
{
	double result = 1;
	for(int i = 0 ; i < degrees.size() ; i++)
		result *= lambda[ degrees[i] ];
	return result;
}

double Orbit::compute_hyper_weight(vector<double> weights, vector<int> sample_index)
{
	double result = 1;
	for(int i = 0 ; i < sample_index.size() ; i++)
		result *= weights[ sample_index[ i ] ];
	return result;
}


void Orbit::progress_bar(int current ,int total)
{
	cout << '\r';
	int progress = (20*(current+1)) / total;
	int leftover = 20 - progress;

	cout << '[';
	for(int i = 0 ; i < progress ; i++)
		cout << '#';
	for(int i = 0 ; i < leftover ; i++)
		cout << ' ';
	cout << "] " << (100*(current+1)) / total << '%';

	if( current == total-1)
		cout << endl;

}

void Orbit::write_matrix(Vec* y, double t, double* Y, double* T, int index, int d)
{
	for(int i = 0 ; i < d ; i++)
		Y[index*d + i] = y->get_element(i);
	T[index] = t;
}

void Orbit::write_matrix_file(double* M, int nrows, int ncols,char* file_name)
{
	ofstream file;
	file.open(file_name);
	for( int i = 0 ; i < nrows ; i++ )
	{
		for(int j = 0 ; j < ncols ; j++ )
			file << M[i*ncols + j] << " ";
		file << endl;
	}
}

void Orbit::write_2dvectord_file(vector< vector<double> > M, char* file_name)
{
	ofstream file;
	file.open(file_name);
	for( int i = 0 ; i < M.size() ; i++ )
	{
		for(int j = 0 ; j < (M[i]).size() ; j++ )
			file << (M[i])[j] << " ";
		file << endl;
	}
	file.close();	
}

void Orbit::write_1dvectord_file(vector<double> M, char* file_name)
{
	ofstream file;
	file.open(file_name);
	for( int i = 0 ; i < M.size() ; i++ )
	{
		file << M[i] << endl;
	}
	file.close();
}

void Orbit::dummy()
{
	cout << "I am dummy and I am here" << endl;
}
