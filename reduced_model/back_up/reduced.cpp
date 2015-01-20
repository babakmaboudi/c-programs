#include "reduced.h"

void reduced::sv2oe( double* pos, double* vel, double* elements )
{
	double r = sqrt( pow(pos[0],2) + pow(pos[1],2) + pow(pos[2],2) );
	double initial_speed = sqrt( pow(vel[0],2) + pow(vel[1],2) + pow(vel[2],2) );
	double xi = pow(initial_speed,2)/2 - mu/r;

	elements[0] = - mu / (2 * xi);
	double* h = new double[3];
	h[0] = pos[1]*vel[2] - vel[1]*pos[2];
	h[1] = vel[0]*pos[2] - pos[0]*vel[2];
	h[2] = pos[0]*vel[1] - vel[0]*pos[1];

	double hxy = sqrt( pow(h[0],2) + pow(h[1],2) );
	elements[2] = acos( h[2]/sqrt( pow(h[0],2) + pow(h[1],2) + pow(h[2],2) ) );
	elements[3] = asin(h[0]/hxy);

	elements[1] = sqrt( 1 + (2 * xi * (pow(h[0],2) + pow(h[1],2) + pow(h[2],2)) ) / pow(mu,2) );

	double p = (pow(h[0],2) + pow(h[1],2) + pow(h[2],2)) / mu;
	elements[5] = acos( (p-r)/(r*elements[1]) );
	
	elements[4] = asin( pos[2]/(r * sin(elements[2])) );
	elements[4] = elements[4] - elements[5];

	delete[] h;
}

void reduced::initiate_from_file(char* file_path)
{
	ifstream file;
	file.open(file_path);

	char line[100];
	
	file.getline(line,100);
	G = atof(line);
	file.getline(line,100);
	M_earth = atof(line);
	file.getline(line,100);
	ae = atof(line);
	file.getline(line,100);
	J2 = atof(line);

	mu = G*M_earth;

	file.getline(line,100);

	X0 = new double[3];
	file.getline(line,100);
	X0[0] = atof(line);
	file.getline(line,100);
	X0[1] = atof(line);
	file.getline(line,100);
	X0[2] = atof(line);
	V0 = new double[3];
	file.getline(line,100);
	V0[0] = atof(line);
	file.getline(line,100);
	V0[1] = atof(line);
	file.getline(line,100);
	V0[2] = atof(line);


	file.getline(line,100);
	X_moon = new double[3];
	file.getline(line,100);
	X_moon[0] = atof(line);
	file.getline(line,100);
	X_moon[1] = atof(line);
	file.getline(line,100);
	X_moon[2] = atof(line);
	file.getline(line,100);
	M_moon = atof(line);

	file.getline(line,100);

	u_sun = new double[3];
	file.getline(line,100);
	u_sun[0] = atof(line);
	file.getline(line,100);
	u_sun[1] = atof(line);
	file.getline(line,100);
	u_sun[2] = atof(line);
	file.getline(line,100);
	P_sun = atof(line);
	file.getline(line,100);
	Cr = atof(line);

	file.getline(line,100);

	file.getline(line,100);
	W_satellite = atof(line);

	orb_vec = new double[6];
	sv2oe( X0, V0, orb_vec );

	file.close();
}

void reduced::dynamic_3d(Vec* dy, double t, Vec* y)
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

void reduced::nonlinear(Vec* dy,double t,Vec* y)
{
	double coef = - mu / pow( pow(y->get_element(0),2) + pow(y->get_element(1),2) + pow(y->get_element(2),2) , 1.5 );
	dy->clear();
	dy->add_element( 0 );
	dy->add_element( 0 );
	dy->add_element( 0 );
	dy->add_element( coef * y->get_element(0) );
	dy->add_element( coef * y->get_element(1) );
	dy->add_element( coef * y->get_element(2) );
}

void reduced::compute_nonlinear_term(void (reduced::*func)(Vec*,double,Vec*), vector<double>* result, double t, vector<double> current)
{
	Vec state, nonlinear_state;
	result->clear();
	for(int i = 0 ; i < current.size() ; i++)
		state.add_element( current[i] );
	(this->*func)(&nonlinear_state,t,&state);
	for(int i = 0 ; i < nonlinear_state.get_size() ; i++)
		result->push_back( nonlinear_state.get_element(i) );
}

void reduced::explicit_rk6(void (reduced::*func)(Vec*,double,Vec*),double* tspan, Vec init_cond, double h, vector<vector<double> >* Y, vector<double>* T,int MAX_ITER)

{
	double t = tspan[0];
	Vec k1,k2,k3,k4,k5,k6,k7,current,temp;
	current = init_cond;

	write_matrix(&current,t,Y,T);

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

		write_matrix(&current,t,Y,T);
	}
}

void reduced::test_reduced_model()
{
	Matrix Phi,Psi,A_tilde;
	int N = 5;
	double tol = 1e-7;

	this->build_reduced_model(&Phi,N,tol);

	cout << Phi << endl;
}

void reduced::build_reduced_model(Matrix* Phi,int N, double tol)
{
	void (reduced::*func)(Vec*,double,Vec*) = &reduced::dynamic_3d;
	double* tspan = new double[2];
	tspan[0] = 0.0;
	tspan[1] = 2*M_PI*sqrt(pow(orb_vec[0],3)/mu);
	double h = tspan[1]/500;
	int MAX_ITER = static_cast<int>( tspan[1] / h );
	vector< vector<double> > Y;
	vector<double> T;

	int d = 3;		//dimension of the random space
	vector<double> linear_space;
	linear_space = this->linspace(-0.5,0.5,N);

	Grid grid;
	for(int i = 0 ; i < d ; i++)
	{
		grid.add_node( linear_space );
	}

	vector< vector<double> > random_space = grid.build_grid();

	vector< vector< vector<double> > > solution;
	Vec init_cond;
	for(int i = 0 ; i < random_space.size() ; i++)
	{
		init_cond.clear();
		init_cond.add_element( X0[0] + 100*(random_space[i])[1] );
		init_cond.add_element( X0[1] + 100*(random_space[i])[2] );
		init_cond.add_element( X0[2] + 100*(random_space[i])[3] );
		init_cond.add_element( V0[0] );
		init_cond.add_element( V0[1] );
		init_cond.add_element( V0[2] );

		Y.clear();
		T.clear();

		explicit_rk6(func,tspan,init_cond,h,&Y,&T,MAX_ITER);

		solution.push_back( Y );
	}

	vector< vector<double> > snapshots;
	vector< vector<double> > snapshots_nonlinear;
	int num_samples = 10;
	vector<int> rand_index;
	srand(1);
	for(int i = 0 ; i < num_samples ; i++)
		rand_index.push_back( rand() % MAX_ITER );
	
	vector<double> nonlinear_state;
	void (reduced::*func_nonlinear)(Vec*,double,Vec*) = &reduced::nonlinear;
	for(int i = 0 ; i < random_space.size() ; i++)
	{
		for(int j = 0 ; j < rand_index.size() ; j++)
		{
			snapshots.push_back( (solution[i])[rand_index[j]] );
			compute_nonlinear_term(func_nonlinear,&nonlinear_state,1,(solution[i])[rand_index[j]]);		// If time is required replace 1 with time
			snapshots_nonlinear.push_back( nonlinear_state );
		}
	}

	Matrix M;
	M.initiate_vector_vector(snapshots);
	M = M.transpose();
	vector< vector<double> > U,V;
	vector<double> S;
	M.svd(&U,&S,&V);

	double total_phi = 0;
	for(int i = 0 ; i < S.size() ; i++)
		total_phi += S[i];
	int iterator_phi = -1;
	double checker = 0;

	while(checker < (1-tol) )
	{
		iterator_phi++;
		checker += S[iterator_phi] / total_phi;
	}

//	Matrix Phi;
	convert_to_matrix( Phi , U , iterator_phi , 'c' );

//	for(int i = 0 ; i < snapshots_nonlinear.size() ; i++)
//	{
//		for(int j = 0 ; j < snapshots_nonlinear[i].size() ; j++)
//			cout << (snapshots_nonlinear[i])[j] << " ";
//		cout << endl;
//	}

	Matrix N_mat;
	N_mat.initiate_vector_vector(snapshots_nonlinear);
	N_mat = N_mat.transpose();
	U.clear();
	V.clear();
	S.clear();
	N_mat.svd(&U,&S,&V);

	double total_psi = 0;
	for(int i = 0 ; i < S.size() ; i++)
		total_psi += S[i];
	int iterator_psi = -1;
	checker = 0;

	while(checker < (1-tol) )
	{
		iterator_psi++;
		checker += S[iterator_psi] / total_psi;
	}

	Matrix Psi;
	convert_to_matrix( &Psi , U , iterator_psi , 'c' );
}

// pass 'c' for column order and 'r' for row order

void reduced::convert_to_matrix(Matrix* mat,vector< vector<double> > M, int sig, char order)
{
	int nrows;
	int ncols;
	vector<double> temp;
	if(order == 'r')
	{
		nrows = sig;
		ncols = M[0].size();
		for(int i = 0 ; i < nrows ; i++)
		{
			for(int j = 0 ; j < ncols ; j++)
				temp.push_back( (M[i])[j] );
		}
		mat->initiate_vector(temp,nrows,ncols);
	}
	else
	{
		ncols = sig;
		nrows = M[0].size();
		for(int i = 0 ; i < nrows ; i++)
		{
			for(int j = 0 ; j < ncols ; j++)
				temp.push_back( (M[j])[i] );
		}
		mat->initiate_vector(temp,nrows,ncols);
	}	
}

void reduced::write_matrix(Vec* y,double t,vector<vector <double> >* mat,vector<double>* time)
{
	mat->push_back( y->get_vec() );
	time->push_back( t );
}

void reduced::write_file(vector<vector<double> >* mat,char* path)
{
	ofstream file;
	file.open(path);

	for(int i = 0 ; i < mat->size() ; i++)
	{
		for(int j = 0 ; j < (*mat)[i].size() ; j++)
			file << ((*mat)[i])[j] << " ";
		file << endl;
	}
	file.close();
}

void reduced::single_sat()
{
	void (reduced::*func)(Vec*,double,Vec*) = &reduced::dynamic_3d;
	double* tspan = new double[2];
	tspan[0] = 0.0;
	tspan[1] = 2*M_PI*sqrt(pow(orb_vec[0],3)/mu);
	double h = tspan[1]/500;
	int MAX_ITER = static_cast<int>( tspan[1] / h );
	vector< vector<double> > Y;
	vector<double> T;

	Vec init_cond;
	init_cond.add_element( X0[0] );
	init_cond.add_element( X0[1] );
	init_cond.add_element( X0[2] );
	init_cond.add_element( V0[0] );
	init_cond.add_element( V0[1] );
	init_cond.add_element( V0[2] );

	explicit_rk6(func,tspan,init_cond,h,&Y,&T,MAX_ITER);

	char path[] = "output.txt";
	write_file(&Y,path);
	delete[] tspan;
}

vector<double> reduced::linspace(double min, double max, int n)
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

