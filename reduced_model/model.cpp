#include "model.h"

void Model::sv2oe( double* pos, double* vel, double* elements )
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

void Model::discrete_empirical_interpolation(Matrix* P, Matrix* mat)
{
	double dummy;
	int arg_max;
	Matrix col = mat->get_col(0);
	
	arg_opt(&dummy,&arg_max,&col,'b');

	int nrows,ncols;
	mat->get_size(&nrows,&ncols);

	Matrix U,c,z,r;
	P->zeros(nrows,1);
	z.zeros(nrows,1);
	P->at(arg_max,0) = 1;
	U = col;

	for(int i = 1 ; i < ncols ; i++)
	{
		col = mat->get_col(i);
		c = (P->tr() * U).inv() * P->tr() * col;
		r = col - U*c;
		arg_opt(&dummy,&arg_max,&r,'b');
		U.add_col(col);
		P->add_col(z);
		P->at(arg_max,i) = 1;
	}
}

void Model::initiate_from_file(char* file_path)
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

//	'c' for column-wise order and 'r' for row-wise order
void Model::convert_to_matrix(Matrix* mat,vector< vector<double> > M, int sig, char order)
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

void Model::write_file(vector<vector<double> >* mat,char* path)
{
	ofstream file;
	file.open(path);

	for(int i = 0 ; i < mat->size() ; i++)
	{
		for(int j = 0 ; j < (*mat)[i].size() ; j++)
			file << std::setprecision(16) << ((*mat)[i])[j] << " ";
		file << endl;
	}
	file.close();
}

void Model::write_matrix(Matrix* y, double t, vector<vector<double> >* mat, vector<double>* time)
{
	mat->push_back( y->get_matrix() );
	time->push_back( t );
}

void Model::dynamic_3d(Matrix* dy, double t, Matrix* y)
{
	vector<double> temp;
	double coef = - mu / pow( pow(y->at(0,0),2) + pow(y->at(1,0),2) + pow(y->at(2,0),2) , 1.5 );
	dy->at(0,0) = y->at(3,0);
	dy->at(1,0) = y->at(4,0);
	dy->at(2,0) = y->at(5,0);
	dy->at(3,0) = coef * y->at(0,0);
	dy->at(4,0) = coef * y->at(1,0);
	dy->at(5,0) = coef * y->at(2,0);
}

void Model::nonlinear(Matrix* dy, double t, Matrix* y)
{
	vector<double> temp;
	double coef = - mu / pow( pow(y->at(0,0),2) + pow(y->at(1,0),2) + pow(y->at(2,0),2) , 1.5 );
	dy->at(0,0) = 0;
	dy->at(1,0) = 0;
	dy->at(2,0) = 0;
	dy->at(3,0) = coef * y->at(0,0);
	dy->at(4,0) = coef * y->at(1,0);
	dy->at(5,0) = coef * y->at(2,0);
}

void Model::dynamic_reduced(Matrix* dy, double t, Matrix* y)
{
	Matrix u,non_lin;
	u = Phi * (*y);

	non_lin.zeros(6,1);
	nonlinear(&non_lin,t,&u);
	(*dy) = A_tilde*(*y) + Phi.tr()*non_lin;
}

void Model::compute_nonlinear_term(void (Model::*func)(Matrix*,double,Matrix*),vector<double>* result, double t, vector<double> current)
{
	Matrix state, nonlinear_state;
	state.zeros(6,1);
	nonlinear_state.zeros(6,1);
	result->clear();

	for(int i = 0 ; i < current.size() ; i++)
		state.at(i,0) = current[i];
	(this->*func)(&nonlinear_state,t,&state);
	int nrows,ncols;
	nonlinear_state.get_size(&nrows,&ncols);
	for(int i = 0 ; i < nrows ; i++)
		result->push_back( nonlinear_state.at(i,0) );
}

void Model::explicit_rk6(void (Model::*func)(Matrix*,double,Matrix*),double* tspan, Matrix init_cond, double h, vector<vector<double> >* Y, vector<double>* T,int MAX_ITER)

{
	double t = tspan[0];
	Matrix k1,k2,k3,k4,k5,k6,k7,current,temp;
	k1.zeros(6,1);
	k2.zeros(6,1);
	k3.zeros(6,1);
	k4.zeros(6,1);
	k5.zeros(6,1);
	k6.zeros(6,1);
	k7.zeros(6,1);
	current.zeros(6,1);
	temp.zeros(6,1);

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

void Model::build_reduced_model(int N, double tol)
{
	void (Model::*func)(Matrix*,double,Matrix*) = &Model::dynamic_3d;
	double* tspan = new double[2];
	tspan[0] = 0.0;
	tspan[1] = 2*M_PI*sqrt(pow(orb_vec[0],3)/mu);
	double h = tspan[1]/500;
	int MAX_ITER = static_cast<int>( tspan[1] / h );
	vector< vector<double> > Y;
	vector<double> T;

	int d = 3;		//dimension of the random space
	vector<double> linear_space;
	linear_space = linspace(-0.5,0.5,N);

	Grid grid;
	for(int i = 0 ; i < d ; i++)
	{
		grid.add_node( linear_space );
	}

	vector< vector<double> > random_space = grid.build_grid();

	vector< vector< vector<double> > > solution;
	Matrix init_cond;
	init_cond.zeros(6,1);
	for(int i = 0 ; i < random_space.size() ; i++)
	{
		init_cond.at(0,0) = X0[0] + 1000*(random_space[i])[1] ;
		init_cond.at(1,0) = X0[1] + 1000*(random_space[i])[2] ;
		init_cond.at(2,0) = X0[2] + 1000*(random_space[i])[3] ;
		init_cond.at(3,0) = V0[0] ;
		init_cond.at(4,0) = V0[1] ;
		init_cond.at(5,0) = V0[2] ;

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
//	srand((time(NULL)%10)*313123141);
	for(int i = 0 ; i < num_samples ; i++)
		rand_index.push_back( rand() % MAX_ITER );
	
	vector<double> nonlinear_state;
	void (Model::*func_nonlinear)(Matrix*,double,Matrix*) = &Model::nonlinear;
	for(int i = 0 ; i < random_space.size() ; i++)
	{
		for(int j = 0 ; j < rand_index.size() ; j++)
		{
			snapshots.push_back( (solution[i])[rand_index[j]] );
			compute_nonlinear_term(func_nonlinear,&nonlinear_state,1,(solution[i])[rand_index[j]]);		// If time is required replace 1 with time
			snapshots_nonlinear.push_back( nonlinear_state );
		}
	}

	char source[] = "./data/all_result.txt";
	write_file(&snapshots,source);

	Matrix M;
	M.initiate_vector_vector(snapshots);
	M = M.tr();
	vector< vector<double> > U,V;
	vector<double> S;
	M.svd(&U,&S,&V);

	char path[] = "./data/U.txt";
	this->write_file(&U,path);

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

	convert_to_matrix( &Phi , U , iterator_phi , 'c' );

	Matrix N_mat;
	N_mat.initiate_vector_vector(snapshots_nonlinear);
	N_mat = N_mat.tr();
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

	convert_to_matrix( &Psi , U , iterator_psi , 'c' );

	discrete_empirical_interpolation(&P,&Psi);
}

void Model::test_reduced_model()
{
	void (Model::*func)(Matrix*,double,Matrix*) = &Model::dynamic_reduced;
	double* tspan = new double[2];
	tspan[0] = 0.0;
	tspan[1] = 2*M_PI*sqrt(pow(orb_vec[0],3)/mu);
	double h = tspan[1]/500;
	int MAX_ITER = static_cast<int>( tspan[1] / h );
	vector< vector<double> > Y;
	vector<double> T;

	int N = 5;
	double tol = 1e-10;

	A.zeros(6,6);
	A.at(0,3) = 1;
	A.at(1,4) = 1;
	A.at(2,5) = 1;

	this->build_reduced_model(N,tol);
	A_tilde = Phi.tr() * A * Phi;

	cout << Phi << endl;
	cout << A_tilde << endl;	

	srand(2);
	double* random_vec = new double[6];
	double r;
	for(int i = 0 ; i < 6 ; i++)
	{
		r = (double) rand() / RAND_MAX;
		random_vec[i] = r - 0.5;
	}

	vector<double> ic;
	ic.push_back( X0[0] + 1000*random_vec[0]);
	ic.push_back( X0[1] + 1000*random_vec[1]);
	ic.push_back( X0[2] + 1000*random_vec[2]);
	ic.push_back( V0[0] );
	ic.push_back( V0[1] );
	ic.push_back( V0[2] );
	Matrix init_cond,init_cond_red;
	init_cond.initiate_vector(ic,6,1);
	init_cond_red = Phi.tr() * init_cond;

// ----------------
//	char path1[] = "./data/phi.txt";
//	char path2[] = "./data/init_cond.txt";
//	char path3[] = "./data/init_cond_red.txt";
//	Phi.save(path1);
//	init_cond.save(path2);
//	init_cond_red.save(path3);
//	Matrix test;
//	test = init_cond - (Phi * init_cond_red);
//
//	cout << "This is the difference :" << endl << test << endl;
// ----------------

	explicit_rk6(func,tspan,init_cond_red,h,&Y,&T,MAX_ITER);

	Matrix result, result_reduced, temp, vec;
	result_reduced.initiate_vector_vector(Y);
	for(int i = 0 ; i < MAX_ITER ; i++)
	{
		temp.clear();
		temp = ( result_reduced.get_row(i) ).tr();

		vec.clear();
		vec = Phi * temp;
		result.add_row(vec.tr());
	}

	char path[] = "./data/result.txt";
	result.save(path);

	void (Model::*func2)(Matrix*,double,Matrix*) = &Model::dynamic_3d;
	vector< vector<double> > exact;
	T.clear();
	explicit_rk6(func2,tspan,init_cond,h,&exact,&T,MAX_ITER);

	char path2[] = "./data/result_exact.txt";
	write_file(&exact,path2);
	
}

void Model::single_sat(double* rand_vec, int len)
{
	void (Model::*func)(Matrix*,double,Matrix*) = &Model::dynamic_3d;
	double* tspan = new double[2];
	tspan[0] = 0.0;
	tspan[1] = 2*M_PI*sqrt(pow(orb_vec[0],3)/mu);
	double h = tspan[1]/500;
	int MAX_ITER = static_cast<int>( tspan[1] / h );
	vector< vector<double> > Y;
	vector<double> T;

	vector<double> ic;
	ic.push_back( X0[0] );
	ic.push_back( X0[1] );
	ic.push_back( X0[2] );
	ic.push_back( V0[0] );
	ic.push_back( V0[1] );
	ic.push_back( V0[2] );
	Matrix init_cond;
	init_cond.initiate_vector(ic,6,1);

	explicit_rk6(func,tspan,init_cond,h,&Y,&T,MAX_ITER);

	char path[] = "./data/output.txt";
	write_file(&Y,path);
	delete[] tspan;	
}
