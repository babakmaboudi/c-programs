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
	P_list.push_back(arg_max);

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
		P_list.push_back(arg_max);
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

void Model::dynamic_3d_param(Matrix* dy, double t, Matrix* y, Parameters* param)
{
	vector<double> temp;
	double coef = - (*param).mu / pow( pow(y->at(0,0),2) + pow(y->at(1,0),2) + pow(y->at(2,0),2) , 1.5 );
	vector<double> polar = get_polar_coord(y);
	double phi = polar[1];
	double coef2 = -param->mu * pow(ae,2) * J2 / (2 * pow(polar[0],5));
	dy->at(0,0) = y->at(3,0);
	dy->at(1,0) = y->at(4,0);
	dy->at(2,0) = y->at(5,0);
	// Earth Gravity 
	dy->at(3,0) = coef * y->at(0,0);
	// Earth Oblateness
	dy->at(3,0) += coef2 * ( -15 * pow( sin(phi) , 2 ) + 3 ) * y->at(0,0);
	// Earth Gravity
	dy->at(4,0) = coef * y->at(1,0);
	// Earth Oblateness
	dy->at(4,0) += coef2 * ( -15 * pow( sin(phi) , 2 ) + 3 ) * y->at(1,0);
	// Earth Gravity
	dy->at(5,0) = coef * y->at(2,0);	
	// Earth Oblateness
	dy->at(5,0) += coef2 * ( -15 * pow( sin(phi) , 2 ) + 9 ) * y->at(2,0);
}

void Model::nonlinear(Matrix* dy, double t, Matrix* y,Parameters* param)
{
	vector<double> temp;
	double coef = - param->mu / pow( pow(y->at(0,0),2) + pow(y->at(1,0),2) + pow(y->at(2,0),2) , 1.5 );
	vector<double> polar = get_polar_coord(y);
	double phi = polar[1];
	double coef2 = -param->mu * pow(ae,2) * J2 / (2*pow(polar[0],5));
	dy->at(0,0) = 0;
	dy->at(1,0) = 0;
	dy->at(2,0) = 0;
	// Earth Gravity
	dy->at(3,0) = coef * y->at(0,0);
	// Earth Oblateness
	dy->at(3,0) += coef2 * (-15 * pow( sin(phi) , 2 ) + 3 ) * y->at(0,0);
	// Earth Gravity
	dy->at(4,0) = coef * y->at(1,0);
	// Earth Oblateness
	dy->at(4,0) += coef2 * ( -15 * pow( sin(phi) , 2 ) + 3 ) * y->at(1,0);
	// Earth Gravity
	dy->at(5,0) = coef * y->at(2,0);
	// Earth Oblateness
	dy->at(5,0) += coef2 * ( -15 * pow( sin(phi) , 2 ) + 9 ) * y->at(2,0);
}

void Model::dynamic_reduced_param(Matrix* dy, double t, Matrix* y,Parameters* param)
{
	Matrix u,non_lin;
	u = Phi * (*y);
/*
	non_lin.zeros(6,1);
	nonlinear(&non_lin,t,&u,param);
	(*dy) = A_tilde*(*y) + Phi.tr()*non_lin;
*/

	double coef = - param->mu / pow( pow(u.at(0,0),2) + pow(u.at(1,0),2) + pow(u.at(2,0),2) , 1.5 );
	vector<double> polar = get_polar_coord(&u);
	double phi = polar[1];
	double coef2 = -param->mu * pow(ae,2) * J2 / (2*pow(polar[0],5));
	non_lin.zeros(P_list.size() , 1);
	for(int i = 0 ; i < P_list.size() ; i++)
	{
		switch( P_list[i] )
		{
			case 0:
				non_lin.at( i , 0 ) = 0;
				break;
			case 1:
				non_lin.at( i , 0 ) = 0;
				break;
			case 2:
				non_lin.at( i , 0 ) = 0;
				break;
			case 3:
				// Earth gravity
				non_lin.at( i , 0 ) = coef * u.at(0,0);
				// Earth Oblateness
				non_lin.at( i , 0 ) += coef2 * ( -15 * pow( sin(phi) , 2 ) + 3 )*u.at(0,0);
				break;
			case 4:
				// Earth gravity
				non_lin.at( i , 0 ) = coef * u.at(1,0);
				// Earth Oblateness
				non_lin.at( i , 0 ) += coef2 * ( -15 * pow( sin(phi) , 2 ) + 3 )*u.at(1,0);
				break;
			case 5:
				// Earth Gravity
				non_lin.at( i , 0 ) = coef * u.at(2,0);
				// Earth Oblateness
				non_lin.at( i , 0 ) += coef2 * ( -15 * pow( sin(phi) , 2 ) + 9 )*u.at(2,0);
				break;
		}
	}
	(*dy) = A_tilde*(*y) + F_tilde * non_lin;

}

void Model::compute_nonlinear_term(void (Model::*func)(Matrix*,double,Matrix*,Parameters*),vector<double>* result, double t, vector<double> current, Parameters* param)
{
	Matrix state, nonlinear_state;
	state.zeros(6,1);
	nonlinear_state.zeros(6,1);
	result->clear();

	for(int i = 0 ; i < current.size() ; i++)
		state.at(i,0) = current[i];
	(this->*func)(&nonlinear_state,t,&state,param);
	int nrows,ncols;
	nonlinear_state.get_size(&nrows,&ncols);
	for(int i = 0 ; i < nrows ; i++)
		result->push_back( nonlinear_state.at(i,0) );
}

void Model::explicit_rk6(void (Model::*func)(Matrix*,double,Matrix*,Parameters*),double* tspan, Matrix init_cond, Parameters* param, double h, vector<vector<double> >* Y, vector<double>* T,int MAX_ITER)

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
		(this->*func)(&k1,t,&current,param);
		k1 = k1.scalar(h);
		
		temp = current + k1;
		(this->*func)(&k2,t+h,&temp,param);
		k2 = k2.scalar(h);

		temp = k1.scalar(3) + k2;
		temp = current + temp.scalar(static_cast<double>(1)/8);
		(this->*func)(&k3,t+h/2,&temp,param);
		k3 = k3.scalar(h);

		temp = k1.scalar(8) + k2.scalar(2) + k3.scalar(8);
		temp = current + temp.scalar(static_cast<double>(1)/27);
		(this->*func)(&k4,t+2*h/3,&temp,param);
		k4 = k4.scalar(h);

		temp = k1.scalar( (3*(3*sqrt(21) - 7) )) + k2.scalar( -8*(7-sqrt(21)) ) + k3.scalar( 48*(7 - sqrt(21)) ) + k4.scalar( -3*(21 - sqrt(21)) );
		temp = current + temp.scalar(static_cast<double>(1)/392);
		(this->*func)(&k5,t + (7 - sqrt(21))*h/14,&temp,param);
		k5 = k5.scalar(h);

		temp = k1.scalar( -5*(231 + 51*sqrt(21)) ) + k2.scalar( -40*(7 + sqrt(21)) ) + k3.scalar( -320*sqrt(21) ) + k4.scalar( 3*(21 + 121*sqrt(21)) ) + k5.scalar( 392*(6 + sqrt(21)) );
		temp = current + temp.scalar(static_cast<double>(1)/1960);
		(this->*func)(&k6,t + (7 + sqrt(21))*h/14,&temp,param);
	
		k6 = k6.scalar(h);
		temp = k1.scalar( 15*(22+7*sqrt(21)) ) + k2.scalar( 120 ) + k3.scalar( 40*( 7*sqrt(21) - 5 ) ) + k4.scalar( -63*(3*sqrt(21) - 2) ) + k5.scalar( - 14*(49+9*sqrt(21)) ) + k6.scalar( 70*( 7 - sqrt(21)) );
		temp = current + temp.scalar(static_cast<double>(1)/180);
		(this->*func)(&k7,t + h,&temp,param);
		k7 = k7.scalar(h);

		temp = k1.scalar(9) + k3.scalar(64) + k5.scalar(49) + k6.scalar(49) + k7.scalar(9);
		current = current + temp.scalar(static_cast<double>(1)/180);
		t = t+h;
		
		write_matrix(&current,t,Y,T);
	}
}

void Model::BRM_params(int N, double tol)
{
	void (Model::*func)(Matrix*,double,Matrix*,Parameters*) = &Model::dynamic_3d_param;
	double* tspan = new double[2];
	tspan[0] = 0.0;
	tspan[1] = 2*M_PI*sqrt(pow(orb_vec[0],3)/mu);
	double h = tspan[1]/500;
	int MAX_ITER = static_cast<int>( tspan[1] / h );
	vector< vector<double> > Y;
	vector<double> T;

	int d = 1;		//dimension of the random space
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

	init_cond.at(0,0) = X0[0];
	init_cond.at(1,0) = X0[1];
	init_cond.at(2,0) = X0[2];
	init_cond.at(3,0) = V0[0];
	init_cond.at(4,0) = V0[1];
	init_cond.at(5,0) = V0[2];

	Parameters param;

	for(int i = 0 ; i < random_space.size() ; i++)
	{
		Y.clear();
		T.clear();

		param.mu = mu + 4e12*(random_space[i])[0];

		explicit_rk6(func,tspan,init_cond,&param,h,&Y,&T,MAX_ITER);

		solution.push_back( Y );
		param_handler.push_back(param);
	}

	vector< vector<double> > snapshots;
	vector< vector<double> > snapshots_nonlinear;
	int num_samples = 300;
	vector<int> rand_index;
	srand(1);
//	srand((time(NULL)%10)*313123141);
	for(int i = 0 ; i < num_samples ; i++)
		rand_index.push_back( rand() % MAX_ITER );
	
	vector<double> nonlinear_state;
	void (Model::*func_nonlinear)(Matrix*,double,Matrix*,Parameters*) = &Model::nonlinear;
	for(int i = 0 ; i < random_space.size() ; i++)
	{
		for(int j = 0 ; j < rand_index.size() ; j++)
		{
			snapshots.push_back( (solution[i])[rand_index[j]] );
			compute_nonlinear_term(func_nonlinear,&nonlinear_state,1,(solution[i])[rand_index[j]],&(param_handler[ i ]) );		// If time is required replace 1 with time
			snapshots_nonlinear.push_back( nonlinear_state );
		}
	}

	Matrix M;
	M.initiate_vector_vector(snapshots);
	M = M.tr();
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

	iterator_phi = 5;
	convert_to_matrix( &Phi , U , iterator_phi , 'c' );

	for(int i = 0 ; i < S.size() ; i++)
		cout << S[i] << endl;

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

	iterator_psi = 3;
	convert_to_matrix( &Psi , U , iterator_psi , 'c' );

	discrete_empirical_interpolation(&P,&Psi);
}

void Model::TRM_params()
{
	void (Model::*func)(Matrix*,double,Matrix*,Parameters*) = &Model::dynamic_reduced_param;
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

	this->BRM_params(N,tol);
	A_tilde = Phi.tr() * A * Phi;
	F_tilde = Phi.tr() * Psi * ( P.tr() * Psi ).inv();

	srand(2);
	double* random_vec = new double[6];
	double r;
	for(int i = 0 ; i < 6 ; i++)
	{
		r = (double) rand() / RAND_MAX;
		random_vec[i] = r - 0.5;
	}

	vector<double> ic;
	ic.push_back( X0[0] );
	ic.push_back( X0[1] );
	ic.push_back( X0[2] );
	ic.push_back( V0[0] );
	ic.push_back( V0[1] );
	ic.push_back( V0[2] );
	Matrix init_cond,init_cond_red;
	init_cond.initiate_vector(ic,6,1);
	init_cond_red = Phi.tr() * init_cond;

	Parameters param;

	param.mu = mu + 4e12*random_vec[0];
//	param.mu = mu;

	explicit_rk6(func,tspan,init_cond_red,&param,h,&Y,&T,MAX_ITER);

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

	char path1[] = "./data/result.txt";
	result.save(path1);

	void (Model::*func2)(Matrix*,double,Matrix*,Parameters*) = &Model::dynamic_3d_param;
	vector< vector<double> > exact;
	T.clear();
	explicit_rk6(func2,tspan,init_cond,&param,h,&exact,&T,MAX_ITER);

	char path2[] = "./data/result_exact.txt";
	write_file(&exact,path2);
}

void Model::single_sat()
{
	void (Model::*func)(Matrix*,double,Matrix*,Parameters*) = &Model::dynamic_3d_param;
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

	Parameters param;
	param.mu = mu;

	explicit_rk6(func,tspan,init_cond,&param,h,&Y,&T,MAX_ITER);

	char path[] = "./data/output.txt";
	write_file(&Y,path);
	delete[] tspan;	
}
