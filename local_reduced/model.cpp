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

void Model::BRM_params()
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
	linear_space = linspace(-0.5,0.5,5);

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

	vector< vector<double> > all_data;
	for(int i = 0 ; i < solution.size() ; i++)
	{
		for(int j = 0 ; j < solution[i].size() ; j++)
			all_data.push_back( (solution[i])[j] );
	}
	Matrix mat;
	mat.initiate_vector_vector(all_data);

	Matrix cluster_centers = kmeans(&mat,10,1000);
	Matrix labels = label_vectors(&mat,&cluster_centers);

	char path1[] = "./data/data.txt";
	char path2[] = "./data/centers.txt";
	char path3[] = "./data/labels.txt";
	mat.save(path1);
	cluster_centers.save(path2);
	labels.save(path3);
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
