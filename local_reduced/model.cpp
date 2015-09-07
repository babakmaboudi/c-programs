#include "model.h"

void Model::sv2oe( vector<double> pos, vector<double> vel, double mu, vector<double>* elements)
{
	double r = sqrt( pow(pos[0],2) + pow(pos[1],2) + pow(pos[2],2) );
	double initial_speed = sqrt( pow(vel[0],2) + pow(vel[1],2) + pow(vel[2],2) );
	double xi = pow(initial_speed,2)/2 - mu/r;

	(*elements)[0] = - mu / (2 * xi);
	double* h = new double[3];
	h[0] = pos[1]*vel[2] - vel[1]*pos[2];
	h[1] = vel[0]*pos[2] - pos[0]*vel[2];
	h[2] = pos[0]*vel[1] - vel[0]*pos[1];

	double hxy = sqrt( pow(h[0],2) + pow(h[1],2) );
	(*elements)[2] = acos( h[2]/sqrt( pow(h[0],2) + pow(h[1],2) + pow(h[2],2) ) );
	(*elements)[3] = asin(h[0]/hxy);

	(*elements)[1] = sqrt( 1 + (2 * xi * (pow(h[0],2) + pow(h[1],2) + pow(h[2],2)) ) / pow(mu,2) );

	double p = (pow(h[0],2) + pow(h[1],2) + pow(h[2],2)) / mu;
	(*elements)[5] = acos( (p-r)/(r*(*elements)[1]) );
	
	(*elements)[4] = asin( pos[2]/(r * sin((*elements)[2])) );
	(*elements)[4] = (*elements)[4] - (*elements)[5];

	delete[] h;
}


int Model::find_closest_center(Matrix vec)
{
	double dist = 0;
	double d;
	Matrix temp;
	vec = vec.tr();
	temp = vec - cluster_centers.get_row(0);
	dist = temp.norm2();
	int center = 0;
	for(int i = 1 ; i < cluster_centers.get_num_rows() ; i++)
	{
		temp = vec - cluster_centers.get_row(i);
		d = temp.norm2();
		if( d < dist )
		{
			dist = d;
			center = i;
		}
	}
	return center;
}

void Model::compute_significant_subspace(Matrix* res, Matrix* mat)
{
	Matrix U, S, V;
	mat->svd(&U,&S,&V);
	cout << S << endl;

	res->clear();
	res->zeros( mat->get_num_rows() , 2 );
	for(int i = 0 ; i < 2 ; i++)
	{
		res->set_col( i , U.get_col(i) );
	}
}

void Model::set_parameters(double mu, Parameters* p)
{
	p->mu = mu;
	p->ae = param.ae;
	p->J2 = param.J2;
	p->X0.push_back( param.X0[0] );
	p->X0.push_back( param.X0[1] );
	p->X0.push_back( param.X0[2] );
	p->V0.push_back( param.V0[0] );
	p->V0.push_back( param.V0[1] );
	p->V0.push_back( param.V0[2] );
	p->X_moon.push_back( param.X_moon[0] );
	p->X_moon.push_back( param.X_moon[1] );
	p->X_moon.push_back( param.X_moon[2] );
	p->M_moon = param.M_moon;
	p->u_sun.push_back( param.u_sun[0] );
	p->u_sun.push_back( param.u_sun[1] );
	p->u_sun.push_back( param.u_sun[2] );
	p->P_sun = param.P_sun;
	p->Cr = param.Cr;
	p->W_satellite = param.W_satellite;
	p->orb_vec = vector<double>( 6, 0.0 );
	sv2oe(p->X0, p->V0, p->mu, &p->orb_vec);
}

void Model::initiate_from_file(char* file_path)
{
	ifstream file;
	file.open(file_path);

	char line[100];
	
	file.getline(line,100);
	param.G = atof(line);
	file.getline(line,100);
	param.M_earth = atof(line);
	file.getline(line,100);
	param.ae = atof(line);
	file.getline(line,100);
	param.J2 = atof(line);

	param.mu = param.G*param.M_earth;

	file.getline(line,100);

	file.getline(line,100);
	param.X0.push_back( atof(line) );
	file.getline(line,100);
	param.X0.push_back( atof(line) );
	file.getline(line,100);
	param.X0.push_back( atof(line) );
	file.getline(line,100);
	param.V0.push_back( atof(line) );
	file.getline(line,100);
	param.V0.push_back( atof(line) );
	file.getline(line,100);
	param.V0.push_back( atof(line) );


	file.getline(line,100);
	file.getline(line,100);
	param.X_moon.push_back( atof(line) );
	file.getline(line,100);
	param.X_moon.push_back( atof(line) );
	file.getline(line,100);
	param.X_moon.push_back( atof(line) );
	file.getline(line,100);
	param.M_moon = atof(line);

	file.getline(line,100);

	file.getline(line,100);
	param.u_sun.push_back( atof(line) );
	file.getline(line,100);
	param.u_sun.push_back( atof(line) );
	file.getline(line,100);
	param.u_sun.push_back( atof(line) );
	file.getline(line,100);
	param.P_sun = atof(line);
	file.getline(line,100);
	param.Cr = atof(line);

	file.getline(line,100);

	file.getline(line,100);
	param.W_satellite = atof(line);

	param.orb_vec = vector<double>( 6, 0.0 );
	sv2oe( param.X0, param.V0, param.mu, &param.orb_vec);

	file.close();
}

double Model::compute_hamiltonian(Matrix q, Matrix p, Parameters param)
{
	double kinetic = ( p.at(0)*p.at(0) + p.at(1)*p.at(1) + p.at(2)*p.at(2) )/2;
	double potential = 0;
	vector<double> polar = get_polar_coord(&q);
	double phi_angle = polar[1];
	double r = polar[0];
	potential += -param.mu / r;
	potential += +(param.mu / r) * pow( param.ae/r , 2 ) * param.J2 * ( 3 * pow( sin(phi_angle) , 2 ) - 1 ) / 2;
	return kinetic + potential;
}

void Model::deriv(Matrix* dy, double t, Matrix* y, Parameters* param)
{
	vector<double> temp;
	double coef = - param->mu / pow( pow(y->at(0,0),2) + pow(y->at(1,0),2) + pow(y->at(2,0),2) , 1.5 );
	vector<double> polar = get_polar_coord(y);
	double phi_angle = polar[1];
	double coef2 = -param->mu * pow(param->ae,2) * param->J2 / (2 * pow(polar[0],5));
	dy->at(0,0) = y->at(3,0);
	dy->at(1,0) = y->at(4,0);
	dy->at(2,0) = y->at(5,0);
	// Earth Gravity 
	dy->at(3,0) = coef * y->at(0,0);
	// Earth Oblateness
	dy->at(3,0) += coef2 * ( -15 * pow( sin(phi_angle) , 2 ) + 3 ) * y->at(0,0);
	// Earth Gravity
	dy->at(4,0) = coef * y->at(1,0);
	// Earth Oblateness
	dy->at(4,0) += coef2 * ( -15 * pow( sin(phi_angle) , 2 ) + 3 ) * y->at(1,0);
	// Earth Gravity
	dy->at(5,0) = coef * y->at(2,0);	
	// Earth Oblateness
	dy->at(5,0) += coef2 * ( -15 * pow( sin(phi_angle) , 2 ) + 9 ) * y->at(2,0);
}

void Model::deriv_symp(Matrix* dy, double t, Matrix* y, Parameters* param,int indicator)
{
	if(indicator == 0)
	{
		dy->at(0) = y->at(0);
		dy->at(1) = y->at(1);
		dy->at(2) = y->at(2);
	}
	else
	{
		double coef = - param->mu / pow( pow(y->at(0,0),2) + pow(y->at(1,0),2) + pow(y->at(2,0),2) , 1.5 );
		vector<double> polar = get_polar_coord(y);
		double phi_angle = polar[1];
		double coef2 = -param->mu * pow(param->ae,2) * param->J2 / (2 * pow(polar[0],5));
		dy->at(0) = coef * y->at(0);
		dy->at(0) += coef2 * ( -15 * pow( sin(phi_angle) , 2 ) + 3 ) * y->at(0,0);
		dy->at(1) = coef * y->at(1);
		dy->at(1) += coef2 * ( -15 * pow( sin(phi_angle) , 2 ) + 3 ) * y->at(1,0);
		dy->at(2) = coef * y->at(2);
		dy->at(2) += coef2 * ( -15 * pow( sin(phi_angle) , 2 ) + 9 ) * y->at(2,0);
	}
}

void Model::deriv_symp_reduced(Matrix* dy, double t, Matrix* y, Parameters* param,int indicator)
{
	if(indicator == 0)
	{
		(*dy) = (*y);
	}
	else
	{
		Matrix x = Phi * (*y);
		double coef = - param->mu / pow( pow(x.at(0),2) + pow(x.at(1),2) + pow(x.at(2),2) , 1.5 );
		vector<double> polar = get_polar_coord(&x);
		double phi_angle = polar[1];
		double coef2 = -param->mu * pow(param->ae,2) * param->J2 / (2 * pow(polar[0],5));

		Matrix f;
		f.zeros(3,1);
		f.at(0) = coef * x.at(0);
		f.at(0) = f.at(0) + coef2 * ( -15 * pow( sin(phi_angle) , 2 ) + 3 ) * x.at(0);
		f.at(1) = coef * x.at(1);
		f.at(1) = f.at(1) + coef2 * ( -15 * pow( sin(phi_angle) , 2 ) + 3 ) * x.at(1);
		f.at(2) = coef * x.at(2);
		f.at(2) = f.at(2) + coef2 * ( -15 * pow( sin(phi_angle) , 2 ) + 9 ) * x.at(2);
		(*dy) = Phi.tr() * f;
	}
}

void Model::integ4_symp(void (Model::*func)(Matrix*,double,Matrix*,Parameters*,int),double* tspan, Matrix init_pos, Matrix init_mom, Parameters* param, double h, Matrix* q, Matrix* p ,Matrix* T, int MAX_ITER)
{
	double t = tspan[0];
	Matrix da,db,a,b;

	a = init_pos;
	b = init_mom;

	da.zeros( a.get_num_rows() , a.get_num_cols() );
	db.zeros( b.get_num_rows() , b.get_num_cols() );

	T->zeros(MAX_ITER,1);
	T->at(0) = 0;

	double c1,c2,c3,c4,d1,d2,d3,d4;
	c1 = 1/( 2 * (2-pow(2,1/3)) );
	c4 = c1;
	c2 = ( 1 - pow(2,1/3) )/( 2 * (2-pow(2,1/3)) ) ;
	c3 = c2;
	d1 = 1 / ( 2 - pow(2,1/3) );
	d3 = d1;
	d2 = - pow(2,1/3) / (2 - pow(2,1/3));
	d4 = 0;

	for(int j = 0 ; j < MAX_ITER ; j++)
	{
		(this->*func)(&da,t,&b,param,0);
		da = da.scalar( c1*h );
		a = a + da;

		(this->*func)(&db,t,&a,param,1);
		db = db.scalar( d1*h );
		b = b + db;

		(this->*func)(&da,t,&b,param,0);
		da = da.scalar( c2*h );
		a = a + da;

		(this->*func)(&db,t,&a,param,1);
		db = db.scalar( d2*h );
		b = b + db;

		(this->*func)(&da,t,&b,param,0);
		da = da.scalar( c3*h );
		a = a + da;

		(this->*func)(&db,t,&a,param,1);
		db = db.scalar( d3*h );
		b = b + db;

		(this->*func)(&da,t,&b,param,0);
		da = da.scalar( c4*h );
		a = a + da;

		(this->*func)(&db,t,&a,param,1);
		db = db.scalar( d4*h );
		b = b + db;

		q->add_row(a.tr());
		p->add_row(b.tr());
		T->at(j) = t;
	}
}

void Model::integ4_symp_reduced(void (Model::*func)(Matrix*,double,Matrix*,Parameters*,int),double* tspan, Matrix init_pos, Matrix init_mom, Parameters* param, double h, Matrix* q, Matrix* p ,Matrix* T, int MAX_ITER)
{
	double t = tspan[0];
	Matrix da,db,a,b;

	a = init_pos;
	b = init_mom;

	da.zeros( a.get_num_rows() , a.get_num_cols() );
	db.zeros( b.get_num_rows() , b.get_num_cols() );

	T->zeros(MAX_ITER,1);
	T->at(0) = 0;

	double c1,c2,c3,c4,d1,d2,d3,d4;
	c1 = 1/( 2 * (2-pow(2,1/3)) );
	c4 = c1;
	c2 = ( 1 - pow(2,1/3) )/( 2 * (2-pow(2,1/3)) ) ;
	c3 = c2;
	d1 = 1 / ( 2 - pow(2,1/3) );
	d3 = d1;
	d2 = - pow(2,1/3) / (2 - pow(2,1/3));
	d4 = 0;

	for(int j = 0 ; j < MAX_ITER ; j++)
	{
		Matrix full_dim = Phi * a;
		full_dim.append( Phi * b , 'r' );
		int check = find_closest_center(full_dim);


		if(check != nearest)
		{
			nearest = check;
			Phi = local_phi[nearest];
			Matrix position = full_dim.get_submat(0,3,0,1);
			Matrix momentum = full_dim.get_submat(3,6,0,1);
			a = Phi.tr() * position;
			b = Phi.tr() * momentum;

		}

		(this->*func)(&da,t,&b,param,0);
		da = da.scalar( c1*h );
		a = a + da;

		(this->*func)(&db,t,&a,param,1);
		db = db.scalar( d1*h );
		b = b + db;

		(this->*func)(&da,t,&b,param,0);
		da = da.scalar( c2*h );
		a = a + da;

		(this->*func)(&db,t,&a,param,1);
		db = db.scalar( d2*h );
		b = b + db;

		(this->*func)(&da,t,&b,param,0);
		da = da.scalar( c3*h );
		a = a + da;

		(this->*func)(&db,t,&a,param,1);
		db = db.scalar( d3*h );
		b = b + db;

		(this->*func)(&da,t,&b,param,0);
		da = da.scalar( c4*h );
		a = a + da;

		(this->*func)(&db,t,&a,param,1);
		db = db.scalar( d4*h );
		b = b + db;

//		q->add_row(a.tr());
		q->add_row(full_dim.tr());
		p->add_row(b.tr());
		T->at(j) = t;
	}
}

void Model::build_reduced_model(int N, double tol)
{
	int d = 1;
	vector<double> linear_space;
	linear_space = linspace(-0.5,0.5,N);

	Grid grid;
	for(int i = 0 ; i < d ; i++)
	{
		grid.add_node( linear_space );
	}
	vector< vector<double> > random_space = grid.build_grid();

	double* tspan = new double[2];
	vector<Parameters> param_handler;
	void (Model::*func)(Matrix*,double,Matrix*,Parameters*,int) = &Model::deriv_symp;
	vector<Matrix> solution_pos, solution_mom;
	int MAX_ITER = 500;
	for(int i=0 ; i<random_space.size() ; i++)
	{
		Parameters temp_param;
		set_parameters(param.mu + 4e12*(random_space[i])[0], &temp_param);
		tspan[0] = 0.0;
		tspan[1] = 2*M_PI*sqrt(pow(temp_param.orb_vec[0],3)/temp_param.mu);
		double h = tspan[1]/MAX_ITER;

		Matrix init_pos(3);
		Matrix init_mom(3);
		init_pos.at(0) = temp_param.X0[0];
		init_pos.at(1) = temp_param.X0[1];
		init_pos.at(2) = temp_param.X0[2];
		init_mom.at(0) = temp_param.V0[0];
		init_mom.at(1) = temp_param.V0[1];
		init_mom.at(2) = temp_param.V0[2];

		Matrix Q,P,T;

		integ4_symp(func,tspan,init_pos,init_mom,&temp_param,h,&Q,&P,&T,MAX_ITER);
		solution_pos.push_back( Q );
		solution_mom.push_back( P );
	}
	
	Matrix all_pos, all_mom, all_data;
	all_pos = solution_pos[0];
	all_mom = solution_mom[0];
	for(int i=1 ; i<solution_pos.size() ; i++)
	{
		all_pos.append(solution_pos[i],'r');
		all_mom.append(solution_mom[i],'r');
	}
	char path[] = "./data/output.txt";
	all_pos.save(path);
	all_data = all_pos;
	all_data.append(all_mom,'c');

	int num_clusters = 10;
	cluster_centers = kmeans(&all_data,num_clusters,1000);
	Matrix labels = label_vectors(&all_data,&cluster_centers);

	vector<int> IDX;
	Matrix section_pos, section_mom;
	for(int i=0 ; i<num_clusters ; i++)
	{
		IDX.clear();
		find(&IDX,&labels,i);
		section_pos.clear();
		section_mom.clear();
		for(int j=0 ; j<IDX.size() ; j++)
		{
			section_pos.add_row( all_pos.get_row( IDX[j] ) );
			section_mom.add_row( all_mom.get_row( IDX[j] ) );
		}
		Matrix U, S, V;
		Matrix mat;
		mat = section_pos;
		mat.append(section_mom,'r');
		mat = mat.tr();

		Matrix dummy;
		compute_significant_subspace(&dummy,&mat);
		local_phi.push_back(dummy);
	}
/*
	Matrix snap_pos;
	Matrix snap_mom;
	int num_samples = 20;
	vector<int> rand_index;
	for(int i = 0 ; i < num_samples ; i++)
		rand_index.push_back( rand() % MAX_ITER );

	for(int i = 0 ; i < random_space.size() ; i++)
	{
		for(int j = 0 ; j < rand_index.size() ; j++)
		{
			double idx = rand_index[j];
			snap_pos.add_row( (solution_pos[i]).get_row(idx) );
			snap_mom.add_row( (solution_mom[i]).get_row(idx) );
		}
	}

	Matrix ext_snap( 2*snap_pos.length() , 3 );
	ext_snap.set_submat(0,0,snap_pos);
	ext_snap.set_submat(snap_pos.length(),0,snap_mom);
	ext_snap = ext_snap.tr();

	Matrix U,S,V;
	ext_snap.svd(&U,&S,&V);

	phi1 = U.get_submat(0,U.get_num_rows(),0,2);
	phi2 = U.get_submat(0,U.get_num_rows(),0,2);
*/
}

void Model::test_reduced_model()
{
	int N = 5;
	double tol = 1e-10;

	this->build_reduced_model(N,tol);

	double* random_vec = new double[6];
	double r;
	for(int i = 0 ; i < 6 ; i++)
	{
		r = (double) rand() / RAND_MAX;
		random_vec[i] = r - 0.5;
	}

	Parameters P0;
	set_parameters(param.mu + random_vec[0] * 4e12,&P0);


	Matrix pos;
	pos.zeros(3,1);
	pos.at(0) = P0.X0[0];
	pos.at(1) = P0.X0[1];
	pos.at(2) = P0.X0[2];

	Matrix mom;
	mom.zeros(3,1);
	mom.at(0) = P0.V0[0];
	mom.at(1) = P0.V0[1];
	mom.at(2) = P0.V0[2];

	Matrix vec = pos;
	vec.append(mom,'r');
	nearest = find_closest_center(vec);
	Phi = local_phi[nearest];

	Matrix a = Phi.tr() * pos;
	Matrix b = Phi.tr() * mom;

	void (Model::*func)(Matrix*,double,Matrix*,Parameters*,int) = &Model::deriv_symp_reduced;
	double* tspan = new double[2];
	tspan[0] = 0.0;
	tspan[1] = 2*M_PI*sqrt(pow(P0.orb_vec[0],3)/P0.mu);
	int MAX_ITER = 500;
	double h = tspan[1]/MAX_ITER;
	Matrix Q,P;
	Matrix T;

	integ4_symp_reduced(func,tspan,a,b,&P0,h,&Q,&P,&T,MAX_ITER);

//	Q = (Phi * Q.tr()).tr();
//	P = (Phi * P.tr()).tr();
	char path[] = "./data/output.txt";
	Q.save(path);
	
	Matrix energy_reduced(Q.length());
	for(int i=0 ; i<Q.length() ; i++)
	{
		Matrix q_vec = Q.get_row(i);
		Matrix p_vec = P.get_row(i);
		energy_reduced.at(i) = compute_hamiltonian(q_vec.tr(),p_vec.tr(),P0);
	}
	
	cout << endl << "------" << endl;
	void (Model::*func2)(Matrix*,double,Matrix*,Parameters*,int) = &Model::deriv_symp;
	Matrix Q_exact,P_exact;
	T.clear();
	integ4_symp(func2,tspan,pos,mom,&P0,h,&Q_exact,&P_exact,&T,MAX_ITER);
	char path2[] = "./data/output2.txt";
	Q_exact.save(path2);

	Matrix energy(Q_exact.length());
	for(int i=0 ; i<Q_exact.length() ; i++)
	{
		Matrix q_vec = Q_exact.get_row(i);
		Matrix p_vec = P_exact.get_row(i);
		energy.at(i) = compute_hamiltonian(q_vec.tr(),p_vec.tr(),P0);
	}

	Matrix temp = energy - energy_reduced;
	cout << temp << endl;
}

void Model::single_sat()
{
	Matrix init_pos(3);
	Matrix init_mom(3);
	init_pos.at(0) = param.X0[0];
	init_pos.at(1) = param.X0[1];
	init_pos.at(2) = param.X0[2];
	init_mom.at(0) = param.V0[0];
	init_mom.at(1) = param.V0[1];
	init_mom.at(2) = param.V0[2];

	Parameters P0;
	set_parameters(param.mu -0.0 * 4e12,&P0);

	void (Model::*func)(Matrix*,double,Matrix*,Parameters*,int) = &Model::deriv_symp;
	double* tspan = new double[2];
	tspan[0] = 0.0;
	tspan[1] = 2*M_PI*sqrt(pow(P0.orb_vec[0],3)/P0.mu);
	double h = tspan[1]/500;
	int MAX_ITER = static_cast<int>( tspan[1] / h );
	Matrix Q;
	Matrix P;
	Matrix T;

	integ4_symp(func,tspan,init_pos,init_mom,&P0,h,&Q,&P,&T,MAX_ITER);

	for(int i=0 ; i<Q.length() ; i++)
	{
		Matrix q_vec = Q.get_row(i);
		Matrix p_vec = P.get_row(i);
		cout << compute_hamiltonian(q_vec.tr(),p_vec.tr(),P0) << endl;
	}

//	char path[] = "./data/output.txt";
//	Q.save(path);

	delete[] tspan;
}
