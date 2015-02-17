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

void Model::set_parameters(double mu, Parameters* p)
{
	p->mu = mu;
	p->ae = param.ae;
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

void Model::deriv(Matrix* dy, double t, Matrix* y, Parameters* param)
{
	vector<double> temp;
	double coef = - param->mu / pow( pow(y->at(0,0),2) + pow(y->at(1,0),2) + pow(y->at(2,0),2) , 1.5 );
	vector<double> polar = get_polar_coord(y);
	double phi = polar[1];
	double coef2 = -param->mu * pow(param->ae,2) * param->J2 / (2 * pow(polar[0],5));
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
		double phi = polar[1];
		double coef2 = -param->mu * pow(param->ae,2) * param->J2 / (2 * pow(polar[0],5));
		dy->at(0) = coef * y->at(0);
		dy->at(0) += coef2 * ( -15 * pow( sin(phi) , 2 ) + 3 ) * y->at(0,0);
		dy->at(1) = coef * y->at(1);
		dy->at(1) += coef2 * ( -15 * pow( sin(phi) , 2 ) + 3 ) * y->at(1,0);
		dy->at(2) = coef * y->at(2);
		dy->at(2) += coef2 * ( -15 * pow( sin(phi) , 2 ) + 9 ) * y->at(2,0);
	}
}

void Model::deriv_symp_reduced(Matrix* dy, double t, Matrix* y, Parameters* param,int indicator)
{
	if(indicator == 0)
	{
		(*dy) = L * (*y);
	}
	else
	{
		Matrix x = phi1 * (*y);
		double coef = - param->mu / pow( pow(x.at(0),2) + pow(x.at(1),2) + pow(x.at(2),2) , 1.5 );
		vector<double> polar = get_polar_coord(y);
		double phi = polar[1];
		double coef2 = -param->mu * pow(param->ae,2) * param->J2 / (2 * pow(polar[0],5));

		Matrix f;
		f.zeros(3,1);
		f.at(0) = coef * x.at(0);
		f.at(0) = f.at(0) + coef2 * ( -15 * pow( sin(phi) , 2 ) + 3 ) * x.at(0);
		f.at(1) = coef * x.at(1);
		f.at(1) = f.at(1) + coef2 * ( -15 * pow( sin(phi) , 2 ) + 3 ) * x.at(1);
		f.at(2) = coef * x.at(2);
		f.at(2) = f.at(2) + coef2 * ( -15 * pow( sin(phi) , 2 ) + 9 ) * x.at(2);
		(*dy) = phi2.tr() * f;
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
	Matrix U,S,V;
	snap_pos = snap_pos.tr();
	snap_pos.svd(&U,&S,&V);
	
	double total = 0;
	for(int i = 0 ; i < S.length() ; i++)
		total += S.at(i);
	int iterator = -1;
	double checker = 0;

	while(checker < (1-tol) )
	{
		iterator++;
		checker += S.at(iterator) / total;
	}

	phi1 = U.get_submat(0,U.get_num_rows(),0,iterator+1);
//	phi1 = U;

	snap_mom = snap_mom.tr();
	snap_mom.svd(&U,&S,&V);

	total = 0;
	for(int i = 0 ; i < S.length() ; i++)
		total += S.at(i);
	iterator = -1;
	checker = 0;

	while(checker < (1-tol) )
	{
		iterator++;
		checker += S.at(iterator) / total;
	}

	phi2 = U.get_submat(0,U.get_num_rows(),0,iterator+1);
}

void Model::test_reduced_model()
{
	int N = 5;
	double tol = 1e-10;

	this->build_reduced_model(N,tol);
	L = phi1.tr() * phi2;

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

	Matrix a = phi1.tr() * pos;
	Matrix b = phi2.tr() * mom;

	cout << phi1 << endl << phi2 << endl;

	void (Model::*func)(Matrix*,double,Matrix*,Parameters*,int) = &Model::deriv_symp_reduced;
	double* tspan = new double[2];
	tspan[0] = 0.0;
	tspan[1] = 2*M_PI*sqrt(pow(P0.orb_vec[0],3)/P0.mu);
	int MAX_ITER = 500;
	double h = tspan[1]/MAX_ITER;
	Matrix Q,P;
	Matrix T;

	integ4_symp(func,tspan,a,b,&P0,h,&Q,&P,&T,MAX_ITER);

	Matrix temp;
	Matrix result;
	for(int i = 0 ; i < MAX_ITER ; i++)
	{
		temp = (phi1 * (Q.get_row(i)).tr()).tr();
		result.add_row(temp);
	}
	char path[] = "./data/output.txt";
	result.save(path);

	void (Model::*func2)(Matrix*,double,Matrix*,Parameters*,int) = &Model::deriv_symp;
	Matrix Q_exact,P_exact;
	T.clear();
	integ4_symp(func2,tspan,pos,mom,&P0,h,&Q_exact,&P_exact,&T,MAX_ITER);
	char path2[] = "./data/output2.txt";
	Q_exact.save(path2);
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
	char path[] = "./data/output.txt";
	Q.save(path);

	delete[] tspan;
}
