#include "model.h"

void Model::sv2oe( double* pos, double* vel, double* elements )
{
	double r = sqrt( pow(pos[0],2) + pow(pos[1],2) + pow(pos[2],2) );
	double initial_speed = sqrt( pow(vel[0],2) + pow(vel[1],2) + pow(vel[2],2) );
	double xi = pow(initial_speed,2)/2 - param.mu/r;

	elements[0] = - param.mu / (2 * xi);
	double* h = new double[3];
	h[0] = pos[1]*vel[2] - vel[1]*pos[2];
	h[1] = vel[0]*pos[2] - pos[0]*vel[2];
	h[2] = pos[0]*vel[1] - vel[0]*pos[1];

	double hxy = sqrt( pow(h[0],2) + pow(h[1],2) );
	elements[2] = acos( h[2]/sqrt( pow(h[0],2) + pow(h[1],2) + pow(h[2],2) ) );
	elements[3] = asin(h[0]/hxy);

	elements[1] = sqrt( 1 + (2 * xi * (pow(h[0],2) + pow(h[1],2) + pow(h[2],2)) ) / pow(param.mu,2) );

	double p = (pow(h[0],2) + pow(h[1],2) + pow(h[2],2)) / param.mu;
	elements[5] = acos( (p-r)/(r*elements[1]) );
	
	elements[4] = asin( pos[2]/(r * sin(elements[2])) );
	elements[4] = elements[4] - elements[5];

	delete[] h;
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

	param.X0 = new double[3];
	file.getline(line,100);
	param.X0[0] = atof(line);
	file.getline(line,100);
	param.X0[1] = atof(line);
	file.getline(line,100);
	param.X0[2] = atof(line);
	param.V0 = new double[3];
	file.getline(line,100);
	param.V0[0] = atof(line);
	file.getline(line,100);
	param.V0[1] = atof(line);
	file.getline(line,100);
	param.V0[2] = atof(line);


	file.getline(line,100);
	param.X_moon = new double[3];
	file.getline(line,100);
	param.X_moon[0] = atof(line);
	file.getline(line,100);
	param.X_moon[1] = atof(line);
	file.getline(line,100);
	param.X_moon[2] = atof(line);
	file.getline(line,100);
	param.M_moon = atof(line);

	file.getline(line,100);

	param.u_sun = new double[3];
	file.getline(line,100);
	param.u_sun[0] = atof(line);
	file.getline(line,100);
	param.u_sun[1] = atof(line);
	file.getline(line,100);
	param.u_sun[2] = atof(line);
	file.getline(line,100);
	param.P_sun = atof(line);
	file.getline(line,100);
	param.Cr = atof(line);

	file.getline(line,100);

	file.getline(line,100);
	param.W_satellite = atof(line);

	param.orb_vec = new double[6];
	sv2oe( param.X0, param.V0, param.orb_vec );

	file.close();
}
/*
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
//	dy->at(3,0) += coef2 * ( -15 * pow( sin(phi) , 2 ) + 3 ) * y->at(0,0);
	// Earth Gravity
	dy->at(4,0) = coef * y->at(1,0);
	// Earth Oblateness
//	dy->at(4,0) += coef2 * ( -15 * pow( sin(phi) , 2 ) + 3 ) * y->at(1,0);
	// Earth Gravity
	dy->at(5,0) = coef * y->at(2,0);	
	// Earth Oblateness
//	dy->at(5,0) += coef2 * ( -15 * pow( sin(phi) , 2 ) + 9 ) * y->at(2,0);
}
*/

void Model::deriv_symp(Matrix* dy, double t, Matrix* y, Parameters* param,int indicator)
{
	vector<double> temp;
	double coef = - param->mu / pow( pow(y->at(0,0),2) + pow(y->at(1,0),2) + pow(y->at(2,0),2) , 1.5 );
	vector<double> polar = get_polar_coord(y);
	double phi = polar[1];
	double coef2 = -param->mu * pow(param->ae,2) * param->J2 / (2 * pow(polar[0],5));
	if(indicator == 0)
	{
		dy->at(0,0) = y->at(3,0);
		dy->at(1,0) = y->at(4,0);
		dy->at(2,0) = y->at(5,0);
		dy->at(3,0) = 0;
		dy->at(4,0) = 0;
		dy->at(5,0) = 0;	
	}
	else
	{
		dy->at(0,0) = 0;
		dy->at(1,0) = 0;
		dy->at(2,0) = 0;
		dy->at(3,0) = coef * y->at(0,0);
		dy->at(3,0) += coef2 * ( -15 * pow( sin(phi) , 2 ) + 3 ) * y->at(0,0);
		dy->at(4,0) = coef * y->at(1,0);
		dy->at(4,0) += coef2 * ( -15 * pow( sin(phi) , 2 ) + 3 ) * y->at(1,0);
		dy->at(5,0) = coef * y->at(2,0);
		dy->at(5,0) += coef2 * ( -15 * pow( sin(phi) , 2 ) + 9 ) * y->at(2,0);
	}
}

void Model::deriv_symp_test(Matrix* dy, double t, Matrix* y, Parameters* param,int indicator)
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

void Model::euler(void (Model::*func)(Matrix*,double,Matrix*,Parameters*),double* tspan, Matrix init_cond, Parameters* param, double h, Matrix* Y, Matrix* T, int MAX_ITER)
{
	double t = tspan[0];
	Matrix current,temp;
	current.zeros(6,1);
	temp.zeros(6,1);

	current = init_cond;
	Y->clear();
	T->zeros(MAX_ITER,1);

	T->at(0) = 0;

	Y->add_row( current.tr() );

	for(int j = 0 ; j < MAX_ITER ; j++)
	{
		(this->*func)(&temp,t,&current,param);
		temp = temp.scalar(h);
		current = temp + current;
		t = t + h;

		Y->add_row( current.tr() );
		T->at(j) = t;
	}
}

void Model::euler_symp(void (Model::*func)(Matrix*,double,Matrix*,Parameters*,int),double* tspan, Matrix init_cond, Parameters* param, double h, Matrix* Y, Matrix* T, int MAX_ITER)
{
	double t = tspan[0];
	Matrix current,temp;
	current.zeros(6,1);
	temp.zeros(6,1);

	current = init_cond;
	Y->clear();
	T->zeros(MAX_ITER,1);

	T->at(0) = 0;

	Y->add_row( current.tr() );

	for(int j = 0 ; j < MAX_ITER ; j++)
	{
		(this->*func)(&temp,t,&current,param,0);
		temp = temp.scalar(h);
		current = temp + current;

		(this->*func)(&temp,t,&current,param,1);
		temp = temp.scalar(h);
		current = temp + current;
		t = t + h;

		Y->add_row( current.tr() );
		T->at(j) = t;
	}
}

void Model::integ2_symp(void (Model::*func)(Matrix*,double,Matrix*,Parameters*,int),double* tspan, Matrix init_cond, Parameters* param, double h, Matrix* Y, Matrix* T, int MAX_ITER)
{
	double t = tspan[0];
	Matrix current,temp,q,p;
	current.zeros(6,1);
	temp.zeros(6,1);
	q.zeros(6,1);
	p.zeros(6,1);

	q.at(0) = init_cond.at(0);
	q.at(1) = init_cond.at(1);
	q.at(2) = init_cond.at(2);
	p.at(3) = init_cond.at(3);
	p.at(4) = init_cond.at(4);
	p.at(5) = init_cond.at(5);

	current = init_cond;
	Y->clear();
	T->zeros(MAX_ITER,1);

	T->at(0) = 0;

	Y->add_row( current.tr() );

	for(int j = 0 ; j < MAX_ITER ; j++)
	{
		(this->*func)(&temp,t,&current,param,0);
		temp = temp.scalar( h/2 );
		current = current + temp;

		(this->*func)(&temp,t,&current,param,1);
		temp = temp.scalar(h/2);
		current = current + temp;

		(this->*func)(&temp,t,&current,param,0);
		temp = temp.scalar( h/2 );
		current = current + temp;

		(this->*func)(&temp,t,&current,param,1);
		temp = temp.scalar(h/2);
		current = current + temp;

		Y->add_row( current.tr() );
		T->at(j) = t;
	}
}

void Model::integ4_symp(void (Model::*func)(Matrix*,double,Matrix*,Parameters*,int),double* tspan, Matrix init_cond, Parameters* param, double h, Matrix* Y, Matrix* T, int MAX_ITER)
{
	double t = tspan[0];
	Matrix current,temp,q,p;
	current.zeros(6,1);
	temp.zeros(6,1);
	q.zeros(6,1);
	p.zeros(6,1);

	q.at(0) = init_cond.at(0);
	q.at(1) = init_cond.at(1);
	q.at(2) = init_cond.at(2);
	p.at(3) = init_cond.at(3);
	p.at(4) = init_cond.at(4);
	p.at(5) = init_cond.at(5);

	current = init_cond;
	Y->clear();
	T->zeros(MAX_ITER,1);

	T->at(0) = 0;

	Y->add_row( current.tr() );

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
		(this->*func)(&temp,t,&current,param,0);
		temp = temp.scalar( c1*h );
		current = current + temp;

		(this->*func)(&temp,t,&current,param,1);
		temp = temp.scalar( d1*h );
		current = current + temp;

		(this->*func)(&temp,t,&current,param,0);
		temp = temp.scalar( c2*h );
		current = current + temp;

		(this->*func)(&temp,t,&current,param,1);
		temp = temp.scalar( d2*h );
		current = current + temp;

		(this->*func)(&temp,t,&current,param,0);
		temp = temp.scalar( c3*h );
		current = current + temp;

		(this->*func)(&temp,t,&current,param,1);
		temp = temp.scalar( d3*h );
		current = current + temp;

		(this->*func)(&temp,t,&current,param,0);
		temp = temp.scalar( c4*h );
		current = current + temp;

		(this->*func)(&temp,t,&current,param,1);
		temp = temp.scalar( d4*h );
		current = current + temp;

		Y->add_row( current.tr() );
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
			
		Matrix temp = a.tr();
		q->add_row(temp);

		q->add_row(a.tr());
		p->add_row(b.tr());
		T->at(j) = t;
	}
}

void Model::build_reduced_model(int N, double tol)
{
	void (Model::*func)(Matrix*,double,Matrix*,Parameters*,int) = &Model::deriv_symp;
	double* tspan = new double[2];
	tspan[0] = 0.0;
	tspan[1] = 2*M_PI*sqrt(pow(param.orb_vec[0],3)/param.mu);
	double h = tspan[1]/500;
	int MAX_ITER = static_cast<int>( tspan[1] / h );
	Matrix Y;
	Matrix T;

	int d = 1;		//dimension of the random space
	vector<double> linear_space;
	linear_space = linspace(-0.5,0.5,N);

	Grid grid;
	for(int i = 0 ; i < d ; i++)
	{
		grid.add_node( linear_space );
	}

	vector< vector<double> > random_space = grid.build_grid();

	vector<double> ic;
	ic.push_back( param.X0[0] );
	ic.push_back( param.X0[1] );
	ic.push_back( param.X0[2] );
	ic.push_back( param.V0[0] );
	ic.push_back( param.V0[1] );
	ic.push_back( param.V0[2] );
	Matrix init_cond;
	init_cond.initiate_vector(ic,6,1);

	Parameters rand_param;
	vector<Matrix> solution;
	vector<Parameters> param_handler;

	for(int i = 0 ; i < random_space.size() ; i++)
	{
		Y.clear();
		T.clear();

		rand_param.mu = param.mu + 4e12*(random_space[i])[0];
		rand_param.J2 = param.J2;
		rand_param.ae = param.ae;

		integ4_symp(func,tspan,init_cond,&rand_param,h,&Y,&T,MAX_ITER);

		solution.push_back( Y );
		param_handler.push_back(param);
	}

	Matrix snap_lin;
	Matrix snap_nonlin;
	int num_samples = 20;
	vector<int> rand_index;
	srand(1);
	for(int i = 0 ; i < num_samples ; i++)
		rand_index.push_back( rand() % MAX_ITER );

	for(int i = 0 ; i < random_space.size() ; i++)
	{
		for(int j = 0 ; j < rand_index.size() ; j++)
		{
			double idx;
			idx = rand_index[j];
			snap_lin.add_row( (solution[i]).get_submat(idx,idx+1,0,3) );
			snap_nonlin.add_row( (solution[i]).get_submat(idx,idx+1,3,6) );
		}
	}
	Matrix U,S,V;
	snap_lin = snap_lin.tr();
	snap_lin.svd(&U,&S,&V);
	
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

	snap_nonlin = snap_nonlin.tr();
	snap_nonlin.svd(&U,&S,&V);

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
//	phi2 = U;
}

void Model::test_reduced_model()
{
	void (Model::*func)(Matrix*,double,Matrix*,Parameters*,int) = &Model::deriv_symp_reduced;
	double* tspan = new double[2];
	tspan[0] = 0.0;
	tspan[1] = 2*M_PI*sqrt(pow(param.orb_vec[0],3)/param.mu);
	double h = tspan[1]/500;
	int MAX_ITER = static_cast<int>( tspan[1] / h );
	Matrix Q,P;
	Matrix T;

	int N = 5;
	double tol = 1e-3;

	this->build_reduced_model(N,tol);

	L = phi1.tr() * phi2;

	srand(2);
	double* random_vec = new double[6];
	double r;
	for(int i = 0 ; i < 6 ; i++)
	{
		r = (double) rand() / RAND_MAX;
		random_vec[i] = r - 0.5;
	}
	Parameters rand_param;
	rand_param.mu = param.mu + 4e12*random_vec[0];
	rand_param.J2 = param.J2;
	rand_param.ae = param.ae;

	Matrix pos;
	pos.zeros(3,1);
	pos.at(0) = param.X0[0];
	pos.at(1) = param.X0[1];
	pos.at(2) = param.X0[2];

	Matrix mom;
	mom.zeros(3,1);
	mom.at(0) = param.V0[0];
	mom.at(1) = param.V0[1];
	mom.at(2) = param.V0[2];

	Matrix a = phi1.tr() * pos;
	Matrix b = phi2.tr() * mom;

	integ4_symp_reduced(func,tspan,a,b,&rand_param,h,&Q,&P,&T,MAX_ITER);

	Matrix temp;
	Matrix result;
	for(int i = 0 ; i < MAX_ITER ; i++)
	{
		temp = (phi1 * (Q.get_row(i)).tr()).tr();
		result.add_row(temp);
	}

	char path[] = "./data/output.txt";
	result.save(path);

	void (Model::*func2)(Matrix*,double,Matrix*,Parameters*,int) = &Model::deriv_symp_test;
	Matrix Q_exact,P_exact;
	T.clear();
	integ4_symp_reduced(func2,tspan,pos,mom,&rand_param,h,&Q_exact,&P_exact,&T,MAX_ITER);
	char path2[] = "./data/output2.txt";
	Q_exact.save(path2);

	delete[] tspan;
	cout << phi1 <<endl;
	cout << phi2 << endl;
}

void Model::single_sat()
{
	void (Model::*func)(Matrix*,double,Matrix*,Parameters*,int) = &Model::deriv_symp_test;
	double* tspan = new double[2];
	tspan[0] = 0.0;
	tspan[1] = 2*M_PI*sqrt(pow(param.orb_vec[0],3)/param.mu);
	double h = tspan[1]/500;
	int MAX_ITER = static_cast<int>( tspan[1] / h );
	Matrix Q,P;
	Matrix T;

	Matrix pos;
	pos.zeros(3,1);
	pos.at(0) = param.X0[0];
	pos.at(1) = param.X0[1];
	pos.at(2) = param.X0[2];

	Matrix mom;
	mom.zeros(3,1);
	mom.at(0) = param.V0[0];
	mom.at(1) = param.V0[1];
	mom.at(2) = param.V0[2];

	double r = -0.5;
	cout << r << endl;
	param.mu = param.mu + 4e12*r;
	integ4_symp_reduced(func,tspan,pos,mom,&param,h,&Q,&P,&T,MAX_ITER);

	char path[] = "./data/output.txt";
	Q.save(path);

	delete[] tspan;
}
