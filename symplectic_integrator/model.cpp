#include "model.h"

void Model::sv2oe( double* pos, double* vel, double mu, double* elements)
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

void Model::set_parameters(double mu, Parameters* p)
{
	p->mu = mu;
	p->ae = param.ae;
	p->X0 = new double[3];
	p->X0[0] = param.X0[0];
	p->X0[1] = param.X0[1];
	p->X0[2] = param.X0[2];
	p->V0 = new double[3];
	p->V0[0] = param.V0[0];
	p->V0[1] = param.V0[1];
	p->V0[2] = param.V0[2];
	p->X_moon = new double[3];
	p->X_moon[0] = param.X_moon[0];
	p->X_moon[1] = param.X_moon[1];
	p->X_moon[2] = param.X_moon[2];
	p->M_moon = param.M_moon;
	p->u_sun = new double[3];
	p->u_sun[0] = param.u_sun[0];
	p->u_sun[1] = param.u_sun[1];
	p->u_sun[2] = param.u_sun[2];
	p->P_sun = param.P_sun;
	p->Cr = param.Cr;
	p->W_satellite = param.W_satellite;
	p->orb_vec = new double[6];
	sv2oe(p->X0, p->V0, p->mu, p->orb_vec);
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
	sv2oe( param.X0, param.V0, param.mu, param.orb_vec);

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

void Model::single_sat()
{
	void (Model::*func)(Matrix*,double,Matrix*,Parameters*,int) = &Model::deriv_symp;
//	void (Model::*func)(Matrix*,double,Matrix*,Parameters*) = &Model::deriv;
	double* tspan = new double[2];
	tspan[0] = 0.0;
	tspan[1] = 2*M_PI*sqrt(pow(param.orb_vec[0],3)/param.mu);
	double h = tspan[1]/500;
	int MAX_ITER = static_cast<int>( tspan[1] / h );
	Matrix Y;
	Matrix T;

	Matrix init_cond(6);
	init_cond.at(0) = param.X0[0];
	init_cond.at(1) = param.X0[1];
	init_cond.at(2) = param.X0[2];
	init_cond.at(3) = param.V0[0];
	init_cond.at(4) = param.V0[1];
	init_cond.at(5) = param.V0[2];

	euler_symp(func,tspan,init_cond,&param,h,&Y,&T,MAX_ITER);
//	euler(func,tspan,init_cond,&param,h,&Y,&T,MAX_ITER);

	char path[] = "./data/output.txt";
	Y.save(path);

	delete[] tspan;
}

void Model::test()
{
	Matrix init_cond(6);
	init_cond.at(0) = param.X0[0];
	init_cond.at(1) = param.X0[1];
	init_cond.at(2) = param.X0[2];
	init_cond.at(3) = param.V0[0];
	init_cond.at(4) = param.V0[1];
	init_cond.at(5) = param.V0[2];

	Parameters P0;
	set_parameters(param.mu -0.5 * 4e12,&P0);

	void (Model::*func)(Matrix*,double,Matrix*,Parameters*,int) = &Model::deriv_symp;
	double* tspan = new double[2];
	tspan[0] = 0.0;
	tspan[1] = 2*M_PI*sqrt(pow(P0.orb_vec[0],3)/P0.mu);
	cout << 2*M_PI*sqrt(pow(P0.orb_vec[0],3)/P0.mu << endl;
			getchar();
	double h = tspan[1]/500;
	int MAX_ITER = static_cast<int>( tspan[1] / h );
	Matrix Y;
	Matrix T;

	euler_symp(func,tspan,init_cond,&P0,h,&Y,&T,MAX_ITER);
	char path[] = "./data/output.txt";
	Y.save(path);

	delete[] tspan;
}
