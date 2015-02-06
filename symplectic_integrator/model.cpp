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
		dy->at(4,0) = coef * y->at(1,0);
		dy->at(5,0) = coef * y->at(2,0);
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

void Model::single_sat()
{
	void (Model::*func)(Matrix*,double,Matrix*,Parameters*,int) = &Model::deriv_symp;
	double* tspan = new double[2];
	tspan[0] = 0.0;
	tspan[1] = 2*M_PI*sqrt(pow(param.orb_vec[0],3)/param.mu);
	double h = tspan[1]/500;
	int MAX_ITER = static_cast<int>( tspan[1] / h );
	Matrix Y;
	Matrix T;

	vector<double> ic;
	ic.push_back( param.X0[0] );
	ic.push_back( param.X0[1] );
	ic.push_back( param.X0[2] );
	ic.push_back( param.V0[0] );
	ic.push_back( param.V0[1] );
	ic.push_back( param.V0[2] );
	Matrix init_cond;
	init_cond.initiate_vector(ic,6,1);

	euler_symp(func,tspan,init_cond,&param,h,&Y,&T,MAX_ITER);

	char path[] = "./data/output.txt";
	Y.save(path);
	delete[] tspan;	
}
