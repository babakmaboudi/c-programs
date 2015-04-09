#include "wave.h"

Wave::Wave()
{
	param.L = 1;
	param.C = 0.1;
}

Matrix Wave::initial_condition(int N)
{
	Matrix X;
	X.linspace(0,1,N);

	Matrix init_cond( X.length() , 1 );

	for(int i=0 ; i<X.length() ; i++)
	{
		double s = X.at(i);
		s = 10*abs(s - 1.0/2.0);
		if(s <= 1)
			init_cond.at(i) = 1.0 - 3.0/2.0*pow(s,2) + 3.0/4.0*pow(s,3);
		else if(s <= 2)
			init_cond.at(i) = 1.0/4.0*pow(2.0-s,3);
		else
			init_cond.at(i) = 0;
	}
	return init_cond;
}

void Wave::compute_significant_subspace(Matrix* res, Matrix* mat, double tol)
{
	Matrix U, S, V;
	mat->svd(&U,&S,&V);

	double total = 0;
	for(int i = 0 ; i < S.length() ; i++)
		total += S.at(i);
	int iterator = -1;
	double checker = 0;
	while( checker < (1-tol) )
	{
		iterator++;
		checker += S.at(iterator) / total;
	}
	res->clear();
	res->zeros( mat->get_num_rows() , iterator );
	for(int i = 0 ; i < iterator ; i++)
	{
		res->set_col( i , U.get_col(i) );
	}
}

void Wave::deriv_symp(Matrix* dy, double t, Matrix* y, Parameters* P0, int indicator)
{
	if(indicator == 0)
	{
		(*dy) = (*y);
	}
	else
	{
		(*dy) = Dxx*(*y);

	}
}

void Wave::integ2_symp(void (Wave::*func)(Matrix*,double,Matrix*,Parameters*,int),double* tspan, Matrix init_pos, Matrix init_mom, Parameters* P0, double h, Matrix* Q, Matrix* P, Matrix* T, int MAX_ITER)
{
	double t = tspan[0];
	Matrix temp,q,p;
	temp.zeros(3,1);
	q = init_pos;
	p = init_mom;

	Q->clear();
	P->clear();
	T->zeros(MAX_ITER,1);
	T->at(0) = 0;

	Q->add_row( q.tr() );
	P->add_row( p.tr() );

	for(int j = 0 ; j < MAX_ITER ; j++)
	{
		(this->*func)(&temp,t,&p,P0,0);
		temp = temp.scalar( h/2 );
		q = q + temp;

		(this->*func)(&temp,t,&q,P0,1);
		temp = temp.scalar(h/2);
		p = p + temp;

		(this->*func)(&temp,t,&p,P0,0);
		temp = temp.scalar( h/2 );
		q = q + temp;

		(this->*func)(&temp,t,&q,P0,1);
		temp = temp.scalar(h/2);
		p = p + temp;

		Q->add_row( q.tr() );
		P->add_row( p.tr() );

		T->at(j) = t;
	}
}

void Wave::build_reduced_model(int N, double tol)
{
	Dxx.zeros(N,N);

	for(int i=1 ; i<N-1 ; i++)
	{
		Dxx.at(i,i-1) = 1;
		Dxx.at(i,i) = -2;
		Dxx.at(i,i+1) = 1;
	}
	Dxx.at(0,0) = -2;
	Dxx.at(0,1) = 1;
	Dxx.at(0,N-1) = 1;

	Dxx.at(N-1,0) = 1;
	Dxx.at(N-1,N-1) = -2;
	Dxx.at(N-1,N-2) = 1;

	double dx =param.L/N;
	Dxx = Dxx/(dx*dx);
	Dxx = Dxx*(pow(param.C,2));

	Matrix init_pos = initial_condition(N);
	Matrix init_mom(N,1);

	void (Wave::*func)(Matrix*,double,Matrix*,Parameters*,int) = &Wave::deriv_symp;
	double* tspan = new double[2];
	tspan[0] = 0.0;
	tspan[1] = 50;
	double h = 0.01;
	int MAX_ITER = static_cast<int>( tspan[1] / h );
	Matrix Q;
	Matrix P;
	Matrix T;
	
	integ2_symp(func,tspan,init_pos,init_mom,&param,h,&Q,&P,&T,MAX_ITER);

	Matrix snap_pos;
	Matrix snap_mom;
	for(int i=0 ; i<100 ; i++)
		snap_pos.add_row( Q.get_row( i*50 ) );
	for(int i=0 ; i<100 ; i++)
		snap_mom.add_row( P.get_row( i*50 ) );

	
	Matrix snapshots = snap_pos;
	int nrows, ncols;
	snapshots.append(snap_mom,'r');
	snapshots = snapshots.tr();

	compute_significant_subspace(&Phi,&snapshots,tol);
}

void Wave::test_reduced_model()
{
	int N = 500;
	double tol = 1e-3;
	build_reduced_model(N,tol);

	Dxx.zeros(N,N);

	for(int i=1 ; i<N-1 ; i++)
	{
		Dxx.at(i,i-1) = 1;
		Dxx.at(i,i) = -2;
		Dxx.at(i,i+1) = 1;
	}
	Dxx.at(0,0) = -2;
	Dxx.at(0,1) = 1;
	Dxx.at(0,N-1) = 1;

	Dxx.at(N-1,0) = 1;
	Dxx.at(N-1,N-1) = -2;
	Dxx.at(N-1,N-2) = 1;

	double dx =param.L/N;
	Dxx = Dxx/(dx*dx);
	Dxx = Dxx*(pow(param.C,2));

	Dxx = Phi.tr() * Dxx * Phi;

	Matrix init_pos = initial_condition(N);
	Matrix init_mom(N,1);

	Matrix init_pos_red = Phi.tr() * init_pos;
	Matrix init_mom_red = Phi.tr() * init_mom;

	void (Wave::*func)(Matrix*,double,Matrix*,Parameters*,int) = &Wave::deriv_symp;
	double* tspan = new double[2];
	tspan[0] = 0.0;
	tspan[1] = 50;
	double h = 0.01;
	int MAX_ITER = static_cast<int>( tspan[1] / h );
	Matrix Q;
	Matrix P;
	Matrix T;
	
	integ2_symp(func,tspan,init_pos_red,init_mom_red,&param,h,&Q,&P,&T,MAX_ITER);
	
	Matrix result = Phi*Q.tr();
	result = result.tr();
	char path[] = "./data/out.txt";
	result.save(path);
}

void Wave::solver(int N)
{
	Dxx.zeros(N,N);

	for(int i=1 ; i<N-1 ; i++)
	{
		Dxx.at(i,i-1) = 1;
		Dxx.at(i,i) = -2;
		Dxx.at(i,i+1) = 1;
	}
	Dxx.at(0,0) = -2;
	Dxx.at(0,1) = 1;
	Dxx.at(0,N-1) = 1;

	Dxx.at(N-1,0) = 1;
	Dxx.at(N-1,N-1) = -2;
	Dxx.at(N-1,N-2) = 1;

	double dx =param.L/N;
	Dxx = Dxx/(dx*dx);
	Dxx = Dxx*(pow(param.C,2));	

	Matrix init_pos = initial_condition(N);
	Matrix init_mom(N,1);

	void (Wave::*func)(Matrix*,double,Matrix*,Parameters*,int) = &Wave::deriv_symp;
	double* tspan = new double[2];
	tspan[0] = 0.0;
	tspan[1] = 50;
	double h = 0.01;
	int MAX_ITER = static_cast<int>( tspan[1] / h );
	Matrix Q;
	Matrix P;
	Matrix T;
	
	integ2_symp(func,tspan,init_pos,init_mom,&param,h,&Q,&P,&T,MAX_ITER);
	char path[] = "./data/out.txt";
	Q.save(path);

	delete[] tspan;
}
