#include "model.h"

void Model::initiate()
{
	L = 50;
	N = 20;
	dx = L/N;
	x0 = 10;

	T = 25;
	dt = 0.0125;
	MAX_ITER = T/dt;
	v = 0.2;

	Dxx.zeros(N,N);
	bd.zeros(N,1);
}

void Model::initial_condition(Matrix* q, Matrix* p)
{
	q->zeros(N,1);
	p->zeros(N,1);

	Matrix X = linspace( 0,L,N );
	for(int i=0 ; i<N ; i++)
	{
		q->at(i) = 4*atan(exp((X.at(i)-x0)/sqrt(1-v*v)));
		p->at(i) = -4*v*exp((X.at(i)-x0)/sqrt(1-v*v))/( sqrt(1-v*v)*(1+exp(2*(X.at(i)-x0)/sqrt(1-v*v))));
	}
}

void Model::sine_Gordon()
{
	for(int i=0 ; i<N ; i++)
		Dxx.at(i,i) = -2;
	for(int i=0 ; i<N-1 ; i++)
	{
		Dxx.at(i,i+1) = 1;
		Dxx.at(i+1,i) = 1;
	}
	Dxx /= (dx*dx);
	bd.at(N-1) = 2*M_PI/(dx*dx);

	Matrix q,p;
	initial_condition(&q,&p);
	void (Model::*func)(Matrix*,double,Matrix*,int) = &Model::deriv_symp;
	double* tspan = new double[2];
	tspan[0] = 0;
	tspan[1] = T;
	Matrix Q,P,T;

	integ2_symp(func,tspan,q,p,dt,&Q,&P,&T);

	Q = Q.tr();

	Matrix state = Q;

	char path[] = "./data/solution.txt";
	state.save(path);
}

void Model::deriv_symp(Matrix* dy, double t, Matrix* y,int indicator)
{
	if(indicator == 0)
	{
		(*dy) = (*y);
	}
	else
	{
		Matrix F(N,1);
		for(int i=0 ; i<N ; i++)
			F.at(i) = sin( y->at(i) );
		(*dy) = Dxx*(*y) - F + bd;
	}
}

void 

void Model::integ2_symp(void (Model::*func)(Matrix*,double,Matrix*,int),double* tspan, Matrix init_pos, Matrix init_mom, double h, Matrix* q, Matrix* p ,Matrix* T)
{
	double t = tspan[0];
	Matrix da,db,a,b;

	a = init_pos;
	b = init_mom;

	da.zeros( a.get_num_rows() , a.get_num_cols() );
	db.zeros( b.get_num_rows() , b.get_num_cols() );

	T->zeros(MAX_ITER,1);
	T->at(0) = 0;

	for(int j = 0 ; j < MAX_ITER ; j++)
	{
		prog_bar(j,MAX_ITER);

		(this->*func)(&da,t,&b,0);
		da = da.scalar( h/2 );
		a = a + da;

		(this->*func)(&db,t,&a,1);
		db = db.scalar( h/2 );
		b = b + db;

		(this->*func)(&da,t,&b,0);
		da = da.scalar( h/2 );
		a = a + da;

		(this->*func)(&db,t,&a,1);
		db = db.scalar( h/2 );
		b = b + db;

		q->add_row(a.tr());
		p->add_row(b.tr());
		T->at(j) = t;
	}

}

void Model::integ4_symp(void (Model::*func)(Matrix*,double,Matrix*,int),double* tspan, Matrix init_pos, Matrix init_mom, double h, Matrix* q, Matrix* p ,Matrix* T)
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
		prog_bar(j,MAX_ITER);

		(this->*func)(&da,t,&b,0);
		da = da.scalar( c1*h );
		a = a + da;

		(this->*func)(&db,t,&a,1);
		db = db.scalar( d1*h );
		b = b + db;

		(this->*func)(&da,t,&b,0);
		da = da.scalar( c2*h );
		a = a + da;

		(this->*func)(&db,t,&a,1);
		db = db.scalar( d2*h );
		b = b + db;

		(this->*func)(&da,t,&b,0);
		da = da.scalar( c3*h );
		a = a + da;

		(this->*func)(&db,t,&a,1);
		db = db.scalar( d3*h );
		b = b + db;

		(this->*func)(&da,t,&b,0);
		da = da.scalar( c4*h );
		a = a + da;

		(this->*func)(&db,t,&a,1);
		db = db.scalar( d4*h );
		b = b + db;

		q->add_row(a.tr());
		p->add_row(b.tr());
		T->at(j) = t;
	}
}

void Model::DEI()
{
	char path[] = "./data/solution.txt";
	Matrix solution;
	solution.open(N,MAX_ITER,path);

	Matrix F;
	F.zeros(N,MAX_ITER);

	for(int i=0 ; i<20 ; i++)
	{
		for(int j=0 ; j<MAX_ITER ; j++)
		{
			solution.at(i,j) = sin(solution.at(i,j));
		}
	}
	F.append(solution,'r');

	Matrix U,S,V;
	F.svd(&U,&S,&V);


	Matrix Psi;
	Psi = U.get_col(0);

	Matrix P;
	P.zeros(2*N,1);
	int idx;
	double dummy;
	arg_opt(&dummy,&idx,&Psi,'b');
	P.at(idx) = 1;

	for(int i=1 ; i<60 ; i++)
	{
		Matrix v = U.get_col(i);
		Matrix c = (P.tr()*Psi).inv*P.tr()*
	}

//	for(int i=0 ; i<60 ; i++)
//	{
//		
//	}

//	char path1[] = "./data/U.txt";
//	char path2[] = "./data/S.txt";
//	char path3[] = "./data/V.txt";
//
//	U.open(2*N,MAX_ITER,path1);
//	S.open(MAX_ITER,1,path2);

//	for(int i=0 ; i<100 ; i++)
//	{
//		cout << S.at(i) << endl;
//	}
}
