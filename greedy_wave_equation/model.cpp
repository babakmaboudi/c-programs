#include "model.h"

Model::Model()
{
	N = 500;
	L = 1;
	dx = L/N;

	T=50;
	dt = 0.01;

	c = 0.1;

	TOL = 1e-4;
}

void Model::initial_condition(Matrix* q_init, Matrix* p_init)
{
	Matrix X;
	X.linspace(0,L,N);
	q_init->zeros(N);
	p_init->zeros(N);

	for(int i=0 ; i<N ; i++)
	{
		double s = 10 * abs( X(i) - 1./2. );
		if( s<=1 )
			(*q_init)(i) = 1.0 - 3.0/2.0*s*s + 3.0/4.0*s*s*s;
		else if( s<=2 )
			(*q_init)(i) = 1.0/4.0*(2.0-s)*(2.0-s)*(2.0-s);
		else
			(*q_init)(i) = 0;
	}
}

void Model::build_stiffness()
{
	Dxx.zeros(N,N);
	Dxx(0,0) = -2.;
	Dxx(0,1) = 1.;
	for(int i=1 ; i<N+1 ; i++)
	{
		Dxx(i,i-1) = 1.;
		Dxx(i,i) = -2.;
		Dxx(i,i+1) = 1.;
	}
	Dxx(N-1,N-2) = 1.;
	Dxx(N-1,N-1) = -2.;
	
	Dxx(0,N-1) = 1;
	Dxx(N-1,0) = 1;
	
	Dxx *= (c*c/dx/dx);
}

void Model::leap_frog(void (Model::*func)(Matrix*,int), Matrix* res_q, Matrix* res_p,Matrix q0, Matrix p0, int MAX_ITER)
{
	Matrix dq,dp;

	q = q0;
	p = p0;

	dq.zeros( q0.get_num_rows() , 1 );
	dp.zeros( p0.get_num_rows() , 1 );

	Matrix snap_q;
	Matrix snap_p;

	snap_q = q;
	snap_p = p;

	(*res_q) = q;
	(*res_p) = p;

	for(int i=0 ; i<MAX_ITER ; i++)
	{
		(this->*func)(&dp,2);
		p += dp*(0.5*dt);

		(this->*func)(&dq,1);
		q += dq*dt;

		(this->*func)(&dp,2);
		p += dp*(0.5*dt);

		snap_q.append(q,'c');
		snap_p.append(p,'c');

		res_q->append(q,'c');
		res_p->append(p,'c');
	}
//	char path1[] = "./data/snap_q_old.txt";
//	snap_q.save(path1);

//	char path2[] = "./data/snap_p_old.txt";
//	snap_p.save(path2);
}

void Model::deriv(Matrix* dz, int ind)
{
	if(ind == 1)
	{
		(*dz) = p;
	}
	else
	{
		(*dz) = Dxx*q;
	}
}

void Model::deriv_reduced(Matrix* dz,int ind)
{
	if(ind == 1)
	{
		(*dz) = Aq*p;
	}
	else
	{
		(*dz) = Ap*q;
	}
}

void Model::build_reduced_basis()
{
	Matrix snap_q, snap_p;
	char path1[] = "./data/snap_q.txt";
	snap_q.open(N,100,path1);

	char path2[] = "./data/snap_p.txt";
	snap_p.open(N,100,path2);

//	double temp = (double) std::rand() / RAND_MAX;
//	temp = (double) std::rand() / RAND_MAX;
//	int index = static_cast<int>(temp*N);
	int index = 65;

	Matrix v1;
	Matrix v2;

	v1 = snap_q.get_submat(0,N,index,index+1);
	v2 = snap_p.get_submat(0,N,index,index+1);

	Matrix E(2*N);
	Matrix F(2*N);

	E.set_submat(0,0,v1);
	F.set_submat(N,0,v2);

	J2n.jay(E.length()/2);

	normalize(&E,&F);
	A = E;
	A.append(F,'c');

	int checker = 0;
	int num_snaps = snap_q.get_num_cols();
	int counter = 1;

	while(checker == 0)
	{
		J2k.jay(counter);
		A_inv = J2k.tr()*A.tr()*J2n;
		
		Matrix samp1(num_snaps);
		Matrix samp2(num_snaps);
		for(int i=0 ; i<num_snaps ; i++)
		{
			Matrix sq = snap_q.get_submat(0,N,i,i+1);
			Matrix sp = snap_p.get_submat(0,N,i,i+1);
			samp1(i) = error_function(&sq,&sp,0);
			samp2(i) = error_function(&sq,&sp,1);
		}
		
		double dummy1 = samp1.max(&index);
		v1 = snap_q.get_submat(0,N,index,index+1);
		Matrix e(2*N);
		e.set_submat(0,0,v1);
		E.append(e,'c');
		

		double dummy2 = samp2.max(&index);
		v2 = snap_p.get_submat(0,N,index,index+1);
		Matrix f(2*N);
		f.set_submat(N,0,v2);
		F.append(f,'c');

		symplectic_qr(&E,&F);
		symplectic_qr(&E,&F);

		A = E;
		A.append(F,'c');
		counter++;

		if( (dummy1 < TOL) && (dummy2 < TOL) )
			checker = 1;

		cout << dummy1 << " " << dummy2 << endl;
	}
	cout << A.get_num_cols() << endl;
	char path[] = "./data/symplectic_basis.txt";
	
	A.save(path);
}

void Model::symplectic_qr(Matrix* mat1, Matrix* mat2)
{
	int size = mat1->get_num_cols();
	Matrix E = mat1->get_submat(0,2*N,0,1);
	Matrix F = mat2->get_submat(0,2*N,0,1);

	normalize(&E,&F);

	for(int i=1 ; i<size ; i++)
	{
		Matrix e = mat1->get_submat(0,2*N,i,i+1);
		Matrix f = mat2->get_submat(0,2*N,i,i+1);

		for(int j=0 ; j<E.get_num_cols() ; j++)
		{
			Matrix vq = E.get_submat(0,2*N,j,j+1);
			Matrix vp = F.get_submat(0,2*N,j,j+1);

			symplectify(&e,&vq,&vp);
			symplectify(&f,&vq,&vp);
		}
		normalize(&e,&f);
		E.append(e,'c');
		F.append(f,'c');
	}
	(*mat1) = E;
	(*mat2) = F;
}

void Model::symplectify(Matrix* vec, Matrix* e, Matrix* f)
{
	double alpha = -vec->bilinear(f);
	double beta = vec->bilinear(e);
	(*vec) = (*vec) + (*e)*(alpha) + (*f)*(beta);
}

void Model::normalize(Matrix* e, Matrix* f)
{
	double alpha = e->bilinear(f);
	(*e) /= sqrt( abs(alpha) );
	(*f) /= sqrt( abs(alpha) );
	if(alpha < 0)
	(*e) *= -1;
}

double Model::error_function(Matrix* sq, Matrix* sp, int ind)
{
	if(ind == 0)
	{
		Matrix snap(2*N);
		snap.set_submat(0,0,(*sq));
		Matrix vec = A*A_inv*snap;
		Matrix vq = vec.get_submat(0,N,0,1);
		Matrix vp = vec.get_submat(N,2*N,0,1);
		return abs(hamiltonian(&vq,&vp,ind)-hamiltonian(sq,sp,ind));
	}
	else
	{
		Matrix snap(2*N);
		snap.set_submat(N,0,(*sp));
		Matrix vec = A*A_inv*snap;
		Matrix vq = vec.get_submat(0,N,0,1);
		Matrix vp = vec.get_submat(N,2*N,0,1);
		return abs(hamiltonian(&vq,&vp,ind)-hamiltonian(sq,sp,ind));		
	}
}

double Model::hamiltonian(Matrix* q, Matrix* p, int ind)
{
	double h = 0;
	if(ind == 0)
	{
		for(int i=0 ; i<N-1 ; i++)
			h += c*c*(q->at(i+1) - q->at(i))*(q->at(i+1) - q->at(i))/4/dx/dx;
		h += c*c*(q->at(N-1) - q->at(0))*(q->at(N-1) - q->at(0))/4/dx/dx;
		for(int i=1 ; i<N ; i++)
			h += c*c*(q->at(i) - q->at(i-1))*(q->at(i) - q->at(i-1))/4/dx/dx;
		h += c*c*(q->at(0) - q->at(N-1))*(q->at(0) - q->at(N-1))/4/dx/dx;
	}
	else if(ind == 1)
	{
		for(int i=0 ; i<N ; i++)
			h += p->at(i) * p->at(i) / 2;
	}
	else
	{
		for(int i=0 ; i<N ; i++)
			h += p->at(i) * p->at(i) / 2;
		for(int i=0 ; i<N-1 ; i++)
			h += c*c*(q->at(i+1) - q->at(i))*(q->at(i+1) - q->at(i))/4/dx/dx;
		h += c*c*(q->at(N-1) - q->at(0))*(q->at(N-1) - q->at(0))/4/dx/dx;
		for(int i=1 ; i<N ; i++)
			h += c*c*(q->at(i) - q->at(i-1))*(q->at(i) - q->at(i-1))/4/dx/dx;
		h += c*c*(q->at(0) - q->at(N-1))*(q->at(0) - q->at(N-1))/4/dx/dx;
	}
	return dx*h;
}

void Model::simulate()
{
	Matrix q0, p0;

	initial_condition(&q0,&p0);
	build_stiffness();

	int MAX_ITER = T / dt;

	Matrix res_q,res_p;
	void (Model::*func)(Matrix*,int) = &Model::deriv;

	clock_t begin = clock();
	leap_frog(func,&res_q,&res_p,q0,p0,1000);
	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout << "Elapsed time for the original system : " << elapsed_secs << "s " << endl;
}

void Model::simulate_reduced()
{
	Matrix q0, p0;

	initial_condition(&q0,&p0);
	build_stiffness();
	
	Matrix L(2*N,2*N);
	Matrix I;
	I.eye(N);

	L.set_submat(0,N,I);
	L.set_submat(N,0,Dxx);

	K = 18;
	J2n.jay(N);
	J2k.jay(K/2);

	char path[] = "./data/symplectic_basis.txt";
	A.open(1000,18,path);

	A_inv = J2k.tr()*A.tr()*J2n;
	Ar = A_inv*L*A;
	Aq = Ar.get_submat(0,K/2,K/2,K);
	Ap = Ar.get_submat(K/2,K,0,K/2);

	Matrix state(2*N);
	state.set_submat(0,0,q0);
	state.set_submat(N,0,p0);

	Matrix y0 = A_inv * state;
	Matrix yq0 = y0.get_submat(0,K/2,0,1);
	Matrix yp0 = y0.get_submat(K/2,K,0,1);

	void (Model::*func)(Matrix*,int) = &Model::deriv_reduced;

	Matrix res_red_q, res_red_p;

	clock_t begin = clock();
	leap_frog(func,&res_red_q,&res_red_p,yq0,yp0,1000);
	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout << "Elapsed time for the reduced system : " << elapsed_secs << "s " << endl;
	
	Matrix y = res_red_q;
	y.append(res_red_p,'r');
	Matrix res = A*y;
	
	char path2[] = "./data/solution_reduced.txt";
	res.save(path2);
}
