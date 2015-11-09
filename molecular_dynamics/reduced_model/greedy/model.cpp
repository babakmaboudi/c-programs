#include "model.h"

Model::Model()
{
	temperature = 2.;
	density = 0.5;
	cut_off = pow (2., 1./6.);

	initUcell.zeros(3);
	initUcell(0)=10.;
	initUcell(1)=10.;
	initUcell(2)=10.;

	double c = 1. / pow (density, 1./3.);
	region.zeros(3);
	region(0) = initUcell(0)*c;
	region(1) = initUcell(1)*c;
	region(2) = initUcell(2)*c;
	
	initUchain.zeros(3);
	initUchain(0) = 1.;
	initUchain(1) = 1.;
	initUchain(2) = 1.;

	num_chains = 2*initUchain(0)*initUchain(1)*initUchain(2);

	if(num_chains == 2)
		num_chains = 1;

	bond_lim = 2.1;
	chain_length = 16;
	dt = 0.0005;
	num_mols = num_chains*chain_length;

	int NDIM = 3;
	vel_mag = sqrt( NDIM * (1. - 1. / num_mols ) * temperature );
	bond_lim = 2.1;
	vtk_count = 0;

	double bz = cut_off * cos(M_PI / 4.);
	mol_dist = sqrt(2*bz*bz);	
}

void Model::initial_condition(Matrix* init_q , Matrix* init_p)
{
	Matrix t;
	t.linspace(0,2*M_PI-2*M_PI/num_mols,num_mols);
	Matrix X(num_mols);
	Matrix Y(num_mols);
	for(int i=0 ; i<num_mols ; i++)
	{
		X(i) = cos(t(i));
		Y(i) = sin(t(i));
	}

	X *= mol_dist / 2 / sin(M_PI/num_mols);
	Y *= mol_dist / 2 / sin(M_PI/num_mols);

	init_q->zeros(3*num_mols);
	init_p->zeros(3*num_mols);

	for(int m=0 ; m<num_mols ; m++)
	{
		Matrix r(3);
		r(1) = X(m);
		r(2) = Y(m);
		(*init_q).set_submat(3*m,0,r);
	}

	int NDIM = 3;
	vel_mag = sqrt( NDIM * (1. - 1. / num_mols ) * temperature );
	Matrix sum(3);

	for(int i=0 ; i<num_mols ; i++)
	{
		Matrix rnd_dir(3);
		double s, x, y;

		s=2.;
		while(s>1.)
		{
			x = 2. * ((double)rand() / (double)RAND_MAX) - 1.;
			y = 2. * ((double)rand() / (double)RAND_MAX) - 1.;
			s = x*x + y*y;
		}
		rnd_dir(2) = 1. - 2. * s;
		s = 2. * sqrt(1. - s);
		rnd_dir(0) = s*x;
		rnd_dir(1) = s*y;
		
		rnd_dir *= vel_mag;

		init_p->set_submat(3*i,0,rnd_dir);
		sum += rnd_dir;
	}
	for(int i=0 ; i<num_mols ; i++)
	{
		Matrix temp = init_p->get_submat( i*3 , i*3+3 , 0 , 1 );
		temp -= sum*(1./num_mols);
		init_p->set_submat( i*3 , 0 , temp );
	}
}

void Model::leap_frog(Matrix init_q, Matrix init_p, int MAX_ITER)
{
	Matrix dq,dp;

	q = init_q;
	p = init_p;

	dq.zeros( init_q.get_num_rows() , 1 );
	dp.zeros( init_p.get_num_rows() , 1 );
	
	if(snaps.get_num_cols() == 0)
	{
		snaps = q;
		snaps.append(p,'c');
	}
	else
	{
		snaps.append(q,'c');
		snaps.append(p,'c');
	}

	for(int i=0 ; i<MAX_ITER ; i++)
	{
		deriv(&dp,2);
		p += dp*(0.5*dt);

		deriv(&dq,1);
		q += dq*dt;

//		boundary_conditions(&q);

		deriv(&dp,2);
		p += dp*(0.5*dt);

//		export_vtk(&q);

		if(i%50 == 0)
		{
			snaps.append(q,'c');	
			snaps.append(p,'c');
		}
	}
}

void Model::deriv(Matrix* dz, int ind)
{
	if(ind == 1)
	{
		(*dz) = p;
	}
	else
	{
		compute_forces();
		compute_bonded_forces();
		(*dz) = f;
	}
}

void Model::compute_forces()
{
	potential = 0;
	double cc = cut_off*cut_off;

	f.zeros(3*num_mols);

	for(int i=0 ; i<num_mols ; i++)
	{
		for(int j=i+1 ; j<num_mols ; j++)
		{
			if( abs(i-j) > 1 )
			{
				if((i==0) && (j==num_mols-1))
					break;
				int idx1 = 3*i;
				int idx2 = 3*j;
				Matrix mol1 = q.get_submat(idx1,idx1+3,0,1);
				Matrix mol2 = q.get_submat(idx2,idx2+3,0,1);

				Matrix dr = mol1 - mol2;
				wrap(&dr);
				
				double dr2 = dr(0)*dr(0)+dr(1)*dr(1)+dr(2)*dr(2);
				
				if( dr2<cc )
				{
					double dr2i = 1./dr2;
					double dr6i = pow(dr2i,3);
					double force = 48. * dr6i * (dr6i - 0.5) * dr2i;
					potential += 4. * dr6i * (dr6i-1.) + 1.;
					Matrix a1 = f.get_submat(idx1,idx1+3,0,1);
					Matrix a2 = f.get_submat(idx2,idx2+3,0,1);
					
					a1 += dr*(force);
					a2 += dr*(-force);

					f.set_submat(idx1,0,a1);
					f.set_submat(idx2,0,a2);
				}
			}

		}
	}
}

void Model::compute_bonded_forces()
{
	double cc = cut_off*cut_off;

	for(int i=0 ; i<num_chains ; i++)
	{
		for(int j=0 ; j<chain_length-1 ; j++)
		{
			int idx1 = (i*chain_length + j)*3;
			int idx2 = idx1+3;

			Matrix mol1 = q.get_submat(idx1,idx1+3,0,1);
			Matrix mol2 = q.get_submat(idx2,idx2+3,0,1);

			Matrix dr = mol1 - mol2;
			wrap(&dr);

			double dr2 = dr(0)*dr(0)+dr(1)*dr(1)+dr(2)*dr(2);

			if( dr2 < cc )
			{
				double dr2i = 1./dr2;
				double dr6i = pow(dr2i,3);
				double force = 48. * dr6i * (dr6i - 0.5) * dr2i;
				potential += 4. * dr6i * (dr6i-1.) + 1.;

				Matrix a1 = f.get_submat(idx1,idx1+3,0,1);
				Matrix a2 = f.get_submat(idx2,idx2+3,0,1);
				
				a1 += dr*(force);
				a2 += dr*(-force);

				f.set_submat(idx1,0,a1);
				f.set_submat(idx2,0,a2);
			}

			double w = 1. - bond_lim / sqrt(dr2);
			assert(w<0.);
			dr2 *= (w*w);

			if( dr2 < cc )
			{
				double dr2i = 1./dr2;
				double dr6i = pow(dr2i,3);
				double force = 48. * w * dr6i * (dr6i - 0.5) * dr2i;

				potential += 4. * dr6i * (dr6i-1.) + 1.;
				Matrix a1 = f.get_submat(idx1,idx1+3,0,1);
				Matrix a2 = f.get_submat(idx2,idx2+3,0,1);
				
				a1 += dr*(force);
				a2 += dr*(-force);

				f.set_submat(idx1,0,a1);
				f.set_submat(idx2,0,a2);
			}
		}
	}
	int idx1 = 0;
	int idx2 = 3*(num_mols-1);

	Matrix mol1 = q.get_submat(idx1,idx1+3,0,1);
	Matrix mol2 = q.get_submat(idx2,idx2+3,0,1);

	Matrix dr = mol1 - mol2;
	wrap(&dr);

	double dr2 = dr(0)*dr(0)+dr(1)*dr(1)+dr(2)*dr(2);

	if( dr2 < cc )
	{
		double dr2i = 1./dr2;
		double dr6i = pow(dr2i,3);
		double force = 48. * dr6i * (dr6i - 0.5) * dr2i;
		potential += 4. * dr6i * (dr6i-1.) + 1.;

		Matrix a1 = f.get_submat(idx1,idx1+3,0,1);
		Matrix a2 = f.get_submat(idx2,idx2+3,0,1);
		
		a1 += dr*(force);
		a2 += dr*(-force);

		f.set_submat(idx1,0,a1);
		f.set_submat(idx2,0,a2);
	}

	double w = 1. - bond_lim / sqrt(dr2);
	assert(w<0.);
	dr2 *= (w*w);

	if( dr2 < cc )
	{
		double dr2i = 1./dr2;
		double dr6i = pow(dr2i,3);
		double force = 48. * w * dr6i * (dr6i - 0.5) * dr2i;

		potential += 4. * dr6i * (dr6i-1.) + 1.;
		Matrix a1 = f.get_submat(idx1,idx1+3,0,1);
		Matrix a2 = f.get_submat(idx2,idx2+3,0,1);
		
		a1 += dr*(force);
		a2 += dr*(-force);

		f.set_submat(idx1,0,a1);
		f.set_submat(idx2,0,a2);
	}
}

void Model::boundary_conditions(Matrix* mols)
{
	for(int i=0 ; i<num_mols ; i++)
	{
		int idx = i*3;
		Matrix mat = mols->get_submat(idx,idx+3,0,1);
		wrap(&mat);
		mols->set_submat( idx,0,mat );
	}
}

void Model::wrap(Matrix* vec)
{
		if( (*vec)(0) > 0.5*region(0) )
			(*vec)(0) -= region(0);
		else if( (*vec)(0) < -0.5*region(0) )
			(*vec)(0) += region(0);

		if( (*vec)(1) > 0.5*region(1) )
			(*vec)(1) -= region(1);
		else if( (*vec)(1) < -0.5*region(1) )
			(*vec)(1) += region(1);

		if( (*vec)(2) > 0.5*region(2) )
			(*vec)(2) -= region(2);
		else if( (*vec)(2) < -0.5*region(2) )
			(*vec)(2) += region(2);
}

void Model::simulate()
{
	Matrix q0, p0;

	initial_condition(&q0, &p0);
//	boundary_conditions(&q0);

	export_vtk(&q0);
	
	leap_frog(q0,p0,2000);
//	evaluate_parameters();
}

void Model::generate_snapshots()
{
	Matrix init_q, init_p;

	for(int i=0 ; i<100 ; i++)
	{
		initial_condition(&init_q, &init_p);
//		boundary_conditions(&init_q);

		leap_frog(init_q,init_p,500);
	}
	cout << snaps.get_num_rows() << " " << snaps.get_num_cols() << endl;
	char path[] = "./data/snapshots.txt";
	snaps.save(path);
}

void Model::generate_reduced_basis()
{
	char path[] = "./data/snapshots.txt";
	snaps.open(48,2200,path);

	Matrix U_mat,S_mat,V_mat;
	snaps.svd(&U_mat,&S_mat,&V_mat);

	Matrix temp = U_mat.get_submat(0 , U_mat.get_num_rows() , 0 , 32);

	int r,c;
	temp.get_size(&r,&c);
	A.zeros(2*r,2*c);
	A.set_submat(0,0,temp);
	A.set_submat(r,c,temp);
	(A.tr()*A)();
}

void Model::greedy()
{
	vector<Matrix> parameters;

	for(int i=0 ; i<1000 ; i++)
	{
		Matrix q, p;
		initial_condition(&q,&p);
		parameters.push_back(p);
	}
	parameters[0]();
}

void Model::evaluate_parameters()
{
	double vel = 0;
	Matrix v_sum(3);
	for(int i=0 ; i<num_mols ; i++)
	{
		Matrix temp = p.get_submat(i*3,i*3+3,0,1);
		v_sum += temp;
		vel += temp(0)*temp(0) + temp(1)*temp(1) + temp(2)*temp(2);
	}

	double kinetic_E = 0.5*vel/num_mols;
	double potential_E = potential/num_mols;
	cout << right << "        System Summary" << endl;
	cout.width(20); cout << right << "Kinetic Energy: " << kinetic_E << endl;
	cout.width(20); cout << right << "Potential Energy: " << potential_E << endl;
	cout.width(20); cout << right << "Total Energy: " << kinetic_E + potential_E << endl;
	cout << "-*--*--*--*--*--*--*--*-" << endl;
}

void Model::export_vtk(Matrix* mols)
{
	ofstream file;
	char path[50];
	sprintf(path,"./output/out%d.vtk",vtk_count);
	file.open(path);

	file << "# vtk DataFile Version 3.0" << endl;
	file << "vtk output" << endl;
	file << "ASCII" << endl;
	file << "DATASET UNSTRUCTURED_GRID" << endl;
	file << endl;
	file << "POINTS " << num_mols << " float" << endl;

	for(int i=0 ; i<num_mols ; i++)
	{
		file << (*mols)(3*i) << " " << (*mols)(3*i+1) << " " << (*mols)(3*i+2) << endl;
	}

	file << endl;
	file << "CELLS " << num_chains << " " << num_chains*chain_length + num_chains << endl;

	for(int i=0 ; i<num_chains ; i++)
	{
		file << chain_length << " ";
		for(int j=0 ; j<chain_length ; j++)
		{
			file << i*chain_length + j << " ";
		}
		file << endl;
	}

	file << endl;
	file << "CELL_TYPES " << num_chains << endl;
	for(int i=0 ; i<num_chains ; i++)
		file << 4 << endl;
	file.close();
	
	vtk_count++;
}
