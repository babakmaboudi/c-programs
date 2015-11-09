#include "model.h"

Model::Model()
{
	/*
	NDIM = 3;
	num_mols_dir.set(10,10,10);
//	num_mols = num_mols_dir.x() * num_mols_dir.y() * num_mols_dir.z();

	num_sites = 4;
	
	dt = 0.0005;
	temperature = 2.;
	density = 0.5;
	cut_off = 2.38;
	region = num_mols_dir * (1./sqrt(density));
	vel_mag = sqrt( NDIM * (1. - 1. / num_mols ) * temperature );

	vtk_count = 0;

	bond_lim = 2.1;
	chain_length = 8;
	initUchain.set(1,1,1);
	initUcell.set(10,10,10);
	num_chains = initUchain.x()*initUchain.y()*initUchain.z();
	if(num_chains == 2)
	num_chains = 1;
	*/

	NDIM = 3;
	temperature = 2.;
	density = 0.5;
	cut_off = pow (2., 1./6.);
	initUcell.set(10,10,10);
	double c = 1. / pow (density, 1./3.);
	region.set(initUcell.x()*c , initUcell.y()*c , initUcell.z()*c);
	initUchain.set(1,1,1);
	num_chains = 2*initUchain.x()*initUchain.y()*initUchain.z();
	if(num_chains == 2)
		num_chains = 1;
	bond_lim = 2.1;
	chain_length = 16;
	vtk_count = 0;
	dt = 0.0005;
	num_mols = initUcell.x()*initUcell.y()*initUcell.z() + num_chains*chain_length;
	vel_mag = sqrt( NDIM * (1. - 1. / num_mols ) * temperature );
	bond_lim = 2.1;


}

void Model::init_coords()
{
	Vec gap = region/num_mols_dir;
	
	for(int i=0 ; i<num_mols_dir.x() ; i++)
	{
		for(int j=0 ; j<num_mols_dir.y() ; j++)
		{
			for(int k=0 ; k<num_mols_dir.z() ; k++)
			{
				Vec c(i+0.5,j+0.5,k+0.5);
				c *= gap;
				c += region*(-.5);

				Mol m;
				m.r = c;
				mols.push_back(m);
			}
		}
	}
}

void Model::init_coords_chain()
{
	double by = cut_off * cos(M_PI / 4.);
	double bz = cut_off * sin(M_PI / 4.);

	num_mols = 0;
	if(num_chains == 1)
	{
		for(int m=0 ; m<chain_length ; m++)
		{
			Mol mol;
			mol.r.set(0.,(m%2)*by,m*bz);
			mol.r += region*(-0.25);
			mols.push_back(mol);
			num_mols++;
		}
	}
	else
	{
		Vec gap( region.x()/initUchain.x() , region.y()/initUchain.y() , region.z()/initUchain.z() );
		for(int nz=0 ; nz<initUchain.z() ; nz++)
		{
			for(int ny=0 ; ny<initUchain.y() ; ny++)
			{
				for(int nx=0 ; nx<initUchain.x() ; nx++)
				{
					Vec c(0.25+nx , 0.25+ny , 0.25+nz);
					c *= gap;
					c += region*(-0.5);
					
					for(int j=0 ; j<2 ; j++)
					{
						for(int m=0 ; m<chain_length ; m++)
						{
							Mol mol;
							mol.r.set(0,(m%2)*by,m*bz);
							mol.r += gap * (0.5*j);
							mol.r += c;
							mols.push_back(mol);
							num_mols++;
						}
					}
				}
			}
		}
	}
	boundary_conditions();
/*
	Vec gap(region.x()/initUcell.x() , region.y()/initUcell.y() , region.z()/initUcell.z());
	for(int nz=0 ; nz<initUcell.z() ; nz++)
	{
		for(int ny=0 ; ny<initUcell.y() ; ny++)
		{
			for(int nx=0 ; nx<initUcell.x() ; nx++)
			{
				Vec c(nx+0.5 , ny+0.5 , nz+0.5);
				c *= gap;
				c += region*(-0.5);
				
				int i;
				for(i=0 ; i<num_chains*chain_length ; i++)
				{
					Vec dr = mols[i].r - c;
					if(dr.quad3() < cut_off*cut_off)
						break;
				}
				if(i == num_chains*chain_length)
				{
					Mol mol;
					mol.r = c;
					mols.push_back(mol);
				}
			}
		}
	}
*/
	num_mols = mols.size();
}

void Model::assign_to_chain()
{
	for(int i=0 ; i<num_chains ; i++)
	{
		for(int j=0 ; j<chain_length ; j++)
		{
			mols[i*chain_length + j].chain = i;
		}
	}
	for(int i=num_chains*chain_length ; i<num_mols ; i++)
		mols[i].chain = -1;
}

void Model::init_vels()
{
	vel_mag = sqrt( NDIM * (1. - 1. / num_mols ) * temperature );
	Vec sum;
	for(int i=0 ; i<num_mols ; i++)
	{
		mols[i].rv.rand_dir_3D();
		mols[i].rv *= vel_mag;
		sum += mols[i].rv;
	}
	for(int i=0 ; i<num_mols ; i++)
	{
		mols[i].rv -= sum*(1./num_mols);
	}
}

void Model::init_accels()
{
	for(int i=0 ; i<num_mols ; i++)
	{
		mols[i].ra.set(0.0,0.0,0.0);
		mols[i].ra1.set(0.0,0.0,0.0);
		mols[i].ra2.set(0.0,0.0,0.0);
	}
}

void Model::initiate_cells()
{
	num_cells_dir = region/cut_off;
	num_cells = num_cells_dir.x() * num_cells_dir.y() * num_cells_dir.z();
	vector<int> dummy;
	for(int i=0 ; i<num_cells ; i++)
		cell_list.push_back(dummy);

	cell_width.set( region.x()/num_cells_dir.x() , region.y()/num_cells_dir.y() , region.z()/num_cells_dir.z() );
	
	Vec_int n;

	n.set(0,0,0);
	neigs.push_back(n);
	n.set(1,0,0);
	neigs.push_back(n);
	n.set(1,1,0);
	neigs.push_back(n);
	n.set(0,1,0);
	neigs.push_back(n);
	n.set(-1,1,0);
	neigs.push_back(n);
	n.set(0,0,1);
	neigs.push_back(n);
	n.set(1,0,1);
	neigs.push_back(n);
	n.set(1,1,1);
	neigs.push_back(n);
	n.set(0,1,1);
	neigs.push_back(n);
	n.set(-1,1,1);
	neigs.push_back(n);
	n.set(-1,0,1);
	neigs.push_back(n);
	n.set(-1,-1,1);
	neigs.push_back(n);
	n.set(0,-1,1);
	neigs.push_back(n);
	n.set(1,-1,1);
	neigs.push_back(n);
}

void Model::define_mol()
{
	Site s;
	s.r.set(0.0,0.0,0.0);

	s.r.z() = -0.0206;
	s.type = 1;
	ref_sites.push_back(s);
	s.r.set(0.0,0.0,0.0);

	s.r.z() = 0.0274;
	s.type = 2;
	ref_sites.push_back(s);
	s.r.set(0.0,0.0,0.0);
	
	s.r.y() = 0.240;
	s.r.z() = 0.165;
	s.type = 3;
	ref_sites.push_back(s);
	s.r.set(0.0,0.0,0.0);

	s.r.y() = -ref_sites[2].r.y();
	s.r.z() = ref_sites[2].r.z();
	s.type = 3;
	ref_sites.push_back(s);

	mol_inert.set(0.00980,0.00340,0.00640);
	bCon = 183.5;

}

void Model::initialize()
{
	init_coords_chain();
	assign_to_chain();
	init_vels();
	init_accels();
	export_vtk();
/*
	define_mol();

	init_coords();
	init_vels();
	init_accels();
	init_ang_coords();
	init_ang_vels();
	init_ang_accels();

//	initiate_cells();
*/
}

void Model::boundary_conditions()
{
	for(int i=0 ; i<num_mols ; i++)
		wrap(&mols[i].r);
}

void Model::compute_forces()
{
	double cc = cut_off*cut_off;

	for(int i=0 ; i<num_mols ; i++)
		mols[i].ra.set(0.0,0.0,0.0);

	potential = 0;
	virSum = 0;

	for(int i=0 ; i<num_mols-1 ; i++)
	{
		for(int j=i+1 ; j<num_mols ; j++)
		{
			if(mols[i].chain == -1 || mols[i].chain != mols[j].chain || abs(i - j) > 1)
			{
				Vec dr = mols[i].r - mols[j].r;
				wrap(&dr);
				double dr2 = dr.quad3();

				if( dr2 < cc )
				{
					double dr2i = 1./dr2;
					double dr6i = pow(dr2i,3);
					double force = 48. * dr6i * (dr6i - 0.5) * dr2i;

					potential += 4. * dr6i * (dr6i-1.) + 1.;
					virSum += force * dr2;

					mols[i].ra += dr*(force);
					mols[j].ra += dr*(-force);
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
			int idx1 = i*chain_length + j;
			int idx2 = idx1+1;

			Vec dr = mols[idx1].r - mols[idx2].r;
			wrap(&dr);
			double dr2 = dr.quad3();

			if( dr2 < cc )
			{
				double dr2i = 1./dr2;
				double dr6i = pow(dr2i,3);
				double force = 48. * dr6i * (dr6i - 0.5) * dr2i;

				potential += 4. * dr6i * (dr6i-1.) + 1.;
				virSum += force * dr2;

				mols[idx1].ra += dr*(force);
				mols[idx2].ra += dr*(-force);
			}

			double w = 1. - bond_lim / sqrt(dr2);

			assert(w<0.); // check if the bond has snapped

			dr2 *= (w*w);
			if( dr2 < cc )
			{
				double dr2i = 1./dr2;
				double dr6i = pow(dr2i,3);
				double force = 48. * w * dr6i * (dr6i - 0.5) * dr2i;

				potential += 4. * dr6i * (dr6i-1.) + 1.;
				virSum += force * dr2;

				mols[idx1].ra += dr*(force);
				mols[idx2].ra += dr*(-force);
			}

		}
	}
}

void Model::leap_frog(int part)
{
	if(part == 1)
	{
		for(int i=0 ; i<num_mols ; i++)
		{
			mols[i].rv += mols[i].ra*(0.5*dt);
			mols[i].r += mols[i].rv*(dt);
		}
	}
	else
	{
		for(int i=0 ; i<num_mols ; i++)
			mols[i].rv += mols[i].ra*(0.5*dt);
	}
}

void Model::predictor_step()
{
	double cr[] = {19., -10., 3.};
	double cv[] = {27., -22., 7.};
	for(int i=0 ; i<num_mols ; i++)
	{
		mols[i].ro = mols[i].r;
		mols[i].rvo = mols[i].rv;

		pc(&mols[i],cr,cv,0);

		mols[i].ra2 = mols[i].ra1;
		mols[i].ra1 = mols[i].ra;
	}
}

void Model::corrector_step()
{
	double cr[] = {3., 10., -1.};
	double cv[] = {7., 6., -1.};
	for(int i=0 ; i<num_mols ; i++)
		pc(&mols[i],cr,cv,1);
}

void Model::pc(Mol* m, double* cr, double* cv, int ind)
{
	double div = 24.;
	double wr = dt*dt / div;
	double wv = dt / div;

	if(ind == 0)
	{
		m->r.x() = m->r.x() + dt*m->rv.x() + wr * ( cr[0]*m->ra.x() + cr[1]*m->ra1.x() + cr[2]*m->ra2.x() );
		m->rv.x() = ( m->r.x() - m->ro.x() ) / dt + wv * ( cv[0]*m->ra.x() + cv[1]*m->ra1.x() + cv[2]*m->ra2.x() );

		m->r.y() = m->r.y() + dt*m->rv.y() + wr * ( cr[0]*m->ra.y() + cr[1]*m->ra1.y() + cr[2]*m->ra2.y() );
		m->rv.y() = ( m->r.y() - m->ro.y() ) / dt + wv * ( cv[0]*m->ra.y() + cv[1]*m->ra1.y() + cv[2]*m->ra2.y() );

		m->r.z() = m->r.z() + dt*m->rv.z() + wr * ( cr[0]*m->ra.z() + cr[1]*m->ra1.z() + cr[2]*m->ra2.z() );
		m->rv.z() = ( m->r.z() - m->ro.z() ) / dt + wv * ( cv[0]*m->ra.z() + cv[1]*m->ra1.z() + cv[2]*m->ra2.z() );

//		m->q.q1() = m->q.q1() + dt*m->qv.q1() + wr * ( cr[0]*m->qa.q1() + cr[1]*m->qa1.q1() + cr[2]*m->qa2.q1() );
//		m->qv.q1() = ( m->q.q1() - m->qo.q1() ) / dt + wv * ( cv[0]*m->qa.q1() + cv[1]*m->qa1.q1() + cv[2]*m->qa2.q1() );

//		m->q.q2() = m->q.q2() + dt*m->qv.q2() + wr * ( cr[0]*m->qa.q2() + cr[1]*m->qa1.q2() + cr[2]*m->qa2.q2() );
//		m->qv.q2() = ( m->q.q2() - m->qo.q2() ) / dt + wv * ( cv[0]*m->qa.q2() + cv[1]*m->qa1.q2() + cv[2]*m->qa2.q2() );

//		m->q.q3() = m->q.q3() + dt*m->qv.q3() + wr * ( cr[0]*m->qa.q3() + cr[1]*m->qa1.q3() + cr[2]*m->qa2.q3() );
//		m->qv.q3() = ( m->q.q3() - m->qo.q3() ) / dt + wv * ( cv[0]*m->qa.q3() + cv[1]*m->qa1.q3() + cv[2]*m->qa2.q3() );

//		m->q.q4() = m->q.q4() + dt*m->qv.q4() + wr * ( cr[0]*m->qa.q4() + cr[1]*m->qa1.q4() + cr[2]*m->qa2.q4() );
//		m->qv.q4() = ( m->q.q4() - m->qo.q4() ) / dt + wv * ( cv[0]*m->qa.q4() + cv[1]*m->qa1.q4() + cv[2]*m->qa2.q4() );
	}
	else
	{
		m->r.x() = m->ro.x() + dt*m->rvo.x() + wr * ( cr[0]*m->ra.x() + cr[1]*m->ra1.x() + cr[2]*m->ra2.x() );
		m->rv.x() = ( m->r.x() - m->ro.x() ) / dt + wv * ( cv[0]*m->ra.x() + cv[1]*m->ra1.x() + cv[2]*m->ra2.x() );

		m->r.y() = m->ro.y() + dt*m->rvo.y() + wr * ( cr[0]*m->ra.y() + cr[1]*m->ra1.y() + cr[2]*m->ra2.y() );
		m->rv.y() = ( m->r.y() - m->ro.y() ) / dt + wv * ( cv[0]*m->ra.y() + cv[1]*m->ra1.y() + cv[2]*m->ra2.y() );

		m->r.z() = m->ro.z() + dt*m->rvo.z() + wr * ( cr[0]*m->ra.z() + cr[1]*m->ra1.z() + cr[2]*m->ra2.z() );
		m->rv.z() = ( m->r.z() - m->ro.z() ) / dt + wv * ( cv[0]*m->ra.z() + cv[1]*m->ra1.z() + cv[2]*m->ra2.z() );		

//		m->q.q1() = m->qo.q1() + dt*m->qvo.q1() + wr * ( cr[0]*m->qa.q1() + cr[1]*m->qa1.q1() + cr[2]*m->qa2.q1() );
//		m->qv.q1() = ( m->q.q1() - m->qo.q1() ) / dt + wv * ( cv[0]*m->qa.q1() + cv[1]*m->qa1.q1() + cv[2]*m->qa2.q1() );

//		m->q.q2() = m->qo.q2() + dt*m->qvo.q2() + wr * ( cr[0]*m->qa.q2() + cr[1]*m->qa1.q2() + cr[2]*m->qa2.q2() );
//		m->qv.q2() = ( m->q.q2() - m->qo.q2() ) / dt + wv * ( cv[0]*m->qa.q2() + cv[1]*m->qa1.q2() + cv[2]*m->qa2.q2() );

//		m->q.q3() = m->qo.q3() + dt*m->qvo.q3() + wr * ( cr[0]*m->qa.q3() + cr[1]*m->qa1.q3() + cr[2]*m->qa2.q3() );
//		m->qv.q3() = ( m->q.q3() - m->qo.q3() ) / dt + wv * ( cv[0]*m->qa.q3() + cv[1]*m->qa1.q3() + cv[2]*m->qa2.q3() );

//		m->q.q4() = m->qo.q4() + dt*m->qvo.q4() + wr * ( cr[0]*m->qa.q4() + cr[1]*m->qa1.q4() + cr[2]*m->qa2.q4() );
//		m->qv.q4() = ( m->q.q4() - m->qo.q4() ) / dt + wv * ( cv[0]*m->qa.q4() + cv[1]*m->qa1.q4() + cv[2]*m->qa2.q4() );
	}
}

void Model::single_step()
{
/*
	predictor_step();
	compute_forces();
	compute_site_coords();
	compute_site_forces();
	compute_torqs();
	compute_quat_accels();
	corrector_step();
	adjust_quats();
	boundary_conditions();
*/
	
	leap_frog(1);
//	predictor_step();
	boundary_conditions();
	compute_forces();
	compute_bonded_forces();
//	corrector_step();
	leap_frog(2);
//	boundary_conditions();

/*
	leap_frog(1);
	boundary_conditions();
	compute_forces();
//	compute_forces_cell();
	leap_frog(2);
*/
}

void Model::simulate()
{
	vector< vector<Mol> > res;
	res.push_back(mols);
	for(int i=0 ; i<5000 ; i++)
	{
		single_step();
		res.push_back(mols);
		evaluate_parameters();
		export_vtk();
	}
//	export_file(&res);
}

void Model::evaluate_parameters()
{
	Vec vel;
	double vel2 = 0;
	double vel3 = 0;


	for(int i=0 ; i<num_mols ; i++)
	{
		vel += mols[i].rv;
		vel2 += mols[i].rv.quad3();

		Vec w;
//		compute_ang_vel(i,&w);
//		vel3 += mol_inert.x()*w.x()*w.x() + mol_inert.y()*w.y()*w.y() + mol_inert.z()*w.z()*w.z();
	}
	
	double kinetic_energy = 0.5 * vel2 / num_mols;
	double rotational_kinetic_energy = 0;//0.5 * vel3 / num_mols;
	double potential_energy = potential / num_mols;
	double total_energy = kinetic_energy + potential_energy + rotational_kinetic_energy;
	double pressure = density * (vel2 + virSum) / (num_mols * NDIM);
	
	cout << "summary of system:" << endl;
	cout.width(20); cout << right<< "kinetic energy : " << kinetic_energy << endl;
//	cout.width(20); cout << right<< " rot kinetic energy : " << rotational_kinetic_energy << endl;
	cout.width(20); cout << right << "potential energy : " << potential_energy << endl;
	cout.width(20); cout << right << "total energy : " << total_energy << endl;
	cout.width(20); cout << right << "pressure : " << pressure << endl;
//	vel.print();
	cout << "-*--*--*--*--*--*--*--*-" << endl;

}

void Model::wrap(Vec* v)
{
	if(NDIM == 2)
	{
		if(v->x() > 0.5*region.x())
			v->x() -= region.x();
		else if(v->x() < -0.5*region.x())
			v->x() += region.x();

		if(v->y() > 0.5*region.y())
			v->y() -= region.y();
		else if(v->y() < -0.5*region.y())
			v->y() += region.y();
	}
	else
	{
		if(v->x() > 0.5*region.x())
			v->x() -= region.x();
		else if(v->x() < -0.5*region.x())
			v->x() += region.x();

		if(v->y() > 0.5*region.y())
			v->y() -= region.y();
		else if(v->y() < -0.5*region.y())
			v->y() += region.y();

		if(v->z() > 0.5*region.z())
			v->z() -= region.z();
		else if(v->z() < -0.5*region.z())
			v->z() += region.z();
	}
}

void Model::wrap_cell(Vec_int* v, Vec* shift)
{
	if( v->x() >= num_cells_dir.x() )
	{
		v->x() = 0;
		shift->x() = region.x();
	}
	else if( v->x() < 0 )
	{
		v->x() = num_cells_dir.x()-1;
		shift->x() = -region.x();
	}

	if( v->y() >= num_cells_dir.y() )
	{
		v->y() = 0;
		shift->y() = region.y();
	}
	else if( v->y() < 0 )
	{
		v->y() = num_cells_dir.y()-1;
		shift->y() = -region.y();
	}

	if( v->z() >= num_cells_dir.z() )
	{
		v->z() = 0;
		shift->z() = region.z();
	}
	else if( v->z() < 0 )
	{
		v->z() = num_cells_dir.z()-1;
		shift->z() = -region.z();
	}
}

void Model::compute_shift(Vec* shift, Vec* v)
{
	shift->set(0.0,0.0,0.0);
	if(v->x() > 0.5*region.x())
		shift->x() -= region.x();
	else if(v->x() < -0.5*region.x())
		shift->x() += region.x();

	if(v->y() > 0.5*region.y())
		shift->y() -= region.y();
	else if(v->y() < -0.5*region.y())
		shift->y() += region.y();

	if(v->z() > 0.5*region.z())
		shift->z() -= region.z();
	else if(v->z() < -0.5*region.z())
		shift->z() += region.z();

}

void Model::export_file(vector< vector<Mol> >* res)
{
	ofstream file;
	char path[] = "./data.txt";
	file.open(path);

	for(int i=0 ; i<num_mols ; i++)
	{
		for(int j=0 ; j<res->size() ; j++)
		{
			file << ((*res)[j])[i].r.x() << " " << ((*res)[j])[i].r.y() << " " << ((*res)[j])[i].r.z() << " ";
		}
		file << endl;
	}
	file.close();
}

void Model::export_vtk()
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
		file << mols[i].r.x() << " " << mols[i].r.y() << " " << mols[i].r.z() << endl;
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

void Model::print()
{
	cout << "-----------------------------" << endl;
	for(int i=0 ; i<num_mols ; i++)
		cout << mols[i].r.x() << " " << mols[i].r.y() << " " << mols[i].r.z() << endl;
	cout << "-----------------------------" << endl;
}
