#include "model.h"

Model::Model()
{
	NDIM = 3;
	num_mols_dir.set(6,6,6);
	num_mols = num_mols_dir.x() * num_mols_dir.y() * num_mols_dir.z();
	
	dt = 0.0005;
	temperature = 3.8;
	density = 0.98;
	cut_off = 2.38;
	region = num_mols_dir * (1./sqrt(density));
	vel_mag = sqrt( NDIM * (1. - 1. / num_mols ) * temperature );
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

void Model::init_vels()
{
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

void Model::initialize()
{
	init_coords();
	init_vels();
	initiate_cells();
}

void Model::boundary_conditions()
{
	for(int i=0 ; i<num_mols ; i++)
		wrap(&mols[i].r);
}

void Model::compute_forces()
{
	real cc = cut_off*cut_off;

	for(int i=0 ; i<num_mols ; i++)
		mols[i].ra.set(0.0,0.0,0.0);

	potential = 0;
	virSum = 0;

	for(int i=0 ; i<num_mols-1 ; i++)
	{
		for(int j=i+1 ; j<num_mols ; j++)
		{
			Vec dr = mols[i].r - mols[j].r;
			wrap(&dr);
			real dr2 = dr.quad3();

			if( dr2 < cc )
			{
				real dr2i = 1./dr2;
				real dr6i = pow(dr2i,3);
				real force = 48. * dr6i * (dr6i - 0.5) * dr2i;

				potential += 4. * dr6i * (dr6i-1.) + 1.;
				virSum += force * dr2;

				mols[i].ra += dr*(force);
				mols[j].ra += dr*(-force);
			}
		}
		mols[i].ra.print();
		getchar();
	}
}

void Model::compute_forces_cell()
{
	cell_list.clear();
	vector<int> dummy;
	for(int i=0 ; i<num_cells ; i++)
		cell_list.push_back(dummy);
	for(int i=0 ; i<num_mols ; i++)
	{
		Vec rs = mols[i].r + region*(0.5);
		Vec_int cc(rs.x()/cell_width.x(),rs.y()/cell_width.y(),rs.z()/cell_width.z());
		cell_list[cc.get_cell_num(num_cells_dir)].push_back(i);
	}

	for(int i=0 ; i<num_mols ; i++)
		mols[i].ra.set(0.0,0.0,0.0);

	real cc = cut_off*cut_off;

	potential = 0;
	virSum = 0;

	for(int m1z=0 ; m1z<num_cells_dir.z() ; m1z++)
	{
		for(int m1y=0 ; m1y<num_cells_dir.y() ; m1y++)
		{
			for(int m1x=0 ; m1x<num_cells_dir.x() ; m1x++)
			{
				Vec_int m1v(m1x,m1y,m1z);
				int cidx1 = m1v.get_cell_num(num_cells_dir);
				for(int offset=0 ; offset<neigs.size() ; offset++)
				{
					Vec_int m2v = m1v + neigs[offset];
					Vec shift;
					wrap_cell(&m2v,&shift);
					int cidx2 = m2v.get_cell_num(num_cells_dir);
					for(int i=0 ; i<cell_list[cidx1].size() ; i++)
					{
						for(int j=0 ; j<cell_list[cidx2].size() ; j++)
						{
							int midx1 = cell_list[cidx1][i];
							int midx2 = cell_list[cidx2][j];
							if(cidx1!=cidx2 || midx2<midx1)
							{
								Vec dr = mols[ midx1 ].r - mols[ midx2 ].r;
								dr -= shift;
								real dr2 = dr.quad3();

								if( dr2<cc )
								{
									real dr2i = 1./dr2;
									real dr6i = pow(dr2i,3);
									real force = 48. * dr6i * (dr6i - 0.5) * dr2i;

									potential += 4. * dr6i * (dr6i-1.) + 1.;
									virSum += force * dr2;

									mols[midx1].ra += dr*(force);
									mols[midx2].ra += dr*(-force);
								}
							}
						}
					}
				}
			}
		}
	}
}

void Model::compute_quat_accels()
{
	for(int i=0 ; i<num_mols ; i++)
	{
		Vec w;
		compute_ang_vel(i,&w);

		Quat qs;
		qs.q1() = mols[i].torq.x() + ( body_inert.y() - body_inert.z() ) * w.y() * w.z() / body_inert.x();
	       	qs.q2() = mols[i].torq.y() + ( body_inert.z() - body_inert.x() ) * w.z() * w.x() / body_inert.y();
		qs.q3() = mols[i].torq.z() + ( body_inert.x() - body_inert.y() ) * w.x() * w.y() / body_inert.z();
		qs.q4() = -2 * mols[i].qv.quad();
		mols[i].qa = ( mols[i].q * qs ) * 0.5;
	}
}

void Model::compute_ang_vel(int n, Vec* w)
{
	Quat qt, qvt;

	qvt = mols[n].qv;
	qvt.q4() *= -1;
	qt = (qvt * mols[n].q) * 2;
	w->set( qt.q1(), qt.q2(), qt.q3());
}

void Model::compute_torqs()
{
	Vec torq;
	for(int i=0 ; i< num_mols ; i++)
	{
		mols[i].ra.set(0.0,0.0,0.0);
		torq.set(0.0,0.0,0.0);
		for(int j=0 ; j<all_sites[i].size() ; j++)
		{
			Site site = all_sites[i][j];
			mols[i].ra += site.f;
			Vec dr = site.r - mols[i].r;
			Vec t = dr.cross(site.f);
			torq += t;
		}
		Matrix rot;
		compute_rotation_mat(&rot,&mols[i].q,0);
		mols[i].torq = rot * torq;
	}
}

void Model::compute_rotation_mat(Matrix* rot, Quat* q, int transpose)
{
	real p[10], tq[4], s;
	
	tq[0] = q->q1();
	tq[1] = q->q2();
	tq[2] = q->q3();
	tq[3] = q->q4();
	for(int k=0, k2=0 ; k2<4 ; k2++)
	{
		for(int k1=k2 ; k1<4 ; k1++, k2++)
			p[k] = 2. * tq[k1] * tq[k2];
	}
	rot->at(0,0) = p[0] + p[9] - 1.;
	rot->at(1,1) = p[4] + p[9] - 1.;
	rot->at(2,2) = p[7] + p[9] - 1.;
	s = transpose ? 1. : -1.;
	rot->at(0,1) = p[1] + s * p[8];
	rot->at(1,0) = p[1] - s * p[8];
	rot->at(0,2) = p[2] - s * p[6];
	rot->at(2,0) = p[2] + s * p[6];
	rot->at(1,2) = p[5] + s * p[3];
	rot->at(1,2) = p[5] - s * p[3];
}

void Model::compute_site_coords()
{
	all_sites.clear();
	vector<Site> dummy;

	for(int i=0 ; i<num_mols ; i++)
	{
		Matrix rot;
		compute_rotation_mat( &rot , &mols[i].q , 1);
		
		all_sites.push_back(dummy);

		for(int j=0 ; j<sites_mol ; j++)
		{
			Vec t = rot * ref_sites[j].r;

			Site s;
			s.r = mols[i].r + t;

			all_sites[i].push_back( s );
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
	real cr[] = {19., -10., 3.};
	real cv[] = {27., -22., 7.};
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
	real cr[] = {3., 10., -1.};
	real cv[] = {7., 6., -1.};
	for(int i=0 ; i<num_mols ; i++)
		pc(&mols[i],cr,cv,1);
}

void Model::pc(Mol* m, real* cr, real* cv, int ind)
{
	real div = 24.;
	real wr = dt*dt / div;
	real wv = dt / div;

	if(ind == 0)
	{
		m->r.x() = m->r.x() + dt*m->rv.x() + wr * ( cr[0]*m->ra.x() + cr[1]*m->ra1.x() + cr[2]*m->ra2.x() );
		m->rv.x() = ( m->r.x() - m->ro.x() ) / dt + wv * ( cv[0]*m->ra.x() + cv[1]*m->ra1.x() + cv[2]*m->ra2.x() );

		m->r.y() = m->r.y() + dt*m->rv.y() + wr * ( cr[0]*m->ra.y() + cr[1]*m->ra1.y() + cr[2]*m->ra2.y() );
		m->rv.y() = ( m->r.y() - m->ro.y() ) / dt + wv * ( cv[0]*m->ra.y() + cv[1]*m->ra1.y() + cv[2]*m->ra2.y() );

		m->r.z() = m->r.z() + dt*m->rv.z() + wr * ( cr[0]*m->ra.z() + cr[1]*m->ra1.z() + cr[2]*m->ra2.z() );
		m->rv.z() = ( m->r.z() - m->ro.z() ) / dt + wv * ( cv[0]*m->ra.z() + cv[1]*m->ra1.z() + cv[2]*m->ra2.z() );
	}
	else
	{
		m->r.x() = m->ro.x() + dt*m->rvo.x() + wr * ( cr[0]*m->ra.x() + cr[1]*m->ra1.x() + cr[2]*m->ra2.x() );
		m->rv.x() = ( m->r.x() - m->ro.x() ) / dt + wv * ( cv[0]*m->ra.x() + cv[1]*m->ra1.x() + cv[2]*m->ra2.x() );

		m->r.y() = m->ro.y() + dt*m->rvo.y() + wr * ( cr[0]*m->ra.y() + cr[1]*m->ra1.y() + cr[2]*m->ra2.y() );
		m->rv.y() = ( m->r.y() - m->ro.y() ) / dt + wv * ( cv[0]*m->ra.y() + cv[1]*m->ra1.y() + cv[2]*m->ra2.y() );

		m->r.z() = m->ro.z() + dt*m->rvo.z() + wr * ( cr[0]*m->ra.z() + cr[1]*m->ra1.z() + cr[2]*m->ra2.z() );
		m->rv.z() = ( m->r.z() - m->ro.z() ) / dt + wv * ( cv[0]*m->ra.z() + cv[1]*m->ra1.z() + cv[2]*m->ra2.z() );		
	}
}

void Model::single_step()
{
/*
	predictor_step();
	boundary_conditions();
	compute_forces_cell();
	corrector_step();
	boundary_conditions();
*/

	leap_frog(1);
	boundary_conditions();
	compute_forces();
	compute_forces_cell();
	leap_frog(2);

}

void Model::simulate()
{
	vector< vector<Mol> > res;
	res.push_back(mols);
	for(int i=0 ; i<500 ; i++)
	{
		single_step();
		res.push_back(mols);
		evaluate_parameters();
//		correct_temperature();
	}
	export_file(&res);
}

void Model::correct_temperature()
{
	real sum = 0;
	for(int i=0 ; i<num_mols ; i++)
		sum += mols[i].rv.quad3();

	real factor = vel_mag / sqrt( sum / num_mols );
	for(int i=0 ; i<num_mols ; i++)
		mols[i].rv *= factor;
}

void Model::evaluate_parameters()
{
	Vec vel;
	real vel2 = 0;

	for(int i=0 ; i<num_mols ; i++)
	{
		vel += mols[i].rv;
		vel2 += mols[i].rv.quad3();
	}
	
	real kinetic_energy = 0.5 * vel2 / num_mols;
	real potential_energy = potential / num_mols;
	real total_energy = kinetic_energy + potential_energy;
	real pressure = density * (vel2 + virSum) / (num_mols * NDIM);
	
	cout << "summary of system:" << endl;
	cout.width(20); cout << right<< "kinetic energy : " << kinetic_energy << endl;
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

void Model::print()
{
	cout << "-----------------------------" << endl;
	for(int i=0 ; i<num_mols ; i++)
		cout << mols[i].r.x() << " " << mols[i].r.y() << " " << mols[i].r.z() << endl;
	cout << "-----------------------------" << endl;
}
