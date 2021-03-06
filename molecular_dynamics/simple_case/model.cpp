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
	real c = 1. / pow (density, 1./3.);
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
	real by = cut_off * cos(M_PI / 4.);
	real bz = cut_off * sin(M_PI / 4.);

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

void Model::init_ang_coords()
{
	for(int i=0 ; i<num_mols ; i++)
	{
		Vec e, e_ang;
		e.rand_dir_3D();
		e_ang.x() = atan2(e.x(),e.y());
		e_ang.y() = acos(e.z());
		e_ang.z() = 2. * M_PI * (real)rand() / (real)RAND_MAX;
		euler_to_quat(&mols[i].q, &e_ang);
	}
}

void Model::euler_to_quat(Quat* q, Vec* v)
{
	real a1 = 0.5 * v->y();
	real a2 = 0.5 * ( v->x() - v->z() );
	real a3 = 0.5 * ( v->x() + v->z() );
	q->set( sin(a1)*cos(a2), sin(a1)*sin(a2), cos(a1)*sin(a3), cos(a1)*cos(a3) );
}

void Model::init_ang_vels()
{
	for(int i=0 ; i<num_mols ; i++)
	{
		Vec e;
		e.rand_dir_3D();
		Quat q( e.x(), e.y(), e.z(), 0.0 );
		mols[i].qv = mols[i].q*q;
		real f = 0.5 * vel_mag / sqrt( mol_inert.x()*e.x()*e.x() + mol_inert.y()*e.y()*e.y() + mol_inert.z()*e.z()*e.z() );
		mols[i].qv *= f;
	}
}

void Model::init_ang_accels()
{
	for(int i=0 ; i<num_mols ; i++)
	{
		mols[i].qa.set(0.0,0.0,0.0,0.0);
		mols[i].qa1.set(0.0,0.0,0.0,0.0);
		mols[i].qa2.set(0.0,0.0,0.0,0.0);
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
//	export_vtk();

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
	real cc = cut_off*cut_off;

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
		}
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

void Model::compute_site_forces()
{
	real cc = cut_off*cut_off;

	for(int i=0 ; i<num_mols ; i++)
		for(int j=0 ; j<num_sites ; j++)
			all_sites[i][j].f.set(0.0,0.0,0.0);

	potential = 0;
	for(int i=0 ; i<num_mols-1 ; i++)
	{
		for(int j=i+1 ; j<num_mols ; j++)
		{
			Vec dr = mols[i].r - mols[j].r;
			Vec shift;
			compute_shift(&shift,&dr);
			dr +=  shift;

			real dr2 = dr.quad3();

			if(dr2 < cc)
			{
				for(int si=0 ; si<num_sites ; si++)
				{
					for(int sj=0 ; sj<num_sites ;sj++)
					{
						int type_sum = ref_sites[si].type + ref_sites[sj].type;
						if(ref_sites[si].type == ref_sites[sj].type || type_sum == 5)
						{
							Vec drs = all_sites[i][si].r - all_sites[j][sj].r;
							drs += shift;

							real drs2 = drs.quad3();
							real drs2i = 1./drs2;
							real force=0.0, pot=0.0, drs6i;
							switch(type_sum)
							{
								case 2:
									drs6i = pow(drs2i,3);
									pot = 4. * drs6i * (drs6i - 1.);
									force = 48. * drs6i * (drs6i - 0.5) * drs2i;
									break;
								case 4:
									pot = 4. * bCon * sqrt( drs2i );
									force = pot * drs2i;
									break;
								case 5:
									pot = -2 * bCon * sqrt( drs2i );
									force = pot * drs2i;
									break;
								case 6:
									pot = bCon * sqrt( drs2i );
									force = pot * drs2i;
									break;
							}
							all_sites[i][si].f += drs * force;
							all_sites[j][sj].f += drs * (-force);
							potential += pot;
						}
					}
				}
			}
		}
	}
}

void Model::compute_bonded_forces()
{
	real cc = cut_off*cut_off;

	for(int i=0 ; i<num_chains ; i++)
	{
		for(int j=0 ; j<chain_length-1 ; j++)
		{
			int idx1 = i*chain_length + j;
			int idx2 = idx1+1;

			Vec dr = mols[idx1].r - mols[idx2].r;
			wrap(&dr);
			real dr2 = dr.quad3();

			if( dr2 < cc )
			{
				real dr2i = 1./dr2;
				real dr6i = pow(dr2i,3);
				real force = 48. * dr6i * (dr6i - 0.5) * dr2i;

				potential += 4. * dr6i * (dr6i-1.) + 1.;
				virSum += force * dr2;

				mols[idx1].ra += dr*(force);
				mols[idx2].ra += dr*(-force);
			}

			real w = 1. - bond_lim / sqrt(dr2);

			assert(w<0.); // check if the bond has snapped

			dr2 *= (w*w);
			if( dr2 < cc )
			{
				real dr2i = 1./dr2;
				real dr6i = pow(dr2i,3);
				real force = 48. * w * dr6i * (dr6i - 0.5) * dr2i;

				potential += 4. * dr6i * (dr6i-1.) + 1.;
				virSum += force * dr2;

				mols[idx1].ra += dr*(force);
				mols[idx2].ra += dr*(-force);
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
		qs.q1() = mols[i].torq.x() + ( mol_inert.y() - mol_inert.z() ) * w.y() * w.z() / mol_inert.x();
	       	qs.q2() = mols[i].torq.y() + ( mol_inert.z() - mol_inert.x() ) * w.z() * w.x() / mol_inert.y();
		qs.q3() = mols[i].torq.z() + ( mol_inert.x() - mol_inert.y() ) * w.x() * w.y() / mol_inert.z();
		qs.q4() = -2 * mols[i].qv.quad();
		mols[i].qa = ( mols[i].q * qs ) * 0.5;
	}
}

void Model::compute_ang_vel(int n, Vec* w)
{
	Quat qt, qvt;

	qvt = mols[n].qv;
	qvt.q4() *= -1;
	qt = (qvt * mols[n].q) * 2.;
	w->set( qt.q1(), qt.q2(), qt.q3());
}

void Model::compute_torqs()
{
	Vec torq;
	for(int i=0 ; i< num_mols ; i++)
	{
		mols[i].ra.set(0.0,0.0,0.0);
		torq.set(0.0,0.0,0.0);
		for(int j=0 ; j<num_sites ; j++)
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
/*
	real p[10], tq[4], s;
	
	tq[0] = q->q1();
	tq[1] = q->q2();
	tq[2] = q->q3();
	tq[3] = q->q4();
	for(int k=0, k2=0 ; k2<4 ; k2++)
	{
		for(int k1=k2 ; k1<4 ; k1++, k++)
			p[k] = 2. * tq[k1] * tq[k2];
	}
	rot->at(0) = p[0] + p[9] - 1.;
	rot->at(4) = p[4] + p[9] - 1.;
	rot->at(8) = p[7] + p[9] - 1.;
	s = transpose ? 1. : -1.;
	rot->at(1) = p[1] + s * p[8];
	rot->at(3) = p[1] - s * p[8];
	rot->at(2) = p[2] - s * p[6];
	rot->at(6) = p[2] + s * p[6];
	rot->at(5) = p[5] + s * p[3];
	rot->at(7) = p[5] - s * p[3];

	rot->print();
*/
	int s = transpose ? 1. : -1.;
	rot->at(0,0) = 2.*q->q1()*q->q1() + 2.*q->q4()*q->q4() - 1.;
	rot->at(0,1) = 2.*q->q1()*q->q2() + 2.*s*q->q3()*q->q4();
	rot->at(0,2) = 2.*q->q1()*q->q3() - 2.*s*q->q2()*q->q4();
	rot->at(1,0) = 2.*q->q1()*q->q2() - 2.*s*q->q3()*q->q4();
	rot->at(1,1) = 2.*q->q2()*q->q2() + 2.*q->q4()*q->q4() - 1.;
	rot->at(1,2) = 2.*q->q2()*q->q3() + 2.*s*q->q1()*q->q4();
	rot->at(2,0) = 2.*q->q1()*q->q3() + 2.*s*q->q2()*q->q4();
	rot->at(2,1) = 2.*q->q2()*q->q3() - 2.*s*q->q1()*q->q4();
	rot->at(0,0) = 2.*q->q3()*q->q3() + 2.*q->q4()*q->q4() - 1.;

}

void Model::adjust_quats()
{
	for(int i=0 ; i<num_mols ; i++)
	{
		real qi = 1./sqrt( mols[i].q.quad() );
		mols[i].q *= qi;
	}
}

void Model::compute_site_coords()
{
	all_sites.clear();
	vector<Site> dummy;

	for(int i=0 ; i<num_mols ; i++)
	{
		Matrix rot;
		compute_rotation_mat( &rot , &mols[i].q , 1 );

		
		all_sites.push_back(dummy);

		for(int j=0 ; j<num_sites ; j++)
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

		mols[i].qo = mols[i].q;
		mols[i].qvo = mols[i].qv;
		
		pc(&mols[i],cr,cv,0);

		mols[i].ra2 = mols[i].ra1;
		mols[i].ra1 = mols[i].ra;

		mols[i].qa2 = mols[i].qa1;
		mols[i].qa1 = mols[i].qa;
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
	apply_thermostat();
	corrector_step();
	adjust_quats();
	boundary_conditions();
*/

	predictor_step();
	boundary_conditions();
	compute_forces();
	compute_bonded_forces();
	corrector_step();
	boundary_conditions();

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
//		evaluate_parameters();

//		export_vtk();

//		correct_temperature();
	}
//	export_file(&res);
}

void Model::gather_snapshots()
{
	for(int i=0 ; i<1000 ; i++)
	{
		progress_bar(i,1000);
		mols.clear();
		num_mols=0;
		initialize();
		simulate();
		add_to_snap_mat();
	}
	export_file();
}

void Model::add_to_snap_mat()
{
	vector<double> snap;
	for(int i=0 ; i<num_mols ; i++)
	{
		snap.push_back( mols[i].r.x() );
		snap.push_back( mols[i].r.y() );
		snap.push_back( mols[i].r.z() );
	}

	for(int i=0 ; i<num_mols ; i++)
	{
		snap.push_back( mols[i].rv.x() );
		snap.push_back( mols[i].rv.y() );
		snap.push_back( mols[i].rv.z() );
	}
	solutions.push_back(snap);
}


void Model::correct_temperature()
{
	real sum = 0;
	for(int i=0 ; i<num_mols ; i++)
		sum += mols[i].rv.quad3();

	real factor = vel_mag / sqrt( sum / num_mols );
	for(int i=0 ; i<num_mols ; i++)
		mols[i].rv *= factor;

	sum = 0;
	Vec w;
	for(int i=0 ; i<num_mols ; i++)
	{
		compute_ang_vel(i,&w);
		sum += mol_inert.x()*w.x()*w.x() + mol_inert.y()*w.y()*w.y() + mol_inert.z()*w.z()*w.z();
	}
	
	factor = vel_mag / sqrt( sum / num_mols );
	for(int i=0 ; i<num_mols ; i++)
		mols[i].qv *= factor;
}

void Model::apply_thermostat()
{
	real s1=0, s2=0;
	for(int i=0 ; i<num_mols ; i++)
	{
		s1 += mols[i].rv.inner( mols[i].ra );
		s2 += mols[i].rv.quad3();
	}

	Vec w;
	for(int i=0 ; i<num_mols ; i++)
	{
		compute_ang_vel(i,&w);
		s1 += w.inner( mols[i].torq );
		s2 += mol_inert.x()*w.x()*w.x() + mol_inert.y()*w.y()*w.y() + mol_inert.z()*w.z()*w.z();
	}
	real term = - s1/s2;
	for(int i=0 ; i<num_mols ; i++)
	{
		mols[i].ra += mols[i].rv * term;
		mols[i].qa += mols[i].qv * term;
	}
}

void Model::evaluate_parameters()
{
	Vec vel;
	real vel2 = 0;
	real vel3 = 0;


	for(int i=0 ; i<num_mols ; i++)
	{
		vel += mols[i].rv;
		vel2 += mols[i].rv.quad3();

		Vec w;
		compute_ang_vel(i,&w);
		vel3 += mol_inert.x()*w.x()*w.x() + mol_inert.y()*w.y()*w.y() + mol_inert.z()*w.z()*w.z();
	}
	
	real kinetic_energy = 0.5 * vel2 / num_mols;
	real rotational_kinetic_energy = 0.5 * vel3 / num_mols;
	real potential_energy = potential / num_mols;
	real total_energy = kinetic_energy + potential_energy + rotational_kinetic_energy;
	real pressure = density * (vel2 + virSum) / (num_mols * NDIM);
	
	cout << "summary of system:" << endl;
	cout.width(20); cout << right<< "kinetic energy : " << kinetic_energy << endl;
	cout.width(20); cout << right<< " rot kinetic energy : " << rotational_kinetic_energy << endl;
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

void Model::export_file()
{
	ofstream file;
	char path[] = "./output/data.txt";
	file.open(path);

	for(int i=0 ; i<solutions.size() ; i++)
	{
		for(int j=0 ; j<solutions[i].size() ; j++)
		{
			file << solutions[i][j] << " ";
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

void Model::progress_bar(int idx, int total)
{
	printf("                                               \r");
	real p = static_cast<double>(idx)/total;
	int p_int = static_cast<int>(p*20);
	cout << "[";
	for(int i=0 ; i<p_int ; i++)
		cout << "#";
	for(int j=p_int ; j<19 ; j++ )
		cout << " ";
	cout << "] - " << static_cast<int>(p*100 + 1) << "%";
	if(idx == total-1)
		cout << endl;

}

void Model::print()
{
	cout << "-----------------------------" << endl;
	for(int i=0 ; i<num_mols ; i++)
		cout << mols[i].r.x() << " " << mols[i].r.y() << " " << mols[i].r.z() << endl;
	cout << "-----------------------------" << endl;
}
