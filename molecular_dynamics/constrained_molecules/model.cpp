#include "model.h"

Model::Model()
{
	NDIM = 3;
	temperature = 4.17;
	density = 0.422;
	cut_off = pow(2., 1./6.);
	real c = 1. / pow (density, 1./3.);
	initUchain.set(3,3,3);
	region.set(initUchain.x()*c , initUchain.y()*c , initUchain.z()*c);
	num_chains = 2*initUchain.x()*initUchain.y()*initUchain.z();
	if(num_chains == 2)
		num_chains = 1;
	bond_lim = 2.1;
	chain_length = 4;
	vtk_count = 0;
	dt = 0.002;
	num_mols = initUcell.x()*initUcell.y()*initUcell.z() + num_chains*chain_length;
	vel_mag = sqrt( NDIM * (1. - 1. / num_mols ) * temperature );
	bond_lim = 2.1;
	bond_angle = 1.91063;
	bond_length = 0.39;
}

void Model::init_coords_chain()
{
	real by = bond_length * cos(bond_angle / 2.);
	real bz = bond_length * sin(bond_angle / 2.);

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
	num_constraints = 2 * chain_length - 3;
	cons_vec = new real[ num_constraints ];

	vel_mag = sqrt( ( NDIM * (1. - 1. / num_mols ) - static_cast<real> (num_constraints) / chain_length) * temperature );
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

void Model::build_constraint_matrix()
{
	m_mat.set(chain_length,num_constraints);
	for(int i=0 ; i<chain_length ; i++)
	{
		int m = 2*i-3;
		if(m >= 0)
			m_mat(i,m) = 2;
		m++;
		if(m >= 0)
			m_mat(i,m) = 2;
		m += 2;
		if(m < num_constraints)
			m_mat(i,m) = -2;
		m++;
		if(m < num_constraints)
			m_mat(i,m) = -2;
	}	
	for(int m=0 ; m<num_constraints ; m++)
	{
		Cons c;
		c.d2 = bond_length*bond_length;
		if(m%2==1)
			c.d2 *= 2. * (1. - cos(bond_angle));
		c.site1 = m/2;
		c.site2 = (m+3)/2;
		cons.push_back(c);
	}
}

void Model::initialize()
{
	init_coords_chain();
	assign_to_chain();
	init_vels();
	init_accels();
	build_constraint_matrix();
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
//			if(mols[i].chain == -1 || mols[i].chain != mols[j].chain || abs(i - j) > 3)
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

void Model::compute_chain_torsion_forces()
{
	real g[6] = {1.000, 1.310, -1.414, -0.330, 2.828, -3.394};
	real tCon = 15.50;

	for(int n=0 ; n<num_chains ; n++)
	{
		for(int i=0 ; i<chain_length-3 ; i++)
		{
			int nn = n*chain_length + i;
			Vec dr1 = mols[nn+1].r - mols[nn].r;
			wrap(&dr1);
			Vec dr2 = mols[nn+2].r - mols[nn+1].r;
			wrap(&dr2);
			Vec dr3 = mols[nn+3].r - mols[nn+2].r;
			wrap(&dr3);

			real c11 = dr1.quad3();
			real c12 = dr1.inner(dr2);
			real c13 = dr1.inner(dr3);
			real c22 = dr2.quad3();
			real c23 = dr2.inner(dr3);
			real c33 = dr3.quad3();

			real ca = c13 * c22 - c12 * c23;
			real cb1 = c11 * c22 - c12 * c12;
			real cb2 = c22 * c33 - c23 * c23;
			real cd = sqrt (cb1 * cb2);
			real c = ca / cd;
			
			real f = - tCon * (g[1] + (2. * g[2] + (3. * g[3] + (4. * g[4] + 5. * g[5] * c) * c) * c) * c);

			real t1 = ca;
			real t2 = c11 * c23 - c12 * c13;
			real t3 = - cb1;
			real t4 = cb2;
			real t5 = c13 * c23 - c12 * c33;
			real t6 = - ca;
			real cr1 = c12 / c22;
			real cr2 = c23 / c22;

			Vec w1 = (dr1*t1 + dr2*t2 + dr3*t3)*(f * c22 / (cd * cb1));
			Vec w2 = (dr1*t4 + dr2*t5 + dr2*t6)*(f * c22 / (cd * cb2));
			mols[nn].ra += w1;
			mols[nn + 1].ra += w1*( - (1. + cr1) );
			mols[nn + 1].ra += w2*( cr2 );
			mols[nn + 2].ra += w1*( cr1 );
			mols[nn + 2].ra += w2*( - (1. + cr2) );
			mols[nn + 3].ra += w2;
//			uSum += tCon * (g[0] + (g[1] + (g[2] + (g[3] + (g[4] + g[5] * c) * c) * c) * c) * c);
		}
	}
}

void Model::compute_chain_angel_forces()
{
	real cCon = cos (M_PI - bond_angle);
	real aCon = 868.6;

	for(int n=0 ; n<num_chains ; n++)
	{
		for(int i=0 ; i<chain_length ; i++)
		{
			int nn = n*chain_length + i;
			Vec dr1 = mols[nn+1].r - mols[nn].r;
			wrap(&dr1);
			Vec dr2 = mols[nn+2].r - mols[nn+1].r;
			wrap(&dr2);

			real c11 = dr1.quad3();
			real c12 = dr1.inner(dr2);
			real c22 = dr2.quad3();
			real cd = sqrt (c11 * c22);
			real c = c12 / cd;
			real f = - aCon * (c - cCon);
			
			Vec w1 = (dr1*(c12/c11) + dr2*(-1))*(f/cd);
			Vec w2 = (dr1 + dr2*(-c12/c22))*(f/cd);
			mols[nn].ra += w1;
			mols[nn + 1].ra -= w1;
			mols[nn + 1].ra -= w2;
			mols[nn + 2].ra += w2;
//			uSum += 0.5 * aCon * Sqr (c - cCon);
		}
	}
}

void Model::compute_constraints()
{
	for(int n=0 ; n<num_chains ; n++)
	{
		int nn = n * chain_length;
		for(int m=0 ; m<num_constraints ; m++)
		{
			cons[m].v = mols[nn + cons[m].site1].r - mols[nn + cons[m].site2].r;
			wrap(&cons[m].v);

		}

		Matrix l_mat(num_constraints,num_constraints);
		for(int m1=0 ; m1<num_constraints ; m1++)
		{
			for(int m2=0 ; m2<num_constraints ; m2++)
			{
				l_mat(m1,m2) = 0;
				int diff = m_mat(cons[m2].site1,m1) - m_mat(cons[m2].site2,m1);
				if(diff != 0)
				{
					l_mat(m1,m2) = diff * cons[m1].v.inner( cons[m2].v );
				}
			}
		}
		
		for(int m=0 ; m<num_constraints ; m++)
		{
			Vec dv = mols[ nn + cons[m].site1 ].rv - mols[ nn + cons[m].site2 ].rv;
			Vec da = mols[ nn + cons[m].site1 ].ra - mols[ nn + cons[m].site2 ].ra;
			cons_vec[m] = - da.inner( cons[m].v ) - dv.quad3();
		}
		lin_eq(l_mat,cons_vec,num_constraints);
		for(int m=0 ; m<num_constraints ; m++)
		{
			for(int i=0 ; i<chain_length ; i++)
			{
				real w = m_mat(i,m);
				if(w != 0.)
					mols[nn+i].ra = cons[m].v * ( w*cons_vec[m] );
			}
		}
	}
}

void Model::boundary_conditions()
{
	for(int i=0 ; i<num_mols ; i++)
		wrap(&mols[i].r);
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
	predictor_step();
	boundary_conditions();
	compute_forces();
//	compute_chain_angel_forces();
//	compute_chain_torsion_forces();
	corrector_step();
	boundary_conditions();
}

void Model::simulate()
{
	vector< vector<Mol> > res;
	res.push_back(mols);
	for(int i=0 ; i<100 ; i++)
	{
		single_step();
		res.push_back(mols);
		evaluate_parameters();
		export_vtk();
//		correct_temperature();
	}
}

void Model::evaluate_parameters()
{
	Vec vel;

	for(int i=0 ; i<num_mols ; i++)
	{
		vel += mols[i].rv;
	}

	cout << vel.quad3() << endl;
}

void Model::lin_eq(Matrix A_mat, real* x, int n)
{
	real* A = new real[ n*n ];
	real* v_max_inv = new real[n];
	int* ptr_max = new int[n];

	for(int i=0 ; i<n ; i++)
	{
		for(int j=0 ; j<n ; j++)
			A[i*n +j] = A_mat(j,i);
	}

	for(int i=0 ; i<n ; i++)
	{
		real v_max = 0.;
		for(int j=0 ; j<n ; j++)
		{
			real v = fabs( A[i+n*j] );
			if( v > v_max )
				v_max = v;
		}
		v_max_inv[i] = 1. / v_max;
	}

	for(int m=0 ; m<n ; m++)
	{
		real v_max = 0.;
		for(int i=m ; i<n ; i++)
		{
			for(int k=0 ; k<m ; k++)
				A[i + n*m] -= A[i + n*k] * A[k + n*m];
			real v = fabs( A[i + n*m] * v_max_inv[i] );
			if( v > v_max )
			{
				v_max = v;
				ptr_max[m] = i;
			}
		}

		if(m != ptr_max[m])
		{
			for(int k=0 ; k<n ; k++)
			{
				real v = A[ m+n*k ];
				A[ m+n*k ] = A[ ptr_max[m] + n*k ];
				A[ ptr_max[m] + n*k ] = v;
			}
			v_max_inv[ ptr_max[m] ] = v_max_inv[m];
		}
		for(int j=m+1 ; j<n ; j++)
		{
			for(int k=0 ; k<m ; k++)
				A[m + n*j] -= A[m+n*k] * A[k+n*j];
			A[m + n*j] /= A[m + n*m];
		}
	}

	for(int i=0 ; i<n ; i++)
	{
		real v = x[ ptr_max[i] ];
		x[ ptr_max[i] ] = x[i];
		x[i] = v;
		for(int j=0 ; j<i ; j++)
			x[i] -= A[i + n*j] * x[j];
		x[i] /= A[i + n*i];
	}
	for(int i=n-2 ; i >= 0 ; i--)
	{
		for(int j=i+1 ; j<n ; j++)
			x[i] -= A[i+n*j] * x[j];
	}

	delete[] A;
	delete[] v_max_inv;
	delete[] ptr_max;
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
