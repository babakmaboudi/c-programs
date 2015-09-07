#include "model.h"

Model::Model()
{
	NDIM = 3;
	num_mols_dir.set(5,5,5);
	num_mols = num_mols_dir.x() * num_mols_dir.y() * num_mols_dir.z();
	
	dt = 0.00125;
	temperature = 1.;
	density = 0.8;
	cut_off = pow(2.,1./6.);
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

void Model::single_step()
{
	leap_frog(1);
	boundary_conditions();
//	compute_forces();
	compute_forces_cell();
	leap_frog(2);
	
}

void Model::simulate()
{
	compute_forces();

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
