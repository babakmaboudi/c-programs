#include "model.h"

Model::Model()
{
	NDIM = 2;
	num_mols_dir.set(20,20);
	num_mols = num_mols_dir.x() * num_mols_dir.y();
	
	dt = 0.005;
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
			Vec c(i+0.5,j+0.5);
			c *= gap;
			c += region*(-.5);

			Mol m;
			m.r = c;
			mols.push_back(m);
		}
	}
}

void Model::init_vels()
{
	Vec sum;
	for(int i=0 ; i<num_mols ; i++)
	{
		mols[i].rv.rand_dir_2D();
		mols[i].rv *= vel_mag;
		sum += mols[i].rv;
	}
	for(int i=0 ; i<num_mols ; i++)
	{
		mols[i].rv -= sum*(1./num_mols);
	}
}

void Model::initialize()
{
	init_coords();
	init_vels();
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
		mols[i].ra.set(0.0,0.0);

	for(int i=0 ; i<num_mols-1 ; i++)
	{
		for(int j=i+1 ; j<num_mols ; j++)
		{
			Vec dr = mols[i].r - mols[j].r;
			wrap(&dr);
			real dr2 = dr.quad2();

			if( dr2 < cc )
			{
				real dr2i = 1./dr2;
				real dr6i = pow(dr2i,3);
				real force = 48. * dr6i * (dr6i - 0.5) * dr2i;
				mols[i].ra += dr*(force);
				mols[j].ra += dr*(-force);
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
	compute_forces();
	leap_frog(2);
	
}

void Model::simulate()
{
	vector< vector<Mol> > res;
	for(int i=0 ; i<1000 ; i++)
	{
		single_step();
		res.push_back(mols);
		evaluate_parameters();
	}
	export_file(&res);
}

void Model::evaluate_parameters()
{
	Vec v_sum;

	for(int i=0 ; i<num_mols ; i++)
	{
		v_sum += mols[i].rv;
	}
	v_sum.print();
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

void Model::export_file(vector< vector<Mol> >* res)
{
	ofstream file;
	char path[] = "./data.txt";
	file.open(path);

	for(int i=0 ; i<num_mols ; i++)
	{
		for(int j=0 ; j<res->size() ; j++)
		{
			file << ((*res)[j])[i].r.x() << " " << ((*res)[j])[i].r.y() << " ";
		}
		file << endl;
	}
	file.close();
}

void Model::print()
{
	cout << "-----------------------------" << endl;
	for(int i=0 ; i<num_mols ; i++)
		cout << mols[i].ra.x() << " " << mols[i].ra.y() << endl;
	cout << "-----------------------------" << endl;
}
