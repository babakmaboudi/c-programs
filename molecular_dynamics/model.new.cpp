#include "model.h"

Model::Model()
{
	NDIM = 2;
	num_mols_dir.zeros(NDIM,1);
	num_mols_dir.at(0) = 20;
	num_mols_dir.at(1) = 20;
	num_mols = num_mols_dir.at(0) * num_mols_dir.at(1);
	
	dt = 0.005;
	temperature = 1.;
	density = 0.8;
	cut_off = pow(2.,1./6.);
	region = num_mols_dir * (1./sqrt(density));
	vel_mag = sqrt( NDIM * (1. - 1. / num_mols ) * temperature );
}

void Model::init_coords()
{
	Matrix gap(NDIM,1);
	gap.at(0) = region.at(0)/num_mols_dir.at(0);
	gap.at(1) = region.at(1)/num_mols_dir.at(1);
	
	for(int i=0 ; i<num_mols_dir.at(0) ; i++)
	{
		for(int j=0 ; j<num_mols_dir.at(1) ; j++)
		{
			Matrix c(NDIM,1);
			c.at(0) = i+0.5;
			c.at(1) = j+0.5;

			c.at(0) += gap.at(0)*c.at(0);
			c.at(1) += gap.at(1)*c.at(1);

			c += (region*(-.5));
		}
	}
}
