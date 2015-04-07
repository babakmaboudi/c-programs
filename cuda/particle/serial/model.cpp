#include "model.h"

Model::Model()
{
	density = 0.0005;
	mass = 0.1;
	cutoff = 0.01;
	min_r = cutoff/100;
	dt = 0.0005;
}

void Model::initiate(int N)
{
	size = sqrt(density*N);
	int sx = (int)ceil(sqrt( static_cast<double>(N) ));
	int sy = (N+sx-1)/sx;


	int *rnd_perm = randperm(N);
	for(int i=0 ; i<N ; i++)
	{
		Particle p;
		int idx = rnd_perm[i];
		
		p.pos_x = sqrt(density*N)*(1.+(idx%sx))/(1+sx);
		p.pos_y = sqrt(density*N)*(1.+(idx/sx))/(1+sy);

		p.vel_x = drand48()*2-1;
		p.vel_y = drand48()*2-1;

		par.push_back(p);
	}
}

void Model::decompose_domain(int npatch)
{
	double patch_size = size / npatch;
	for(int i=0 ; i<npatch ; i++)
	{
		for(int j=0 ; j<npatch ; j++)
		{
			Patch pat;
			pat.index = i*npatch + j;
			pat.south = i*patch_size;
			pat.north = (i+1)*patch_size;
			pat.west = j*patch_size;
			pat.east = (j+1)*patch_size;
			pat.neighbours = find_patch_neighbors(i,j,npatch);
			patches.push_back(pat);
		}
	}

	for(int i=0 ; i<par.size() ; i++)
	{
		for(int j=0 ; j<patches.size() ; j++)
		{
			if( (par[i].pos_x<=patches[j].east) && (par[i].pos_x>patches[j].west) && (par[i].pos_y<=patches[j].north) && (par[i].pos_y>patches[j].south) )
				par[i].patch_number = j;
		}
	}

	for(int i=0 ; i<patches.size() ; i++)
	{
		for(int j=0 ; j<par.size() ; j++)
		{
			if( par[j].patch_number == i )
				patches[i].particles.push_back(j);
		}
	}
}

vector<int> Model::find_patch_neighbors(int i, int j, int npatch)
{
	vector<int> res;
	if( (i!=0) && (i!=npatch-1) && (j!=0) && (j!=npatch-1) )
	{
		res.push_back( (i-1)*npatch + j );
		res.push_back( (i+1)*npatch + j );
		res.push_back( i*npatch + (j-1) );
		res.push_back( i*npatch + (j+1) );
		res.push_back( (i-1)*npatch + (j-1) );
		res.push_back( (i-1)*npatch + (j+1) );
		res.push_back( (i+1)*npatch + (j-1) );
		res.push_back( (i+1)*npatch + (j+1) );
	}
	if( (i==0) && (j!=0) && (j!=npatch-1) )
	{
		res.push_back( i*npatch + (j-1) );
		res.push_back( i*npatch + (j+1) );
		res.push_back( (i+1)*npatch + (j-1) );
		res.push_back( (i+1)*npatch + (j+1) );
		res.push_back( (i+1)*npatch + j );
	}
	if( (i==npatch-1) && (j!=0) && (j!=npatch-1) )
	{
		res.push_back( i*npatch + (j-1) );
		res.push_back( i*npatch + (j+1) );
		res.push_back( (i-1)*npatch + (j-1) );
		res.push_back( (i-1)*npatch + (j+1) );
		res.push_back( (i-1)*npatch + j );
	}
	if( (j==0) && (i!=0) && (i!=npatch-1) )
	{
		res.push_back( (i-1)*npatch + j );
		res.push_back( (i+1)*npatch + j );
		res.push_back( (i+1)*npatch + (j+1) );
		res.push_back( (i-1)*npatch + (j+1) );
		res.push_back( i*npatch + (j+1) );
	}
	if( (j==npatch-1) && (i!=0) && (i!=npatch-1) )
	{
		res.push_back( (i-1)*npatch + j );
		res.push_back( (i+1)*npatch + j );
		res.push_back( (i+1)*npatch + (j-1) );
		res.push_back( (i-1)*npatch + (j-1) );
		res.push_back( i*npatch + (j-1) );
	}
	if( (i==0) && (j==0) )
	{
		res.push_back( (i+1)*npatch + j );
		res.push_back( i*npatch + (j+1) );
		res.push_back( (i+1)*npatch + (j+1) );
	}
	if( (i==0) && (j==npatch-1) )
	{
		res.push_back( i*npatch + (j-1) );
		res.push_back( (i+1)*npatch + j );
		res.push_back( (i+1)*npatch + (j-1) );
	}	
	if( (i==npatch-1) && (j==0) )
	{
		res.push_back( (i-1)*npatch + j );
		res.push_back( i*npatch + (j+1) );
		res.push_back( (i-1)*npatch + (j+1) );
	}
	if( (i==npatch-1) && (j==npatch-1) )
	{
		res.push_back( (i-1)*npatch + j );
		res.push_back( i*npatch + (j-1) );
		res.push_back( (i-1)*npatch + (j-1) );
	}
	return res;
}

void Model::apply_force(int idx, int nei)
{
	double dx = par[nei].pos_x - par[idx].pos_x;
	double dy = par[nei].pos_y - par[idx].pos_y;
	double r2 = dx * dx + dy * dy;
	if( r2 > cutoff*cutoff )
		return;

		
	r2 = fmax( r2, min_r*min_r );
	double r = sqrt( r2 );
 
    
	
	//
	//  very simple short-range repulsive force
	//
	double coef = ( 1 - cutoff / r ) / r2 / mass;
	par[idx].acc_x += coef * dx;
	par[idx].acc_y += coef * dy;	
}

void Model::move( int idx )
{
	//
	//  symplectic euler method to preserve eregy
	//
	par[idx].vel_x += par[idx].acc_x * dt;
	par[idx].vel_y += par[idx].acc_y * dt;
	par[idx].pos_x += par[idx].vel_x * dt;
	par[idx].pos_y += par[idx].vel_y * dt;

	//
	//  bounce from walls
	//
	while( par[idx].pos_x < 0 || par[idx].pos_x > size )
	{
		par[idx].pos_x  = (par[idx].pos_x < 0) ? -par[idx].pos_x : 2*size-par[idx].pos_x;
		par[idx].vel_x = -par[idx].vel_x;
	}
	while( par[idx].pos_y < 0 || par[idx].pos_y > size )
	{
		par[idx].pos_y  = (par[idx].pos_y < 0) ? -par[idx].pos_y : 2*size-par[idx].pos_y;
		par[idx].vel_y = -par[idx].vel_y;
	}
}

void Model::simulate(int MAX_ITER)
{
	for(int iter=0 ; iter<MAX_ITER ; iter++)
	{
		for(int i=0 ; i<patches.size() ; i++)
		{
			for(int j=0 ; j<patches[i].particles.size() ; j++)
			{
				int idx = patches[i].particles[j];
				par[idx].acc_x = 0;
				par[idx].acc_y = 0;

				for(int k=0 ; k<patches[i].neighbours.size() ; k++)
				{
					int neighbor_patch = patches[i].neighbours[k];
					for(int l=0 ; l<patches[neighbor_patch].particles.size() ; l++)
					{
						int nidx = patches[neighbor_patch].particles[l];
						apply_force(idx,nidx);
					}
				}

				for(int k=0 ; k<patches[i].particles.size() ; k++)
				{
					int nidx = patches[i].particles[k];
					apply_force(idx,nidx);
				}
			}
		}
		for(int i=0 ; i<par.size() ; i++)
			move(i);
		add_to_mat();
		for(int i=0 ; i<par.size() ; i++)
		{
			for(int j=0 ; j<patches.size() ; j++)
			{
				if( (par[i].pos_x<=patches[j].east) && (par[i].pos_x>patches[j].west) && (par[i].pos_y<=patches[j].north) && (par[i].pos_y>patches[j].south) )
					par[i].patch_number = j;
			}
		}
	
		for(int i=0 ; i<patches.size() ; i++)
		{
			patches[i].particles.clear();
			for(int j=0 ; j<par.size() ; j++)
			{
				if( par[j].patch_number == i )
					patches[i].particles.push_back(j);
			}
		}
	}
/*
	for(int iter=0 ; iter<MAX_ITER ; iter++)
	{
		for(int i=0 ; i<par.size() ; i++)
		{
			par[i].acc_x = 0;
			par[i].acc_y = 0;
			for(int j=0 ; j<par.size() ; j++)
			{
				apply_force(i,j);
			}
		}

		for(int i=0 ; i<par.size() ; i++)
			move(i);
		add_to_mat();
	}
*/
}

void Model::add_to_mat()
{
	Matrix row(1,2);
	for(int i=0 ; i<par.size() ; i++)
	{
		row.at(0) = par[i].pos_x;
		row.at(1) = par[i].pos_y;
		result.add_row(row);
	}
}

void Model::save()
{
	char path[] = "./data/positions.txt";
	result.save(path);
}
