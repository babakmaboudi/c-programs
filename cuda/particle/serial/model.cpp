#include "model.h"

Model::Model()
{
	density = 0.0005;
}

void Model::initiate(int N)
{
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

void Model::save()
{
	char path[] = "./data/positions.txt";
	ofstream file;
	file.open(path);
	for(int i=0 ; i<par.size() ; i++)
	{
		file << par[i].pos_x << " " << par[i].pos_y << endl;
	}
}
