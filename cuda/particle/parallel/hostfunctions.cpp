#include "hostfunctions.h"

void dummy()
{
	cout << "hi" << endl;
}

Particle* initiate(int N)
{
	float size = sqrt( density * N );
	int sx = (int)ceil(sqrt( static_cast<float>(N) ));
	int sy = (N+sx-1)/sx;

	Particle* par = new Particle[N];
	int *rnd_perm = randperm(N);
	for(int i=0 ; i<N ; i++)
	{
		Particle p;
		int idx = rnd_perm[i];

		p.pos_x = sqrt(density*N)*(1.+(idx%sx))/(1+sx);
		p.pos_y = sqrt(density*N)*(1.+(idx/sx))/(1+sy);

		p.vel_x = drand48()*2-1;
		p.vel_y = drand48()*2-1;

		par[i] = p;
	}
	delete[] rnd_perm;
	return par;

}

Cell* decompose_domain(int num_particles, int num_cells, int**neighbour_cells)
{
	float size = sqrt( density * num_particles );
	float cell_width = size / num_cells;

	(*neighbour_cells) = new int[num_cells*num_cells*8];
	Cell* domain = new Cell[num_cells*num_cells];
	for(int i=0 ; i<num_cells ; i++)
	{
		for(int j=0 ; j<num_cells ; j++)
		{
			Cell cell;
			cell.north = (i+1)*cell_width;
			cell.east = (j+1)*cell_width;
			cell.south = i*cell_width;
			cell.west = j*cell_width;
			domain[i*num_cells + j] = cell;
			int* p = (*neighbour_cells) + (i*num_cells+j)*8;
			find_neighbours(i,j,num_cells,p);
		}
	}
	return domain;
}

void find_neighbours(int i, int j, int num_cells, int* neighbours)
{
	if( (i!=0) && (i!=num_cells-1) && (j!=0) && (j!=num_cells-1) )
	{
		neighbours[0] = (i+1)*num_cells + j;	// north
		neighbours[1] = (i+1)*num_cells + (j+1);	// north east
		neighbours[2] = i*num_cells + (j+1);	// east
		neighbours[3] = (i-1)*num_cells + (j+1);	// south east
		neighbours[4] = (i-1)*num_cells + j;	// south
		neighbours[5] = (i-1)*num_cells + (j-1);	// south west
		neighbours[6] = i*num_cells + (j-1);	// west
		neighbours[7] = (i+1)*num_cells + (j-1);	// north west
	}
	else if( (i==0) && (j!=0) && (j!=num_cells-1) )
	{
		neighbours[0] = (i+1)*num_cells + j;
		neighbours[1] = (i+1)*num_cells + (j+1);
		neighbours[2] = i*num_cells + (j+1);
		neighbours[3] = (num_cells-1)*num_cells + (j+1);
		neighbours[4] = (num_cells-1)*num_cells + j;;
		neighbours[5] = (num_cells-1)*num_cells + (j-1);;
		neighbours[6] = i*num_cells + (j-1);
		neighbours[7] = (i+1)*num_cells + (j-1);
	}
	else if( (i==num_cells-1) && (j!=0) && (j!=num_cells-1) )
	{
		neighbours[0] = j;
		neighbours[1] = (j+1);
		neighbours[2] = i*num_cells + (j+1);
		neighbours[3] = (i-1)*num_cells + (j+1);
		neighbours[4] = (i-1)*num_cells + j;
		neighbours[5] = (i-1)*num_cells + (j-1);
		neighbours[6] = i*num_cells + (j-1);
		neighbours[7] = (j-1);
	}
	else if( (j==0) && (i!=0) && (i!=num_cells-1) )
	{
		neighbours[0] = (i+1)*num_cells + j;
		neighbours[1] = (i+1)*num_cells + (j+1);
		neighbours[2] = i*num_cells + (j+1);
		neighbours[3] = (i-1)*num_cells + (j+1);
		neighbours[4] = (i-1)*num_cells + j;
		neighbours[5] = (i-1)*num_cells + (num_cells-1);
		neighbours[6] = i*num_cells + (num_cells-1);
		neighbours[7] = (i+1)*num_cells + (num_cells-1);
	}
	else if( (j==num_cells-1) && (i!=0) && (i!=num_cells-1) )
	{
		neighbours[0] = (i+1)*num_cells + j;
		neighbours[1] = (i+1)*num_cells;
		neighbours[2] = i*num_cells;
		neighbours[3] = (i-1)*num_cells;
		neighbours[4] = (i-1)*num_cells + j;
		neighbours[5] = (i-1)*num_cells + (j-1);
		neighbours[6] = i*num_cells + (j-1);
		neighbours[7] = (i+1)*num_cells + (j-1);
	}
	else if( (i==0) && (j==0) )
	{
		neighbours[0] = (i+1)*num_cells + j;
		neighbours[1] = (i+1)*num_cells + (j+1);
		neighbours[2] = i*num_cells + (j+1);
		neighbours[3] = (num_cells-1)*num_cells + (j+1);
		neighbours[4] = (num_cells-1)*num_cells + j;
		neighbours[5] = (num_cells-1)*num_cells + (num_cells-1);
		neighbours[6] = i*num_cells + (num_cells-1);
		neighbours[7] = (i+1)*num_cells + (num_cells-1);	
	}
	else if( (i==0) && (j==num_cells-1) )
	{
		neighbours[0] = (i+1)*num_cells + j;
		neighbours[1] = (i+1)*num_cells;
		neighbours[2] = i*num_cells;
		neighbours[3] = (num_cells-1)*num_cells;
		neighbours[4] = (num_cells-1)*num_cells + j;
		neighbours[5] = (num_cells-1)*num_cells + (j-1);
		neighbours[6] = i*num_cells + (j-1);
		neighbours[7] = (i+1)*num_cells + (j-1);	
	
	}
	else if( (i==num_cells-1) && (j==0) )
	{
		neighbours[0] = j;
		neighbours[1] = (j+1);
		neighbours[2] = i*num_cells + (j+1);
		neighbours[3] = (i-1)*num_cells + (j+1);
		neighbours[4] = (i-1)*num_cells + j;
		neighbours[5] = (i-1)*num_cells + (num_cells-1);
		neighbours[6] = i*num_cells + (num_cells-1);
		neighbours[7] = (num_cells-1);

	}
	else//( (i==num_cells-1) && (j==num_cells-1) )
	{	
		neighbours[0] = j;
		neighbours[1] = 0;
		neighbours[2] = i*num_cells;
		neighbours[3] = (i-1)*num_cells;
		neighbours[4] = (i-1)*num_cells + j;
		neighbours[5] = (i-1)*num_cells + (j-1);
		neighbours[6] = i*num_cells + (j-1);
		neighbours[7] = (j-1);
	}
}

int* randperm(int n)
{
	int* perm = new int[n];
	for(int i = 0 ; i < n ; i++)
	{
		perm[i] = i;
	}

	int temp,r;
	for(int i = 0 ; i < n ; i++)
	{
		r = rand()%(n-i) + i;
		temp = perm[r];
		perm[r] = perm[i];
		perm[i] = temp;
	}
	return perm;
}


