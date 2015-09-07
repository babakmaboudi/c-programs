#include "devicefunctions.h"

__global__ void find_cell(Particle* par, Cell* cell, int num_particles, int num_cells)
{
	int idx = threadIdx.x + blockIdx.x*blockDim.x;
	if(idx < num_particles)
	{
		for(int i=0 ; i<num_cells ; i++)
		{
			if( (par[idx].pos_x<=cell[i].east) && (par[idx].pos_x>cell[i].west) && (par[idx].pos_y<=cell[i].north) && (par[idx].pos_y>cell[i].south) )
					par[idx].cell_number = i;
		}
	}
	__syncthreads();
}

__global__ void find_particles(Particle* par, Cell* cell, int* tracker, int num_particles, int num_cells, int width)
{
	int idx = threadIdx.x + blockIdx.x*blockDim.x;
	if(idx == blockIdx.x*blockDim.x)
	{
		int cell_number = blockIdx.x;
		int ptr = 0;
		for(int i=0 ; i<num_particles ; i++)
		{
			if( par[i].cell_number == cell_number )
			{
				tracker[ cell_number*width + ptr ] = i;
				ptr++;
			}
		}
		cell[cell_number].num_particles = ptr;
	}
	__syncthreads();
}

__global__ void compute_force(Particle* par, Cell* cell, int* tracker, int num_particles, int num_cells, int width, int* neighbours)
{

	__shared__ float temp_part[TILTE_SIZE];
	__shared__ float acc[TILTE_SIZE];
	int idx = threadIdx.x;
	int num_cell = blockIdx.x;
	int par_num = tracker[ num_cell*width + idx ];
	float xx;
	float yy;

	if(idx < cell[num_cell].num_particles)
	{
		temp_part[idx*2] = par[ par_num ].pos_x;
		temp_part[idx*2+1] = par[ par_num ].pos_y;
	}

	__syncthreads();

	if(idx < cell[num_cell].num_particles)
	{
		acc[idx*2] = 0.0;
		acc[idx*2+1] = 0.0;
		xx = temp_part[idx*2];
		yy = temp_part[idx*2+1];
		for(int i=0 ; i<cell[num_cell].num_particles ; i++)
		{
			float ax = 0;
			float ay = 0;
			apply_force( &ax , &ay , xx , yy , temp_part[i*2] , temp_part[i*2+1] );

			acc[idx*2] += ax;
			acc[idx*2+1] += ay;
		}
	}
	__syncthreads();

	for(int nei=0 ; nei<8 ; nei++)
	{
		int nei_cell_num = neighbours[ num_cell*8 + nei ];
		if(idx < cell[nei_cell_num].num_particles)
		{
			int nei_num = tracker[ nei_cell_num*width + idx ];
			temp_part[idx*2] = par[ nei_num ].pos_x;
			temp_part[idx*2+1] = par[ nei_num ].pos_y;	
		}
		
		__syncthreads();

		if(idx < cell[num_cell].num_particles)
		{

			for(int i=0 ; i<cell[nei_cell_num].num_particles ; i++)
			{
				float ax = 0;
				float ay = 0;
				apply_force( &ax , &ay , xx , yy , temp_part[i*2] , temp_part[i*2+1] );

				acc[idx*2] += ax;
				acc[idx*2+1] += ay;
			}
//				acc[idx*2] += nei_cell_num;
//				acc[idx*2+1] += nei_cell_num;

		}	
	}

	if(idx < cell[num_cell].num_particles)
	{
		par[ par_num ].acc_x = acc[idx*2];
		par[ par_num ].acc_y = acc[idx*2+1];
	}
	__syncthreads();
}

__global__ void move(Particle* par ,int num_particles)
{
	int idx = threadIdx.x + blockIdx.x*blockDim.x;
	if(idx < num_particles)
	{
		par[idx].vel_x += par[idx].acc_x * dt;
		par[idx].vel_y += par[idx].acc_y * dt;
		par[idx].pos_x += par[idx].vel_x * dt;
		par[idx].pos_y += par[idx].vel_y * dt;


		// ------------------------------------------
		//
		//		CHANGE THIS PART
		//
		// ------------------------------------------

		float dsize = sqrt(density*num_particles);
		while( par[idx].pos_x < 0 || par[idx].pos_x > dsize )
		{
			par[idx].pos_x  = (par[idx].pos_x < 0) ? -par[idx].pos_x : 2*dsize-par[idx].pos_x;
			par[idx].vel_x = -par[idx].vel_x;
		}
		while( par[idx].pos_y < 0 || par[idx].pos_y > dsize )
		{
			par[idx].pos_y  = (par[idx].pos_y < 0) ? -par[idx].pos_y : 2*dsize-par[idx].pos_y;
			par[idx].vel_y = -par[idx].vel_y;
		}
	}
	__syncthreads();
}

__device__ void apply_force(float* acc_x, float* acc_y, float pos_x, float pos_y, float nei_x, float nei_y)
{
	float dx = nei_x - pos_x;
	float dy = nei_y - pos_y;
	float r2 = dx*dx + dy*dy;
	if( r2 > cutoff*cutoff )
		return;

	r2 = (r2 > min_r*min_r ) ? r2 : min_r*min_r;
	float r = sqrt(r2);

	float coef = ( 1 - cutoff/r ) / r2 / mass;
	(*acc_x) += coef * dx;
	(*acc_y) += coef * dy;
}
