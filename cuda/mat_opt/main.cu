#include <stdio.h>
#include <stdlib.h>

#define TILE_WIDTH 32

__global__ void mat_mul(float* Md, float* Nd, float* Pd)
{
	__shared__ float Mds[TILE_WIDTH*TILE_WIDTH];
	__shared__ float Nds[TILE_WIDTH*TILE_WIDTH];

	int bx = blockIdx.x;
	int by = blockIdx.y;
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	int Row = bx*TILE_WIDTH + tx;
	int Col = by*TILE_WIDTH + ty;

	int width = blockDim.x * gridDim.x;

	float temp = 0;
	for( int i = 0 ; i < width/TILE_WIDTH ; i++ )
	{
		Mds[tx*TILE_WIDTH+ty] = Md[ Row*width + (i*TILE_WIDTH + ty) ];
		Nds[tx*TILE_WIDTH+ty] = Nd[ (i*TILE_WIDTH + tx)*width + Col ];
		__syncthreads();

		for(int j = 0 ; j < TILE_WIDTH ; j++)
			temp += Mds[ tx*TILE_WIDTH + j ] * Nds[ j*TILE_WIDTH + ty ];
		__syncthreads();
	}
	Pd[Row*width+Col] = temp;
}

void print_mat(float* mat, int nrows, int ncols)
{
	for(int i = 0 ; i < nrows ; i++)
	{
		for(int j = 0 ; j < ncols ; j++)
			printf("%f ", mat[i*ncols+j]);
		printf("\n");
	}
}

void sequential(float* result, float* mat1, float* mat2, int size)
{
	float temp = 0;
	for(int i = 0 ; i < size ; i++)
	{
		for(int j = 0 ; j < size ; j++)
		{
			temp = 0;
			for (int k = 0 ; k < size ; k++)
				temp += mat1[i*size + k] * mat2[j*size+k];
			result[i*size + j] = temp;
		}
	}
}

int main()
{
	int nrows = 1024;
	int ncols = 1024;
	int num_bytes = nrows*ncols*sizeof(float);

	float* M = (float*)malloc(num_bytes);
	float* Md;
	cudaMalloc( (void**)&Md , num_bytes );

	float* N = (float*)malloc(num_bytes);
	float* Nd;
	cudaMalloc( (void**)&Nd , num_bytes );

	for(int i = 0 ; i < nrows ; i++)
	{
		for(int j = 0 ; j < ncols ; j++)
		{
			M[i*ncols + j] = 1.0;
			N[i*ncols + j] = 1.0;
		}
	}

	float* P = (float*)malloc(num_bytes);
	float* Pd;
	cudaMalloc( (void**)&Pd , num_bytes );

	cudaMemcpy(Md,M,num_bytes,cudaMemcpyHostToDevice);
	cudaMemcpy(Nd,N,num_bytes,cudaMemcpyHostToDevice);

	dim3 grid_size;
	grid_size.x = nrows/TILE_WIDTH;
	grid_size.y = ncols/TILE_WIDTH;

	dim3 block_size;
	block_size.x = TILE_WIDTH;
	block_size.y = TILE_WIDTH;

	cudaEvent_t start, stop;
	float elapsedTime;

	cudaEventCreate(&start);
	cudaEventRecord(start,0);

	mat_mul<<<grid_size,block_size>>>(Md,Nd,Pd);

	cudaEventCreate(&stop);
	cudaEventRecord(stop,0);
	cudaEventSynchronize(stop);

	cudaEventElapsedTime(&elapsedTime, start,stop);
	printf("Elapsed time gpu : %f ms\n" ,elapsedTime);

	cudaMemcpy(P,Pd,num_bytes,cudaMemcpyDeviceToHost);

	float* Ph = (float*)malloc(num_bytes);

	cudaEventCreate(&start);
	cudaEventRecord(start,0);

	sequential(Ph, M, N, nrows);

	cudaEventCreate(&stop);
	cudaEventRecord(stop,0);
	cudaEventSynchronize(stop);

	cudaEventElapsedTime(&elapsedTime, start,stop);
	printf("Elapsed time cpu : %f ms\n" ,elapsedTime);


	return 0;
}
