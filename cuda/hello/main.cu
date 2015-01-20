#include<stdio.h>

__global__ void kernel(int* array)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;

	array[index] = index;
}

int main()
{
	int num_element = 1024;

	int* host_array = (int*)malloc( num_element * sizeof(int) );
	int* device_array;

	cudaMalloc( (void**)&device_array , num_element * sizeof(int) );

	int block_size = 128;
	int grid_size = num_element / block_size;

	kernel<<<grid_size,block_size>>>(device_array);

	cudaMemcpy(host_array,device_array,num_element * sizeof(int),cudaMemcpyDeviceToHost);

	for(int i = 0 ; i < num_element ; i++)
	{
		printf("%d ",host_array[i]);
	}
	printf("\n");
	
	free(host_array);
	cudaFree(device_array);
	return 0;
}
