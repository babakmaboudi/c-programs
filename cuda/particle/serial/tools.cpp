#include "tools.h"

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
