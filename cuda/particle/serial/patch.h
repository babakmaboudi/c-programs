#ifndef PATCH_H
#define PATCH_H

#include <iostream>
#include <vector>

using namespace std;

struct Patch
{
	int index;
	double north;
	double south;
	double east;
	double west;
	vector<int> neighbours;
	vector<int> particles;
};

#endif
