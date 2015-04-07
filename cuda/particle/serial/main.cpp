#include <iostream>
#include "model.h"

using namespace std;

int main(int argcc, char** argv)
{
	srand48( time( NULL ) );

	int num_part = 100;
	int MAX_ITER = 1000;

	Model model;
	model.initiate(num_part);
	model.decompose_domain(10);
	model.simulate(MAX_ITER);
	model.save();
	return 0;
}
