#include <iostream>
#include <ctime>
#include "model.h"

using namespace std;

int main(int argcc, char** argv)
{
	srand48( time( NULL ) );

	int num_part = 1000;
	int MAX_ITER = 1000;

	Model model;
	model.initiate(num_part);
	model.decompose_domain(5);

	clock_t begin = clock();
	model.simulate(MAX_ITER);

	model.save();

	clock_t end = clock();

	double elapsed_time = double(end - begin) / CLOCKS_PER_SEC;
	cout << elapsed_time << endl;
	return 0;
}
