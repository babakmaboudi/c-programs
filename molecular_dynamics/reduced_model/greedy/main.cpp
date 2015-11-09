#include <iostream>
#include "model.h"

using namespace std;

int main(int argc, char** argv)
{

	srand(0);
	Model model;
//	model.generate_snapshots();
//	model.generate_reduced_basis();
//	model.simulate();
	model.greedy();

//	model.evaluate_parameters();
	return 0;
}
