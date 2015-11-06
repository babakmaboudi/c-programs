#include <iostream>
#include "model.h"

using namespace std;

int main()
{
	srand(1);
	Model model;
//	model.build_reduced_basis();
//	model.simulate();
	model.simulate_reduced();
	return 0;
	
}
