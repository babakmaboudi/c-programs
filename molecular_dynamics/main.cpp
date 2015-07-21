#include <iostream>
#include "vec.h"
#include "model.h"
#include "mol.h"

using namespace std;

int main(int argc, char** argv)
{
	srand(1);
	Model model;
	model.initialize();
	model.simulate();
//	model.evaluate_parameters();
	
	return 0;
}
