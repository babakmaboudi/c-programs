#include <iostream>
#include "model.h"

using namespace std;

int main()
{
	Model model;
	model.initialize();
	model.export_vtk();
	model.evaluate_parameters();
	model.simulate();
//	model.export_vtk();
//	model.dummy();
	return 0;
}
