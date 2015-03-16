#include <iostream>
#include "model.h"

using namespace std;

int main(int argc, char** argv)
{
	srand(0);
	char path[] = "param.txt";
	Model model;
	model.initiate_from_file(path);
	model.test_reduced_model();
//	model.build_reduced_model(5);
//	model.single_sat();
	return 0;
}
