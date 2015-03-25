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
//	model.single_sat();
/*
	Matrix_Complex mat;
	mat.rand(17,17);
	cout << mat << endl;
	mat.print();
*/
	return 0;
}
