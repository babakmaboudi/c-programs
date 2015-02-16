#include <iostream>
#include "model.h"
#include "../libraries/math/matrix.h"

using namespace std;

int main(int argc, char** argv)
{
/*
	char path[] = "param.txt";
	Model model;
	model.initiate_from_file(path);
//	model.test_reduced_model();
	model.single_sat();
*/
	Matrix mat1,mat2,mat3;
	mat1.ones(3,4);
	mat2.ones(1,4);
	mat2 = mat2*2;
	mat2.at(0,0) = 3;
	mat2.at(0,1) = 10;
	mat1.at(2,3) = 9;

	mat3 = mat1 * mat2.tr();
	cout << mat2 << endl << mat1 << endl;
	cout << mat3 << endl;
	return 0;
}
