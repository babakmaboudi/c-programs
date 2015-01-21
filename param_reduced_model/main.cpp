#include <iostream>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include "matrix.h"
#include "model.h"
#include <cmath>
#include "tools.h"

using namespace std;

int main(int argc, char** argv)
{
	char file_path[] = "param.txt";
	Model model;
	model.initiate_from_file(file_path);
//	model.test_reduced_model();
	model.TRM_params();
/*
	vector<double> x;
	x.push_back(1);
	x.push_back(0);
	x.push_back(0);
	Matrix X;
	X.initiate_vector(x,3,1);

	vector<double> res = get_polar_coord(&X);
	for(int i = 0 ; i < 3 ; i++)
		cout << res[i] << " ";
	cout << endl;
*/
	return 0;
}
