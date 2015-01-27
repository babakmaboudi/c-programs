#include <iostream>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include "model.h"
#include "../libraries/math/matrix.h"
#include <cmath>

using namespace std;

int main(int argc, char** argv)
{
	char file_path[] = "param.txt";
	Model model;
	model.initiate_from_file(file_path);
	model.test_reduced_model();
	return 0;
}
