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
	model.BRM_params();
	return 0;
}
