#include <iostream>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include "model.h"
#include <cmath>

using namespace std;

int main(int argc, char** argv)
{
	char file_path[] = "param.txt";
	Model model;
	model.initiate_from_file(file_path);
	model.TRM_params();
	model.single_sat();
	return 0;
}
