#include <iostream>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include "reduced.h"
#include <vector>
#include "grid.h"
#include "node.h"
#include "matrix.h"

using namespace std;

int main(int argc, char** argv)
{
	char file_path[] = "param.txt";
	reduced model;
	model.initiate_from_file(file_path);
	model.test_reduced_model();
	return 0;
}
