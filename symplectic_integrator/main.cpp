#include <iostream>
#include "model.h"

using namespace std;

int main(int argc, char** argv)
{
	char path[] = "param.txt";
	Model model;
	model.initiate_from_file(path);
	model.single_sat();
	return 0;
}
