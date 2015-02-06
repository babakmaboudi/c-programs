#include <iostream>
#include "model.h"
#include "polynomial.h"

using namespace std;

int main(int argc, char** argv)
{
	char path[] = "param.txt";
	Model sat;
	sat.initiate_from_file(path);
	sat.SC_gaussian(5,4);
	return 0;
}
