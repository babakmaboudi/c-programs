#include <iostream>
#include "wave.h"

using namespace std;

int main()
{
	Wave wave;
//	wave.build_reduced_model(500,1e-3);
	wave.test_reduced_model();
//	wave.solver(500);
	return 0;
}
