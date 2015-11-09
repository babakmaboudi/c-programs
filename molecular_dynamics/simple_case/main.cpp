#include <iostream>
#include "model.h"

using namespace std;

int main(int argc, char** argv)
{
	srand(1);
	Model model;
	model.initialize();
//	model.simulate();
	model.gather_snapshots();
	return 0;
}
