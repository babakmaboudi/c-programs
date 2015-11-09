#include <iostream>
#include "model.h"

using namespace std;

int main(int argc, char** argv)
{
	srand(0);
	Model model;
	model.initiate();
//	model.sine_Gordon();
	model.DEI();
	return 0;
}
