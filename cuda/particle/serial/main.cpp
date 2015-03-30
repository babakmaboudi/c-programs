#include <iostream>
#include "model.h"

using namespace std;

int main(int argcc, char** argv)
{
	srand48( time( NULL ) );

	Model model;
	model.initiate(1000);
	model.save();
	return 0;
}
