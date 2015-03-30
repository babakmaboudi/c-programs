#ifndef MODEL_H
#define MODEL_H

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "particle.h"
#include "tools.h"

using namespace std;

class Model
{
	private:
		vector<Particle> par;

		double density;
	public:
		Model();

		void initiate(int);

		void save();
};

#endif
