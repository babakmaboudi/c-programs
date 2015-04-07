#ifndef MODEL_H
#define MODEL_H

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "particle.h"
#include "patch.h"
#include "tools.h"
#include "../../../libraries/math/matrix.h"

using namespace std;

class Model
{
	private:
		vector<Particle> par;
		vector<Patch> patches;

		double density;
		double cutoff;
		double mass;
		double min_r;
		double size;

		double dt;

		Matrix result;
	public:
		Model();

		void initiate(int);
		void decompose_domain(int);
		vector<int> find_patch_neighbors(int,int,int);

		void apply_force(int,int);
		void move(int);

		void simulate(int);

		void add_to_mat();
		void save();
};

#endif
