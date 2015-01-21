#ifndef GRID_H
#define GRID_H

#include <iostream>
#include "node.h"
#include <vector>

using namespace std;

class Grid
{
	private:
		Node* root;
		int size;
	public:
		Grid();
		void add_node(vector<double>);

		vector< vector<double> > build_grid();
		vector< vector<double> > kron_product(vector< vector<double> > , vector<double>);

		int get_size();

		void print_grid();
		void print_vec(vector<double>);
};

#endif
