#ifndef NODE_H
#define NODE_H

#include <iostream>
#include <vector>

using namespace std;

class Node
{
	private:
		vector<double> points;
		vector<double> weights;
		Node* next;
		Node* pre;
	public:
		Node();
		void set_points(vector<double>);
		void set_weights(vector<double>);
		void set_next(Node*);
		void set_pre(Node*);

		vector<double> get_points();
		vector<double> get_weights();
		Node* get_next();
		Node* get_pre();
};

#endif
