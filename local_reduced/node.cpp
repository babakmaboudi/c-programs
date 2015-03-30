#include "node.h"

Node::Node()
{
	next = NULL;
	pre = NULL;
}

void Node::set_points(vector<double> p)
{
	points = p;
}

void Node::set_weights(vector<double> w)
{
	weights = w;
}

void Node::set_next(Node* n)
{
	next = n;
}

void Node::set_pre(Node* p)
{
	pre = p;
}

vector<double> Node::get_points()
{
	return points;
}

vector<double> Node::get_weights()
{
	return weights;
}

Node* Node::get_next()
{
	return next;
}

Node* Node::get_pre()
{
	return pre;
}
