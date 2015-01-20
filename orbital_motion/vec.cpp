#include "vec.h"

Vec::~Vec()
{
	this->clear();
}

void Vec::initiate(double* v,int size)
{
	for(int i = 0 ; i < size ; i++)
		elements.push_back( v[i] );
}

void Vec::add_element(double v)
{
	elements.push_back(v);
}

void Vec::clear()
{
	while (!elements.empty())
		elements.pop_back();
}

int Vec::get_size()
{
	return elements.size();
}

double Vec::get_element(int index)
{
	return elements[ index ];
}

Vec Vec::operator+(Vec vec)
{
	Vec result;
	if(elements.size() != vec.get_size())
		cout << "Error, vector dimesions must agree!" << endl;
	else
	{
		for(int i = 0 ; i < elements.size() ; i ++)
			result.add_element( elements[i] + vec.get_element(i) );
	}
	return result;
}

Vec Vec::operator-(Vec vec)
{
	Vec result;
	if(elements.size() != vec.get_size())
		cout << "Error, vector dimesions must agree!" << endl;
	else
	{
		for(int i = 0 ; i < elements.size() ; i ++)
			result.add_element( elements[i] - vec.get_element(i) );
	}
	return result;
}

void Vec::operator=(Vec vec)
{
	this->clear();
	for(int i = 0 ; i < vec.get_size() ; i++)
		elements.push_back(vec.get_element(i));
}

Vec Vec::scalar(double s)
{
	Vec result;
	for(int i = 0 ; i < elements.size() ; i++)
	{
		result.add_element( s * elements[i] );
	}
	return result;
}

void Vec::print()
{
	for(int i = 0 ; i < elements.size() ; i++)
		cout << elements[i] << " ";
	cout << endl;
}


