#ifndef VEC_H
#define VEC_H

#include <iostream>
#include <vector>

using namespace std;

class Vec
{
	private:
		vector<double> elements;
	public:
		~Vec();
		void initiate(double*,int);
		void add_element(double);
		void clear();

		int get_size();
		double get_element(int);

		Vec operator+(Vec);
		Vec operator-(Vec);
		void operator=(Vec);

		Vec scalar(double);

		void print();


};

#endif
