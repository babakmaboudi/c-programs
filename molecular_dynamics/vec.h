#ifndef VEC_H
#define VEC_H

#include <iostream>
#include <cmath>

using namespace std;

typedef double real;

class Vec
{
	private:
		real x_val, y_val, z_val;
	public:
		Vec();
		Vec(real,real);

		void set(real,real);
		real& x();
		real& y();
		real& z();

		void operator=(Vec);
		Vec operator+(Vec);
		void operator+=(Vec);
		Vec operator+(real);
		void operator+=(real);
		Vec operator-(Vec);
		void operator-=(Vec);
		Vec operator-(real);
		void operator-=(real);
		Vec operator*(Vec);
		void operator*=(Vec);
		Vec operator*(real);
		void operator*=(real);
		Vec operator/(Vec);
		void operator/=(Vec);
		Vec operator/(real);
		void operator/=(real);

		void rand_dir_2D();
		real norm2();
		real norm3();
		real quad2();
		real quad3();

		void print();
};

#endif
