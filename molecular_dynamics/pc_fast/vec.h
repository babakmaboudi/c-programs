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
		Vec(real,real,real);

		void set(real,real);
		void set(real,real,real);
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
		Vec cross(Vec);
		Vec operator*(real);
		void operator*=(real);
		Vec operator/(Vec);
		void operator/=(Vec);
		Vec operator/(real);
		void operator/=(real);

		void rand_dir_2D();
		void rand_dir_3D();
		real norm2();
		real norm3();
		real quad2();
		real quad3();

		void print();
};

class Vec_int
{
	private:
		int x_val, y_val, z_val;
	
	public:
		Vec_int();
		Vec_int(int,int);
		Vec_int(int,int,int);

		void set(int,int);
		void set(int,int,int);
		int& x();
		int& y();
		int& z();

		void operator=(Vec_int);
		void operator=(Vec);
		Vec_int operator+(Vec_int);
		void operator+=(Vec_int);
		Vec_int operator+(real);
		void operator+=(real);
		Vec_int operator-(Vec_int);
		void operator-=(Vec_int);
		Vec_int operator-(real);
		void operator-=(real);
		Vec_int operator*(Vec_int);
		void operator*=(Vec_int);
		Vec_int operator*(real);
		void operator*=(real);
		Vec_int operator/(Vec_int);
		void operator/=(Vec_int);
		Vec_int operator/(real);
		void operator/=(real);

		int get_cell_num(Vec_int);

		void print();
};

#endif
