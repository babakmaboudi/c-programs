#ifndef VEC_H
#define VEC_H

#include <iostream>
#include <cmath>

using namespace std;

class Vec
{
	private:
		double x_val, y_val, z_val;
	public:
		Vec();
		Vec(double,double);
		Vec(double,double,double);

		void set(double,double);
		void set(double,double,double);
		double& x();
		double& y();
		double& z();

		void operator=(Vec);
		Vec operator+(Vec);
		void operator+=(Vec);
		Vec operator+(double);
		void operator+=(double);
		Vec operator-(Vec);
		void operator-=(Vec);
		Vec operator-(double);
		void operator-=(double);
		Vec operator*(Vec);
		void operator*=(Vec);
		Vec cross(Vec);
		Vec operator*(double);
		void operator*=(double);
		Vec operator/(Vec);
		void operator/=(Vec);
		Vec operator/(double);
		void operator/=(double);
		double inner(Vec);

		void rand_dir_2D();
		void rand_dir_3D();
		double norm2();
		double norm3();
		double quad2();
		double quad3();

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
		Vec_int operator+(double);
		void operator+=(double);
		Vec_int operator-(Vec_int);
		void operator-=(Vec_int);
		Vec_int operator-(double);
		void operator-=(double);
		Vec_int operator*(Vec_int);
		void operator*=(Vec_int);
		Vec_int operator*(double);
		void operator*=(double);
		Vec_int operator/(Vec_int);
		void operator/=(Vec_int);
		Vec_int operator/(double);
		void operator/=(double);

		int get_cell_num(Vec_int);

		void print();
};

#endif
