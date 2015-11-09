#ifndef QUAT_H
#define QUAT_H

#include <iostream>
#include <cmath>
#include "vec.h"

class Quat
{
	private:
		real q1_val, q2_val, q3_val, q4_val;
	public:
		Quat();
		Quat(real,real,real,real);

		void set(real,real,real,real);

		real& q1();
		real& q2();
		real& q3();
		real& q4();

		void operator=(Quat);
		Quat operator+(Quat);
		void operator+=(Quat);
		Quat operator+(real);
		void operator+=(real);
		Quat operator-(Quat);
		void operator-=(Quat);
		Quat operator-(real);
		void operator-=(real);
		Quat operator*(Quat);
		void operator*=(Quat);
		Quat operator*(real);
		void operator*=(real);
		Quat operator/(real);
		void operator/=(real);

		real quad();
		real norm();

		void print();
};

#endif
