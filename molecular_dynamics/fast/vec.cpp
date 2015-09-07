#include "vec.h"

Vec::Vec()
{
	x_val = 0;
	y_val = 0;
	z_val = 0;
}

Vec::Vec(real x, real y)
{
	x_val = x;
	y_val = y;
}

Vec::Vec(real x, real y, real z)
{
	x_val = x;
	y_val = y;
	z_val = z;
}

void Vec::set(real x, real y)
{
	x_val = x;
	y_val = y;
}

void Vec::set(real x, real y, real z)
{
	x_val = x;
	y_val = y;
	z_val = z;
}

real& Vec::x()
{
	return x_val;
}

real& Vec::y()
{
	return y_val;
}

real& Vec::z()
{
	return z_val;
}

void Vec::operator=(Vec v)
{
	x_val = v.x_val;
	y_val = v.y_val;
	z_val = v.z_val;
}

Vec Vec::operator+(Vec v)
{
	Vec res;
	res.x_val = x_val + v.x_val;
	res.y_val = y_val + v.y_val;
	res.z_val = z_val + v.z_val;
	return res;	
}

void Vec::operator+=(Vec v)
{
	x_val += v.x_val;
	y_val += v.y_val;
	z_val += v.z_val;
}

Vec Vec::operator+(real c)
{
	Vec res;
	res.x_val = x_val + c;
	res.y_val = y_val + c;
	res.z_val = z_val + c;	
	return res;
}

void Vec::operator+=(real c)
{
	x_val += c;
	y_val += c;
	z_val += c;	
}

Vec Vec::operator-(Vec v)
{
	Vec res;
	res.x_val = x_val - v.x_val;
	res.y_val = y_val - v.y_val;
	res.z_val = z_val - v.z_val;
	return res;
}

void Vec::operator-=(Vec v)
{
	x_val -= v.x_val;
	y_val -= v.y_val;
	z_val -= v.z_val;
}

Vec Vec::operator-(real c)
{
	Vec res;
	res.x_val = x_val - c;
	res.y_val = y_val - c;
	res.z_val = z_val - c;	
	return res;	
}

void Vec::operator-=(real c)
{
	x_val -= c;
	y_val -= c;
	z_val -= c;
}

Vec Vec::operator*(Vec v)
{
	Vec res;
	res.x_val = x_val * v.x_val;
	res.y_val = y_val * v.y_val;
	res.z_val = z_val * v.z_val;
	return res;
}

void Vec::operator*=(Vec v)
{
	x_val *= v.x_val;
	y_val *= v.y_val;
	z_val *= v.z_val;
}

Vec Vec::operator*(real c)
{
	Vec res;
	res.x_val = x_val * c;
	res.y_val = y_val * c;
	res.z_val = z_val * c;	
	return res;	
}

void Vec::operator*=(real c)
{
	x_val *= c;
	y_val *= c;
	z_val *= c;
}

Vec Vec::operator/(Vec v)
{
	Vec res;
	res.x_val = x_val / v.x_val;
	res.y_val = y_val / v.y_val;
	res.z_val = z_val / v.z_val;
	return res;
}

void Vec::operator/=(Vec v)
{
	x_val /= v.x_val;
	y_val /= v.y_val;
	z_val /= v.z_val;	
}

Vec Vec::operator/(real c)
{
	Vec res;
	res.x_val = x_val / c;
	res.y_val = y_val / c;
	res.z_val = z_val / c;	
	return res;
}

void Vec::operator/=(real c)
{
	x_val /= c;
	y_val /= c;
	z_val /= c;
}

void Vec::rand_dir_2D()
{
	real s;
	s = 2. * M_PI * (double)rand() / (double)RAND_MAX;
	x_val = cos(s);
	y_val = sin(s);	
}

void Vec::rand_dir_3D()
{
	real s, x, y;

	s=2.;
	while(s>1.)
	{
		x = 2. * ((double)rand() / (double)RAND_MAX) - 1.;
		y = 2. * ((double)rand() / (double)RAND_MAX) - 1.;
		s = x*x + y*y;
	}
	z_val = 1. - 2. * s;
	s = 2. * sqrt(1. - s);
	x_val = s*x;
	y_val = s*y;
}

real Vec::norm2()
{
	return sqrt(x_val*x_val + y_val*y_val);
}

real Vec::norm3()
{
	return sqrt(x_val*x_val + y_val*y_val + z_val*z_val);
}

real Vec::quad2()
{
	return x_val*x_val + y_val*y_val;
}

real Vec::quad3()
{
	return x_val*x_val + y_val*y_val + z_val*z_val;
}

void Vec::print()
{
	cout << x_val << endl << y_val << endl << z_val << endl;
}

Vec_int::Vec_int()
{
	x_val = 0;
	y_val = 0;
	z_val = 0;
}

Vec_int::Vec_int(int x, int y)
{
	x_val = x;
	y_val = y;
}

Vec_int::Vec_int(int x, int y, int z)
{
	x_val = x;
	y_val = y;
	z_val = z;
}

void Vec_int::set(int x, int y)
{
	x_val = x;
	y_val = y;
}

void Vec_int::set(int x, int y, int z)
{
	x_val = x;
	y_val = y;
	z_val = z;
}

int& Vec_int::x()
{
	return x_val;
}

int& Vec_int::y()
{
	return y_val;
}

int& Vec_int::z()
{
	return z_val;
}

void Vec_int::operator=(Vec_int v)
{
	x_val = v.x_val;
	y_val = v.y_val;
	z_val = v.z_val;
}

void Vec_int::operator=(Vec v)
{
	x_val = v.x();
	y_val = v.y();
	z_val = v.z();
}

Vec_int Vec_int::operator+(Vec_int v)
{
	Vec_int res;
	res.x_val = x_val + v.x_val;
	res.y_val = y_val + v.y_val;
	res.z_val = z_val + v.z_val;
	return res;
}

void Vec_int::operator+=(Vec_int v)
{
	x_val += v.x_val;
	y_val += v.y_val;
	z_val += v.z_val;
}

Vec_int Vec_int::operator+(real c)
{
	Vec_int res;
	res.x_val = x_val + c;
	res.y_val = y_val + c;
	res.z_val = z_val + c;
	return res;
}

void Vec_int::operator+=(real c)
{
	x_val += c;
	y_val += c;
	z_val += c;
}

Vec_int Vec_int::operator-(Vec_int v)
{
	Vec_int res;
	res.x_val = x_val - v.x_val;
	res.y_val = y_val - v.y_val;
	res.z_val = z_val - v.z_val;
	return res;
}

void Vec_int::operator-=(Vec_int v)
{
	x_val -= v.x_val;
	y_val -= v.y_val;
	z_val -= v.z_val;
}

Vec_int Vec_int::operator-(real c)
{
	Vec_int res;
	res.x_val = x_val - c;
	res.y_val = y_val - c;
	res.z_val = z_val - c;
	return res;
}

void Vec_int::operator-=(real c)
{
	x_val -= c;
	y_val -= c;
	z_val -= c;
}

Vec_int Vec_int::operator*(Vec_int v)
{
	Vec_int res;
	res.x_val = x_val * v.x_val;
	res.y_val = y_val * v.y_val;
	res.z_val = z_val * v.z_val;
	return res;
}

void Vec_int::operator*=(Vec_int v)
{
	x_val *= v.x_val;
	y_val *= v.y_val;
	z_val *= v.z_val;
}

Vec_int Vec_int::operator*(real c)
{
	Vec_int res;
	res.x_val = x_val * c;
	res.y_val = y_val * c;
	res.z_val = z_val * c;
	return res;
}

void Vec_int::operator*=(real c)
{
	x_val *= c;
	y_val *= c;
	z_val *= c;
}

Vec_int Vec_int::operator/(Vec_int v)
{
	Vec_int res;
	res.x_val = x_val / v.x_val;
	res.y_val = y_val / v.y_val;
	res.z_val = z_val / v.z_val;
	return res;
}

void Vec_int::operator/=(Vec_int v)
{
	x_val /= v.x_val;
	y_val /= v.y_val;
	z_val /= v.z_val;
}

Vec_int Vec_int::operator/(real c)
{
	Vec_int res;
	res.x_val = x_val / c;
	res.y_val = y_val / c;
	res.z_val = z_val / c;
	return res;
}

void Vec_int::operator/=(real c)
{
	x_val /= c;
	y_val /= c;
	z_val /= c;
}

int Vec_int::get_cell_num(Vec_int v)
{
	return (z_val*v.y() + y_val)*v.x() + x_val;
}

void Vec_int::print()
{
	cout << x_val << endl << y_val << endl << z_val << endl;
}
