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

void Vec::set(real x, real y)
{
	x_val = x;
	y_val = y;
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
