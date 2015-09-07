#include "quat.h"

Quat::Quat()
{
	q1_val = 0.0;
	q2_val = 0.0;
	q3_val = 0.0;
	q4_val = 0.0;
}

Quat::Quat(real q1, real q2, real q3,real q4)
{
	q1_val = q1;
	q2_val = q2;
	q3_val = q3;
	q4_val = q4;
}

void Quat::set(real q1, real q2, real q3,real q4)
{
	q1_val = q1;
	q2_val = q2;
	q3_val = q3;
	q4_val = q4;
}

real& Quat::q1()
{
	return q1_val;
}

real& Quat::q2()
{
	return q2_val;
}

real& Quat::q3()
{
	return q3_val;
}

real& Quat::q4()
{
	return q4_val;
}

void Quat::operator=(Quat q)
{
	q1() = q.q1();
	q2() = q.q2();
	q3() = q.q3();
	q4() = q.q4();
}

Quat Quat::operator+(Quat q)
{
	Quat res;
	res.q1() = q1() + q.q1();
	res.q2() = q2() + q.q2();
	res.q3() = q3() + q.q3();
	res.q4() = q4() + q.q4();
	return res;
}

void Quat::operator+=(Quat q)
{
	q1() += q.q1();
	q2() += q.q2();
	q3() += q.q3();
	q4() += q.q4();
}

Quat Quat::operator+(real c)
{
	Quat res;
	res.q1() = q1() + c;
	res.q2() = q2() + c;
	res.q3() = q3() + c;
	res.q4() = q4() + c;
	return res;
}

void Quat::operator+=(real c)
{
	q1() += c;
	q2() += c;
	q3() += c;
	q4() += c;
}

Quat Quat::operator-(Quat q)
{
	Quat res;
	res.q1() = q1() - q.q1();
	res.q2() = q2() - q.q2();
	res.q3() = q3() - q.q3();
	res.q4() = q4() - q.q4();
	return res;
}

void Quat::operator-=(Quat q)
{
	q1() -= q.q1();
	q2() -= q.q2();
	q3() -= q.q3();
	q4() -= q.q4();
}

Quat Quat::operator-(real c)
{
	Quat res;
	res.q1() = q1() - c;
	res.q2() = q2() - c;
	res.q3() = q3() - c;
	res.q4() = q4() - c;
	return res;
}

void Quat::operator-=(real c)
{
	q1() -= c;
	q2() -= c;
	q3() -= c;
	q4() -= c;
}

Quat Quat::operator*(Quat q)
{
	Quat res;
	res.q1() = q4()*q.q1() - q3()*q.q2() + q2()*q.q3() + q1()*q.q4();
	res.q2() = q3()*q.q1() + q4()*q.q2() - q1()*q.q3() + q2()*q.q4();
	res.q3() = - q2()*q.q1() + q1()*q.q2() + q4()*q.q3() + q3()*q.q4();
	res.q4() = - q1()*q.q1() - q2()*q.q2() - q3()*q.q3() + q4()*q.q4();
	return res;
}

void Quat::operator*=(Quat q)
{
	q1() = q4()*q.q1() - q3()*q.q2() + q2()*q.q3() + q1()*q.q4();
	q2() = q3()*q.q1() + q4()*q.q2() - q1()*q.q3() + q2()*q.q4();
	q3() = - q2()*q.q1() + q1()*q.q2() + q4()*q.q3() + q3()*q.q4();
	q4() = - q1()*q.q1() - q2()*q.q2() - q3()*q.q3() + q4()*q.q4();	
}

Quat Quat::operator*(real c)
{
	Quat res;
	res.q1() = q1() * c;
	res.q2() = q2() * c;
	res.q3() = q3() * c;
	res.q4() = q4() * c;
	return res;
}

void Quat::operator*=(real c)
{
	q1() *= c;
	q2() *= c;
	q3() *= c;
	q4() *= c;
}

Quat Quat::operator/(real c)
{
	Quat res;
	res.q1() = q1() / c;
	res.q2() = q2() / c;
	res.q3() = q3() / c;
	res.q4() = q4() / c;
	return res;
}

void Quat::operator/=(real c)
{
	q1() /= c;
	q2() /= c;
	q3() /= c;
	q4() /= c;
}

real Quat::quad()
{
	return q1()*q1() + q2()*q2() + q3()*q3() + q4()*q4();
}

real Quat::norm()
{
	return sqrt(q1()*q1() + q2()*q2() + q3()*q3() + q4()*q4());
}

void Quat::print()
{
	cout << q1() << endl << q2() << endl << q3() << endl << q4() << endl;
}
