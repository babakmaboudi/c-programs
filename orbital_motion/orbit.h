#ifndef ORBIT_H
#define ORBIT_H

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>
#include <random>
#include "vec.h"
#include "basis.h"
#include "polynomial.h"

using namespace std;

class Orbit
{
	private:
		double mu;
		double G;
		double M;

		double a;
		double e;
		double i;
		double omega;
		double w;
		double nu;

		double* X0;
		double* V0;

	public:
		Orbit();
		void set_orbit_constants(double,double);
		void set_initial_position(double,double,double,double,double,double);
		void oe2sv();
		double* get_state_vec();

		void dynamic_2d(Vec*,double,Vec*);
		void dynamic_3d(Vec*,double,Vec*);

		void explicit_rk6(void (Orbit::*func)(Vec*,double,Vec*),double*,Vec,double,double*,double*,int);

		void monte_carlo(int);
		void SC_gaussian(int, int);

		double compute_hyper_phi(Basis*,vector<int>,vector<double>);
		double compute_hyper_lambda(vector<double>,vector<int>);
		double compute_hyper_weight(vector<double>,vector<int>);

		static void progress_bar(int,int);
		static void write_matrix(Vec*,double,double*,double*,int,int);
		static void write_matrix_file(double*,int,int,char*);
		static void write_2dvectord_file(vector< vector<double> >,char*);
		static void write_1dvectord_file(vector<double>,char*);
		void dummy();

};

#endif
