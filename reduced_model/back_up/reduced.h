#ifndef REDUCED_H
#define REDUCED_H

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <cmath>
#include <vector>
#include "vec.h"
#include "grid.h"
#include "matrix.h"
#include "matrix_operators.h"

using namespace std;

class reduced
{
	private: 
		double G; 		// Gravitational constant
		double M_earth; 	// Mass of the Earth
		double mu;
		double ae;		// Elipsoidal paramtere of the Earth
		double J2;		// The second hoarmonic of the Earth

		double* X0;		// Initial position of the satellite
		double* V0;		// Initial velocity of the satellite
		double* orb_vec;	// Orbital elements [a,e,i,omega,w,nu]

		double* X_moon;		// Posiion of the moon with respect to Geocenter
		double M_moon;		// Mass of the Moon

		double* u_sun;		// Direction of the radiation of the Sun
		double P_sun;		// Momentum flux of the sun
		double Cr;		// Reflectivity coefficient of the satellite
		double W_satellite;	// Weight of the satellite
	public:
		void sv2oe(double*,double*,double*);		// State vector to orbital elements 

		void initiate_from_file(char*);

		void dynamic_3d(Vec*,double,Vec*);
		void nonlinear(Vec*,double,Vec*);
		void compute_nonlinear_term(void (reduced::*func)(Vec*,double,Vec*),vector<double>*,double,vector<double>);
		void explicit_rk6(void (reduced::*func)(Vec*,double,Vec*),double*,Vec,double,vector<vector<double> >*,vector<double>*,int);

		void test_reduced_model();
		void build_reduced_model(Matrix*,int,double);
		vector< vector<double> >* build_random_space(int,int);

		void convert_to_matrix(Matrix*,vector< vector<double> >,int,char);
		void write_matrix(Vec*,double,vector<vector<double> >*,vector<double>*);
		void write_file(vector<vector<double> >*,char*);

		void single_sat();

		static vector<double> linspace(double,double,int);
};

#endif
