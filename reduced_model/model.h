#ifndef MODEL_H
#define MODEL_H

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <cmath>
#include <vector>
#include "grid.h"
#include "matrix.h"
#include "tools.h"

using namespace std;

class Model
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

		Matrix A;
		Matrix A_tilde;
		Matrix Phi;
		Matrix Psi;
		Matrix P;
		Matrix F_tilde;
		vector<int> P_list;
	public:
		void sv2oe(double*,double*,double*);		// State vector to orbital elements 
		void discrete_empirical_interpolation(Matrix*,Matrix*);

		void initiate_from_file(char*);

		void convert_to_matrix(Matrix*,vector< vector<double> >,int,char);
		void write_file(vector<vector<double> >*,char*);
		void write_matrix(Matrix*,double,vector<vector<double> >*,vector<double>*);

		void dynamic_3d(Matrix*,double,Matrix*);
		void nonlinear(Matrix*,double,Matrix*);
		void dynamic_reduced(Matrix*,double,Matrix*);
		void compute_nonlinear_term(void (Model::*func)(Matrix*,double,Matrix*),vector<double>*,double,vector<double>);
		void explicit_rk6(void (Model::*func)(Matrix*,double,Matrix*),double*,Matrix,double,vector<vector<double> >*,vector<double>*,int);

		void build_reduced_model(int,double);
		void test_reduced_model();

		void single_sat(double*,int);
};

#endif
