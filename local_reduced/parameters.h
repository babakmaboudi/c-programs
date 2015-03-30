#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <iostream>
#include <vector>

using namespace std;

struct Parameters
{
	double G; 		// Gravitational constant
	double M_earth; 	// Mass of the Earth
	double mu;
	double ae;		// Elipsoidal paramtere of the Earth
	double J2;		// The second hoarmonic of the Earth

	vector<double> X0;		// Initial position of the satellite
	vector<double> V0;		// Initial velocity of the satellite
	vector<double> orb_vec;	// Orbital elements [a,e,i,omega,w,nu]

	vector<double> X_moon;		// Posiion of the moon with respect to Geocenter
	double M_moon;		// Mass of the Moon

	vector<double> u_sun;		// Direction of the radiation of the Sun
	double P_sun;		// Momentum flux of the sun
	double Cr;		// Reflectivity coefficient of the satellite
	
	double W_satellite;	// Weight of the satellite
};

#endif
