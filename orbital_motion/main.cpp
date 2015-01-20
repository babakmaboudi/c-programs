#include <iostream>
#include <cmath>
#include "orbit.h"
#include "polynomial.h"
#include "basis.h"
#include "vec.h"

using namespace std;

int main()
{
	Orbit orbit;
	orbit.set_orbit_constants(6.67384e-11,5.972e24);
	orbit.set_initial_position(35000000,0.2,0,0,0,0);
	orbit.oe2sv();
	double* X0;
	X0 = orbit.get_state_vec();

//	orbit.monte_carlo(10000);

	orbit.SC_gaussian(5,4);
	
//	double p[11] = { -945, 0, 4725, -0, -3150, 0, 630, -0, -45, 0, 1 };
//	Polynomial pol;
//	pol.initiate_with_reference(10,p);
//	cout << "other : " << pol.evaluate(-4.85946282833231) << endl;
/*	
	int degree = 10;
	double* c = new double[11];
	for(int i = 0 ; i < 11 ; i++)
		c[i] = i;
	Polynomial pol1;
	pol1.initiate_with_reference(degree,c);
	pol1.set_coef(10,22);
	pol1.print();

	degree = 5;
	double* c2 = new double[11];
	for(int i = 0 ; i < 11 ; i++)
		c2[i] = i;
	Polynomial pol2;
	pol2.initiate_with_reference(degree,c2);
	pol2.set_coef(0,13);
	pol2.print();

	pol1 = pol2 + pol1;

	pol1.print();
	pol2.print();

	Basis b;
	b.initialize_hermite(10);
	Polynomial* pol3 = b.get_element(10);
	pol3->extract_roots();

	cout << "-------------------------------------" <<endl;	
*/
	return 0;
}
