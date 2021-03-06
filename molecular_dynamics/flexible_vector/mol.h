#ifndef MOL_H
#define MOL_H

#include "vec.h"

struct Mol
{
	Vec r, rv, ra, ra1, ra2, ro, rvo;		// ra# are the values in two prevous time steps of correction predictor method and ro and rvo are temporary storage for the old coordinates and velocities 
	int chain;
};

#endif
