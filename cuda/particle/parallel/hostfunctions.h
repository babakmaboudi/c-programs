#ifndef HOSTFUNCTIONS_H
#define HOSTFUNCTIONS_H

#include <iostream>
#include <cmath>
#include "particle.h"
#include "cell.h"

#define density 0.0005
#define mass 0.1
#define cutoff 0.01
#define min_r cutoff/100
#define dt 0.0005

void dummy();

Particle* initiate(int);
Cell* decompose_domain(int,int,int**);
void find_neighbours(int,int,int,int*);

int* randperm(int n);

#endif
