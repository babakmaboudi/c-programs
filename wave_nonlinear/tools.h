#ifndef TOOLS_H
#define TOOLS_H

#include <iostream>
#include <vector>
#include <cmath>
#include "../libraries/math/matrix.h"

using namespace std;

Matrix linspace(double,double,int);
void arg_opt(double*,int*,Matrix*,char);

void prog_bar(int,int);

vector<double> get_polar_coord(Matrix*);

Matrix compute_mixed_expansion(int,int); // computes the att the expansion degrees of a mixed polynomial of N th degree in d dimensions

#endif
