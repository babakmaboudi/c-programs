#ifndef MATRIX_OPERATORS_H
#define MATRIX_OPERATORS_H

#include <iostream>
#include <Accelerate/Accelerate.h>
#include <vector>
#include <cassert>

typedef __CLPK_integer integer;

using namespace std;

void svd(vector< vector<double> >,vector< vector<double> >*,vector<double>*,vector< vector<double> >*);

#endif
