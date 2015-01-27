#ifndef TOOLS_H
#define TOOLS_H

#include <iostream>
#include <vector>
#include <cmath>
#include "../libraries/math/matrix.h"

using namespace std;

vector<double> linspace(double,double,int);
void arg_opt(double*,int*,Matrix*,char);

void prog_bar(int,int);

#endif
