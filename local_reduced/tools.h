#ifndef TOOLS_H
#define TOOLS_H

#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include "../libraries/math/matrix.h"

using namespace std;

vector<double> linspace(double,double,int);
void arg_opt(double*,int*,Matrix*,char);

void prog_bar(int,int);

vector<double> get_polar_coord(Matrix*);

Matrix kmeans(Matrix*,int,int);
Matrix label_vectors(Matrix*,Matrix*);
int* randperm(int);

void find(vector<int>*,Matrix*,int);

#endif
