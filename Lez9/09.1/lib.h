#include<cmath>
#include<iostream>
#include<fstream>
#include<cstdlib>
#include"random.h"
#include<algorithm>

using namespace std;

double error(double* , double* , int );
void MC_MeanProg(int, int , double*, double* );
void MC_ErrProg(int , int , double*, double* );
void DataPrint(double* , int , string );
