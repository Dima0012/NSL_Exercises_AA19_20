#include<cmath>
#include<iostream>
#include<fstream>
#include<cstdlib>
#include"random.h"

using namespace std;

void MC_Ave_Err_Block(int n_throws, int n_blocks, double* mean, double* error, double* data);
void DataPrint(double* data, int n, string str);
