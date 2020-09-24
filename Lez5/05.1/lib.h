#include<cmath>
#include<iostream>
#include<fstream>
#include<cstdlib>
#include"random.h"
#include<algorithm>

using namespace std;

void MC_Ave_Err_Block(int n_throws, int n_blocks, double* mean, double* error, double* data);   //Calcolo media e erroe a blocchi
void DataPrint(double* data, double* error,  int n, string str);    //Stampa i dati
