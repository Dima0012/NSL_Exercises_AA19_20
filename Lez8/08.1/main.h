#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "random.h"

using namespace std;

//random gen
int seed[4];
Random rnd;

//observables
const int m_props=100;
int n_props, iv, ik, ie;
double walker[m_props];

//averages
double blk_average[m_props], glob_av[m_props], glob_av2[m_props], blk_norm;
int accepted, attempted, accept_opt;
double acc_tot;
double stima_epot, stima_ekin, stima_etot;
double err_epot, err_ekin, err_etot;
double H, H_old;

//position
double x, delta_x;

//parameters
double mu, sigma, mu_old, sigma_old, delta_m, delta_s;

//simulation
int nstep, nblock, nopt, print;

//functions
void Initialize(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
void Move(void);
void Measure(void);
double Error(double, double, int);
double Psi(double);
double K(double);
double V(double);