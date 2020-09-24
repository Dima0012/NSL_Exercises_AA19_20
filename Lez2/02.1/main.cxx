#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include"integraleMC.h"
#include"funzione.h"
#include"lib.h"

using namespace std;
 
int main (int argc, char *argv[]){

	//===========Settaggio Random Generator ================

	Random rnd;
	int seed[4];
	int p1, p2;
	ifstream Primes("Primes");
	if (Primes.is_open()){
		Primes >> p1 >> p2 ;
	} else cerr << "PROBLEM: Unable to open Primes" << endl;
	Primes.close();

	ifstream input("seed.in");
	string property;
	if (input.is_open()){
		while ( !input.eof() ){
			input >> property;
			if( property == "RANDOMSEED" ){
				input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
				rnd.SetRandom(seed,p1,p2);
			}
		}
	input.close();
	} else cerr << "PROBLEM: Unable to open seed.in" << endl;


	//=======================================================

	
	function *f = new function;

	IntegraleMC intg(rnd);

	int n_step = 1E4;		//Numero step
	int n_blocks = 100;		//Numero blocchi
	int n_val = 1E3;		//Valutazioni integrale col metodo della media

	double IUnif[n_step];
	double ISamp[n_step];

	for(int i = 0; i<n_step; i++){
		IUnif[i] = intg.IntegraleAVE(0., 1., f, n_val);  //Sampling Unif.
		ISamp[i] = intg.IntegraleSAMP(0., 1., f, n_val); //Importance sampling
	}

	double* mean = new double[n_blocks];
	double* error = new double[n_blocks];

	MC_Ave_Err_Block(n_step, n_blocks, mean, error, IUnif);
	DataPrint(mean, n_blocks, "integral_uniform.out");		//Dati Uniforme
	DataPrint(error, n_blocks, "integral_uniform_err.out");


	MC_Ave_Err_Block(n_step, n_blocks, mean, error, ISamp);
	DataPrint(mean, n_blocks, "integral_sampling.out");		//Dati Importance Sampling
	DataPrint(error, n_blocks, "integral_sampling_err.out");

	delete[] mean;
	delete[] error;
	
	rnd.SaveSeed();

	return 0;
}
