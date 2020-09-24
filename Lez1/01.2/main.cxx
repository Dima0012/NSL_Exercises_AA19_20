#include"lib.h"


 
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


	int N[4]; //Ripetizioni somme
	N[0] = 1;
	N[1] = 2;
	N[2] = 10;
	N[3] = 100; 

	int n = 1E4; //Numero realizzazioni somme variabili aleatorie


	for(int k = 0; k<4; k++){	//Ciclo per ogni N della somma

		double* W = new double[n];

		//STD				
		for(int j = 0; j<n; j++)	//Dati per il dado standard (6 facce)
			W[j] = 0;

		for(int i = 0; i<n; i++){

			for(int j = 0; j<N[k]; j++)
				W[i] += (rnd.Rannyu() / N[k]);	//Distribuzione uniforme tra [0,1)
		}

		string s = "dataD6_";
		s.append(to_string(k));
		s.append(".out");
	
		DataPrint(W, n, s);	//Dati per analisi


		//EXP
		for(int j = 0; j<n; j++)	//Dado esponenziale
			W[j] = 0;

		for(int i = 0; i<n; i++){

			for(int j = 0; j<N[k]; j++)
				W[i] += (rnd.Exp(1.) / N[k]);	//Distribuzione Exp con media = 1
		}

		s = "dataEXP_";
		s.append(to_string(k));
		s.append(".out");
						//Dati per analisi
		DataPrint(W, n, s);


		//LORENTZ
		for(int j = 0; j<n; j++)	//Dado lorenziano
			W[j] = 0;

		for(int i = 0; i<n; i++){

			for(int j = 0; j<N[k]; j++)
				W[i] += (rnd.Lorentz(0., 1.) / N[k]);	//Distribuzione Lorentz con mu = 0, gamma = 1
		}

		s = "dataLOR_";
		s.append(to_string(k));
		s.append(".out");
						//Dati per analisi
		DataPrint(W, n, s);	

		delete[] W;
	}

	rnd.SaveSeed();

	return 0;
}
