#include"lib.h"

int main() {

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


	//Parametri esperimento
	
	double d = 7.0;	//Distanza reticolo
	double L = 4.5;	//Lunghezza ago


	int N_step = 1E4;  // Ripetizioni esperimento
	int N_throws = 1E5; // Lanci in ogni esperimento
	int N_blocks = 100; //# Blocchi
	int N_hit = 0;	//Contatore per i successi

	double PI[N_step];

	 
	for(int i = 0; i<N_step; i++){	
		for(int	j = 0; j<N_throws; j++)  //Lanci
		if( rnd.Rannyu(0., d/2.) <= (L/2.)*sin(rnd.Rannyu(0., 0.5*M_PI)) ) 	//Se l'ago interseca la linea
				N_hit++;

	PI[i] = (2.*L*N_throws)/(N_hit*d); //La stima di Pi_greco è data da N_throws lanci, per un totale di N_step stime
	N_hit = 0;		
	}	

	double* mean = new double[N_blocks];
	double* error = new double[N_blocks];

	MC_Ave_Err_Block(N_step, N_blocks, mean, error, PI);

	DataPrint(mean, N_blocks, "pi.out");
	DataPrint(error, N_blocks, "error_pi.out");

	delete[] mean;
	delete[] error;

	rnd.SaveSeed();

	return 0;
}
