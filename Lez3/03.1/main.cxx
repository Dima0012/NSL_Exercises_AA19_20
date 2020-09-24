#include"lib.h"

int main(int argc, char* argv[]) {

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


	//Parametri 

	double S0 = 100.; //Prezzo dell'asset
	double T = 1.;	//Scadenza contratto
	double K = 100.;//Prezzo di strike
	double r = 0.1;	//Tasso di interesse senza rischio
	double sigma = 0.25; //Volatilità

	int n_inter = 100; //Numero di step per la discretizzazione


	int n_blocks = 100; // Blocchi
	int n_step = 1E4;   // Ripetizioni esperimento
	
	double CALL[n_step];
	double PUT[n_step];


	//Sampling diretto di S(T)

	for(int j = 0; j<n_step; j++){		

		double random = rnd.Gauss(0., T);
		double S_j = S0*exp( (r - sigma*sigma*0.5)*T + sigma*random*sqrt(T) );
		
		CALL[j] = exp(-r*T)*max(0., S_j - K);
		PUT[j] = exp(-r*T)*max(0., K - S_j);
	}
	

	double* mean = new double[n_blocks];
	double* error = new double[n_blocks];

	MC_Ave_Err_Block(n_step, n_blocks, mean, error, CALL);
	DataPrint(mean, n_blocks, "call_1.out");
	DataPrint(error, n_blocks, "call_1_err.out");			//Dati sampling diretto

	MC_Ave_Err_Block(n_step, n_blocks, mean, error, PUT);
	DataPrint(mean, n_blocks, "put_1.out");
	DataPrint(error, n_blocks, "put_1_err.out");


	//=========================================================


	//Sampling discreto 

	for(int j = 0; j<n_step; j++){

		double S_j = S0;	// Si parte da S0

		double help = 0.;

		for( int k = 1; k<n_inter; k++){	//Ciclo sulla discretizzazione per raggiungere S(T)
				
			double tk = double(k)/n_inter;
			double tk1 = double(k+1)/n_inter;

			help = S_j;
			S_j = help*exp( (r - sigma*sigma*0.5)*(tk1-tk) + sigma*rnd.Gauss(0.,1.)*sqrt(tk1-tk) ); // S_j -> S_(j+1)
		}
			
		CALL[j] = exp(-r*T)*max(0., S_j - K);
		PUT[j] = exp(-r*T)*max(0., K - S_j);
	}


	MC_Ave_Err_Block(n_step, n_blocks, mean, error, CALL);
	DataPrint(mean, n_blocks, "call_2.out");
	DataPrint(error, n_blocks, "call_2_err.out");			//Dati sampling discreto

	MC_Ave_Err_Block(n_step, n_blocks, mean, error, PUT);
	DataPrint(mean, n_blocks, "put_2.out");
	DataPrint(error, n_blocks, "put_2_err.out");

	delete[] mean;
	delete[] error;

	rnd.SaveSeed();

	return 0;
}
