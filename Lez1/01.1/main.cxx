#include"random.h"

#include<cmath>
#include<fstream>
#include<iostream>

using namespace std;

double error(double* AV, double* AV2, int n) {	// Funzione che calcola l'errore statistico

	if (n==0 || n==1) return 0;				
	else return sqrt((AV2[n] - AV[n]*AV[n])/(n-1));
}

int main(){

	Random rnd;		//Set del Random Generator
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

	

//============Esercizio 01.1==========================

	int M = 1E6;		// Totale lanci
	int N = 100;		// Totale blocchi
	int L = int(M/N);	// Numero di lanci in ogni blocco


	double* r = new double[M];
	for(int i = 0; i<M; i++){	//Estrazione casuale del lancio
		r[i] = rnd.Rannyu();
	}

	double* ave = new double[N];
	double* ave2 = new double[N];

	double* sum_prog = new double[N]; 
	double* sum2_prog = new double[N];
	double* err_prog = new double[N];

	for (int i = 0; i<N; i++){
		ave[i] = 0.;
		ave2[i] = 0.;
		sum_prog[i] = 0.;			
		sum2_prog[i] = 0.;
		err_prog[i] = 0.;
	}	

	for (int i = 0; i<N; i++){	//Ciclo per calcolare la stima di r_i per ogni blocco
		double sum = 0;
		for(int j = 0; j<L; j++){
			int k = j + i*L;
			sum += r[k];
		}
		ave[i] = sum/L; //r_i per ogni blocco
		ave2[i] = ave[i]*ave[i]; //r_i^2 per ogni blocco

	}


	for (int i=0; i<N; i++){			//Ciclo per calcolare la media progressiva sui blocchi con errore
		for (int j=0; j<i+1; j++){
			sum_prog[i] += ave[j];		//Somma progressiva i-esima	
			sum2_prog[i] += ave2[j];	//Somma quadra progressiva i-esima	
		}
		sum_prog[i] /= (i+1); 	//Media cumulata fino ad i
		sum2_prog[i] /= (i+1); 	//Media quadra cumulata fino ad i
		err_prog[i] = error(sum_prog,sum2_prog,i); //Incertezza statistica
	}

	
	ofstream out1("Data/sum_prog_1.out");
	for(int i = 0; i<N; i++) out1<<sum_prog[i]<<endl;

	ofstream out2("Data/err_prog_1.out");				//Dati per l'analisi
	for(int i = 0; i<N; i++) out2<<err_prog[i]<<endl;

	out1.close();
	out2.close();


//============Esercizio 01.2==========================

	//Per questo seconda parte di esercizio si tengono le variabili definite precedentemente, opportunamente azzerate

	for (int i = 0; i<N; i++){
		ave[i] = 0.;
		ave2[i] = 0.;
		sum_prog[i] = 0.;			
		sum2_prog[i] = 0.;
		err_prog[i] = 0.;
	}	

	for(int i = 0; i<M; i++)	
		r[i] = rnd.Rannyu();	

	for (int i = 0; i<N; i++){		//Ciclo per calcolare la stima di sigma_i per ogni blocco
		double sum = 0;
		for(int j = 0; j<L; j++){
			int k = j + i*L;
			sum += (r[k] - 0.5)*(r[k]-0.5);
		}
		ave[i] = sum/L; //sigma_i stimata per ogni blocco
		ave2[i] = ave[i]*ave[i]; //sigma_i^2 per ogni blocco
	}

	for (int i=0; i<N; i++){			//Ciclo per calcolare la media progressiva sui blocchi con errore
		for (int j=0; j<i+1; j++){
			sum_prog[i] += ave[j];			
			sum2_prog[i] += ave2[j];		
		}
		sum_prog[i] /= (i+1); 	
		sum2_prog[i] /= (i+1); 	
		err_prog[i] = error(sum_prog,sum2_prog,i);
	}
	

	ofstream out3("Data/sum_prog_2.out");
	for(int i = 0; i<N; i++) out3<<sum_prog[i]<<endl;

	ofstream out4("Data/err_prog_2.out");				//Dati per l'analisi
	for(int i = 0; i<N; i++) out4<<err_prog[i]<<endl;

	out3.close();
	out4.close();


//============Esercizio 01.3==========================

	int a = 100;	// Sottointervalli di [0,1]
	int n = 1E4;	// Numero di lanci
	
	double* chi_2 = new double[a];

	for(int i = 0; i < a; i++){

		double* V = new double[n];	// Set degli n lanci casuali
		for(int j = 0; j<n; j++)
			V[j] = rnd.Rannyu();

		double* W = new double[a];	// Contatore per le frequenze osservate

		for(int j = 0; j<a; j++){
			W[j] = 0;
			for(int k = 0; k<n; k++)
				if( V[k] >= double(j)/a && V[k] < double(j+1)/a)	// Calcolo delle frequenze osservate
					W[j]++;

		} 		


		double chi_sum = 0;

		for(int k = 0; k<a; k++)
			chi_sum += pow( W[k] - (n/a), 2);		 

		chi_2[i] = chi_sum / (double(n)/a); 		//Calcolo del Chi_2 
		delete[] V;
		delete[] W;
	}

	ofstream out5("Data/chi_2.out");
	for(int i = 0; i<a; i++)		//Dati per l'analisi
		out5<<chi_2[i]<<endl;

	out5.close();

	rnd.SaveSeed();

return 0;
}
