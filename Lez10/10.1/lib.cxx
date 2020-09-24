#include"lib.h"

void DataPrint(double* data, int n, string str) {

	string s = "Data/";
	s.append(str);

	ofstream out(s);

	for(int i = 0; i<n; i++)
		out<<data[i]<<endl;
	
	out.close();	
}

double error(double* AV, double* AV2, int n) {

	if (n==0 || n==1) return 0;				// Funzione che calcola l'errore statistico
	else return sqrt((AV2[n] - AV[n]*AV[n])/(n-1));
}


void MC_MeanProg(int n_throws, int n_blocks, double* mean, double* data) {	//Funzione per calcolare il valore medio a blocchi

	if(n_throws%n_blocks != 0){
		cerr<<"Simulation must be divided into n_blocks: please use a multiple of n_blocks as n_throws"<<endl;
		exit(-1);
	}

	int L = int(n_throws/n_blocks);

	double ave[n_blocks];

	for(int i = 0; i<n_blocks; i++)
		mean[i] = 0.;

	for (int i = 0; i<n_blocks; i++){	//Ciclo per calcolare la stima di data_i per ogni blocco
		double sum = 0;
		for(int j = 0; j<L; j++){
			int k = j + i*L;
			sum += data[k];
		}

		ave[i] = sum/L;
	}

	for(int i = 0; i<n_blocks; i++){
		for(int j = 0; j<i+1; j++)
			mean[i] += ave[j];
		mean[i] /= (i+1);
	}
}



void MC_ErrProg(int n_throws, int n_blocks, double* error, double* data){

	if(n_throws % n_blocks != 0){
		cerr<<"Simulation must be divided into n_blocks: please use a multiple of n_blocks as n_throws"<<endl;
		exit(-1);
	}

	int L = int(n_throws/n_blocks);	//Step in ogni blocco

	double ave[n_blocks];
	double ave2[n_blocks];
	double sum_prog[n_blocks];
	double sum2_prog[n_blocks];

	for(int i = 0; i<n_blocks; i++){
		ave[i] = 0.;
		ave2[i] = 0.;
		sum_prog[i] = 0.;
		sum2_prog[i] = 0.;
		error[i] = 0.;
	}

	for(int i = 0; i<n_blocks; i++){
		double sum = 0;
		for(int j = 0; j<L; j++){
			int pos = j + i*L;
			sum += data[pos];
		}
		ave[i] = sum/L;			//Media in ongi blocco
		ave2[i] = pow( ave[i], 2);	//Media quadra in ogni blocco
	}

	for(int i = 0; i<n_blocks; i++){
		for(int j = 0; j<i+1; j++){
			sum_prog[i] += ave[j];
			sum2_prog[i] += ave2[j];
		}

	sum_prog[i] /= (i+1);
	sum2_prog[i] /= (i+1);

	if(i != 0)
		error[i] =  sqrt( (sum2_prog[i] - pow( sum_prog[i], 2 )) / i );	//errore statistico progressivo

	error[0] = 0;

	}
}











