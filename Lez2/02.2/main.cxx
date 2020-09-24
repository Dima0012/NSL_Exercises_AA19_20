#include"lib.h"
#include"walker.h"	//Classe per rappresentare il camminatore in 3D

int main(int argc, char* argv[]) {

	//=========== Settaggio Random Generator ================

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

	int n_walk = 1E4; //Ripetizioni RW

	double a = 1.; //Lunghezza reticolo
	int T = 100; //Step temporali

	Walker *W = new Walker[n_walk];		//Camminatori reticolo
	Walker *W1 = new Walker[n_walk];	//Camminatori continuo

	double walk = 0.;

	double* value = new double[T];
	double* error = new double[T];	

	double* value1 = new double[T];
	double* error1 = new double[T];

	for(int i = 0; i<T; i++){	//Muovo ad ogni passo tutti i camminatori

		double* M = new double[n_walk];		//Modulo^2 vettore posizione reticolo
		double* M_sqrd = new double[n_walk];

		double* M1 = new double[n_walk];	//Modulo^2 vettore posizione continuo
		double* M1_sqrd = new double[n_walk];

		for(int j = 0; j<n_walk; j++){

			double A = rnd.Rannyu();
			double B = rnd.Rannyu();

			double theta = rnd.Rannyu(0., M_PI);
			double phi = rnd.Rannyu(0., 2.*M_PI);

			//Passo nel reticolo cubico

			if ( 0 <= A && A <= 0.5 )
				walk = -a;
			else walk = a;

			if ( 0. <= B && B < 1./3.)  W[j].Set_X( W[j].Get_X() + walk );

			if ( 1./3. <= B && B < 2./3.) W[j].Set_Y( W[j].Get_Y() + walk );	

			if ( 2./3. <= B && B < 1.) W[j].Set_Z( W[j].Get_Z() + walk );

			//Passo nel continuo

			W1[j].Set_X( W1[j].Get_X() + a*sin(theta)*cos(phi) ); 
			W1[j].Set_Y( W1[j].Get_Y() + a*sin(theta)*sin(phi) ); 
			W1[j].Set_Z( W1[j].Get_Z() + a*cos(theta) ); 

			M[j] = W[j].Get_X()*W[j].Get_X() + W[j].Get_Y()*W[j].Get_Y() + W[j].Get_Z()*W[j].Get_Z();			
			M1[j] = W1[j].Get_X()*W1[j].Get_X() + W1[j].Get_Y()*W1[j].Get_Y() + W1[j].Get_Z()*W1[j].Get_Z();

			M_sqrd[j] = M[j]*M[j];	
			M1_sqrd[j] = M1[j]*M1[j];
		}
	
		double sum = 0.;
		double sum1 = 0.;
		double sum_sqrd = 0.;
		double sum1_sqrd = 0.;		

		for(int j = 0; j<n_walk; j++){
			sum += M[j];
			sum1 += M1[j];
			sum_sqrd += M_sqrd[j];
			sum1_sqrd += M1_sqrd[j];
		}


		double mean = sum / n_walk;
		double mean1 = sum1 / n_walk;

		value[i] = sqrt( mean );
		value1[i] = sqrt( mean1 );

		error[i] = sqrt( sum_sqrd/n_walk - pow(mean, 2) );
		error1[i] = sqrt( sum1_sqrd/n_walk - pow(mean1, 2) );		//errore = dev std media

		error[i] *= 0.5 * (1./sqrt(mean));	//Per avere errore sulla radice del valore medio, usiamo la propagazione degli errori
		error1[i] *= 0.5 * (1./sqrt(mean1));


		delete[] M;
		delete[] M1;
		delete[] M_sqrd;
		delete[] M1_sqrd;
	}	

	DataPrint(value, T, "cubic_RW.out");
	DataPrint(error, T, "cubic_RW_err.out");

	DataPrint(value1, T, "continuum_RW.out");
	DataPrint(error1, T, "continuum_RW_err.out");

	delete[] value;
	delete[] value1;
	delete[] error;
	delete[] error1;

	rnd.SaveSeed();

	return 0;
}
