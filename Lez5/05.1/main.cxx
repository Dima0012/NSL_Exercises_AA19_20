#include"lib.h"
#include"funzione.h"	//Classe per il calcolo delle funzioni d'onda
#include"position.h"	//Classe per la posizione nello spazio

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


	int n_step = 1E6;
	int n_block = 100;

	Position *P = new Position[n_step];	//Si parte dall'origine

	double r = 1.13; //Raggio per lo step Metropolis in unità a0
	int n_accept = 0; //Contattore per acceptance rate 50%

	WF *WaveFunction = new WF;
	
	// 1s = | 1 0 0 >

	for( int i = 1; i<n_step; i++){	//Metropolis

		double x1 = P[i-1].Get_X();
		double y1 = P[i-1].Get_Y();
		double z1 = P[i-1].Get_Z();

		double theta = rnd.Rannyu(0., M_PI);
		double phi = rnd.Rannyu(0., 2.*M_PI);	//Sampling uniforme nello spazio

		double x2 = x1 + r*sin(theta)*cos(phi);
		double y2 = y1 + r*sin(theta)*sin(phi);
		double z2 = z1 + r*cos(theta);

		double a = min(1., pow( WaveFunction->WF_1s(x2, y2, z2), 2) / pow( WaveFunction->WF_1s(x1,y1,z1), 2) );

		if(rnd.Rannyu()<=a){
			P[i].Set_X(x2);
			P[i].Set_Y(y2);
			P[i].Set_Z(z2);
			n_accept++;
		}
		else {
			P[i].Set_X(x1);
			P[i].Set_Y(y1);
			P[i].Set_Z(z1);
		}
	}


	ofstream out("Data/1s_pos.out");
	for (int i = 0; i < n_step; i++)	
		out << P[i].Get_X() << " "<< P[i].Get_Y() << " " << P[i].Get_Z() << endl;
	out.close();
	

	double a_rate = double(n_accept)/n_step * 100.;
	cout<<"Acceptance rate for 1s is "<<a_rate<<"%; best is around 50%"<<endl;

	double *R = new double[n_step];

	for(int i = 0; i<n_step; i++){
		R[i] = sqrt( P[i].Get_X()*P[i].Get_X() + P[i].Get_Y()*P[i].Get_Y() + P[i].Get_Z()*P[i].Get_Z());
	}

	double mean[n_block];
	double error[n_block];

	MC_Ave_Err_Block(n_step, n_block, mean, error, R);

	DataPrint(mean, error, n_block, "r_1s.out");


	// 2p = | 2 1 0 >

	r = 2.75;  
	n_accept = 0;

	for (int i = 0; i < n_step; i++){	//Partiamo da uno dei due lobi
		P[i].Reset(0,0,5.);
	}

	for( int i = 1; i<n_step; i++){

		double x1 = P[i-1].Get_X();
		double y1 = P[i-1].Get_Y();
		double z1 = P[i-1].Get_Z();

		double theta = rnd.Rannyu(0., M_PI);
		double phi = rnd.Rannyu(0., 2.*M_PI);

		double x2 = x1 + r*sin(theta)*cos(phi);
		double y2 = y1 + r*sin(theta)*sin(phi);
		double z2 = z1 + r*cos(theta);

		double a = min(1., pow( WaveFunction->WF_2p(x2, y2, z2), 2) / pow( WaveFunction->WF_2p(x1,y1,z1), 2) );

		if(rnd.Rannyu()<=a){
			P[i].Set_X(x2);
			P[i].Set_Y(y2);
			P[i].Set_Z(z2);
			n_accept++;
		}
		else {
			P[i].Set_X(x1);
			P[i].Set_Y(y1);
			P[i].Set_Z(z1);
		}
	}
	
	out.open("Data/2p_pos.out");
	for (int i = 0; i < n_step; i++)	
		out << P[i].Get_X() << " "<< P[i].Get_Y() << " " << P[i].Get_Z() << endl;
	out.close();

	a_rate = double(n_accept)/n_step * 100.;
	cout<<"Acceptance rate for 2p is "<<a_rate<<"%; best is around 50%"<<endl;

	for(int i = 0; i<n_step; i++){
		R[i] = sqrt( P[i].Get_X()*P[i].Get_X() + P[i].Get_Y()*P[i].Get_Y() + P[i].Get_Z()*P[i].Get_Z());
	}

	MC_Ave_Err_Block(n_step, n_block, mean, error, R);

	DataPrint(mean, error, n_block, "r_2p.out");


	//Multivariate distribution (Gaussian)

	// 1s = | 1 0 0 >

	n_accept = 0;
	double sigma = 0.78; //Sigma for Gauss Distribution

	for (int i = 0; i < n_step; i++)
	{
		P[i].Reset();
	}
	

	for( int i = 1; i<n_step; i++){	//Metropolis

		double x1 = P[i-1].Get_X();
		double y1 = P[i-1].Get_Y();
		double z1 = P[i-1].Get_Z();

		double x2 = rnd.Gauss(x1, sigma);
		double y2 = rnd.Gauss(y1, sigma);	//Sampling con Gauss
		double z2 = rnd.Gauss(z1, sigma);

		double a = min(1., pow( WaveFunction->WF_1s(x2, y2, z2), 2) / pow( WaveFunction->WF_1s(x1,y1,z1), 2) );

		if(rnd.Rannyu()<=a){
			P[i].Set_X(x2);
			P[i].Set_Y(y2);
			P[i].Set_Z(z2);
			n_accept++;
		}
		else {
			P[i].Set_X(x1);
			P[i].Set_Y(y1);
			P[i].Set_Z(z1);
		}
	}


	out.open("Data/1s_pos_gauss.out");
	for (int i = 0; i < n_step; i++)	
		out << P[i].Get_X() << " "<< P[i].Get_Y() << " " << P[i].Get_Z() << endl;
	out.close();
	

	a_rate = double(n_accept)/n_step * 100.;
	cout<<"Acceptance rate for 1s Gauss is "<<a_rate<<"%; best is around 50%"<<endl;


	for(int i = 0; i<n_step; i++){
		R[i] = sqrt( P[i].Get_X()*P[i].Get_X() + P[i].Get_Y()*P[i].Get_Y() + P[i].Get_Z()*P[i].Get_Z());
	}

	MC_Ave_Err_Block(n_step, n_block, mean, error, R);

	DataPrint(mean, error, n_block, "r_1s_gauss.out");


	// 2p = | 2 1 0 >

	for(int i = 0; i<n_step; i++){  //Partiamo dal centro di uno dei due lobi
		P[i].Reset(0, 0, 5);
	}

	n_accept = 0;
	sigma = 1.85; //Sigma for Gauss Distribution

	for( int i = 1; i<n_step; i++){	//Metropolis

		double x1 = P[i-1].Get_X();
		double y1 = P[i-1].Get_Y();
		double z1 = P[i-1].Get_Z();

		double x2 = rnd.Gauss(x1, sigma);
		double y2 = rnd.Gauss(y1, sigma);	//Sampling con Gauss
		double z2 = rnd.Gauss(z1, sigma);

		double a = min(1., pow( WaveFunction->WF_2p(x2, y2, z2), 2) / pow( WaveFunction->WF_2p(x1,y1,z1), 2) );

		if(rnd.Rannyu()<=a){
			P[i].Set_X(x2);
			P[i].Set_Y(y2);
			P[i].Set_Z(z2);
			n_accept++;
		}
		else {
			P[i].Set_X(x1);
			P[i].Set_Y(y1);
			P[i].Set_Z(z1);
		}
	}


	out.open("Data/2p_pos_gauss.out");
	for (int i = 0; i < n_step; i++)	
		out << P[i].Get_X() << " "<< P[i].Get_Y() << " " << P[i].Get_Z() << endl;
	out.close();
	

	a_rate = double(n_accept)/n_step * 100.;
	cout<<"Acceptance rate for 2p Gauss is "<<a_rate<<"%; best is around 50%"<<endl;

	for(int i = 0; i<n_step; i++){
		R[i] = sqrt( P[i].Get_X()*P[i].Get_X() + P[i].Get_Y()*P[i].Get_Y() + P[i].Get_Z()*P[i].Get_Z());
	}

	MC_Ave_Err_Block(n_step, n_block, mean, error, R);

	DataPrint(mean, error, n_block, "r_2p_gauss.out");

	delete[] R;

	rnd.SaveSeed();

	return 0;
}
