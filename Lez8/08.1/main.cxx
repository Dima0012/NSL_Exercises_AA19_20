#include"main.h"

int main(int argc, char* argv[])
{
	Initialize();

	for(int iopt=0;iopt<nopt;iopt++)
	{
		if(iopt%100==0) cout<<"Optimization step: "<<iopt<<endl;

		if(iopt==nopt-1) print=1;	//Print when optimitazione ended

		if(iopt!=0)	//Generating random parameters
		{
			mu_old = mu;
			sigma_old = sigma;
			mu = rnd.Rannyu(mu-delta_m,mu+delta_m);
			sigma = rnd.Rannyu(sigma-delta_s,sigma+delta_s);
		}

		for(int iblk=1;iblk<=nblock;iblk++)
		{
			Reset(iblk);
			for(int istep=1;istep<=nstep;istep++)
			{
				Move();
				Measure();
				Accumulate();
			}
			Averages(iblk);
		}

		H = glob_av[ie]/nblock;

		if(iopt==0) H_old = H;
		
		if(iopt!=0)	//Minimization
		{
			if(H<H_old)
			{
				accept_opt++;
				H_old = H;
				delta_m/=2.;
				delta_s/=2.;
			}
			else
			{
				mu = mu_old;
				sigma = sigma_old;
			}
		}
	}

	cout<<endl<<"Optimized parameters for trial wavefunction"<<endl;
	cout<<"   mu   ="<<mu<<endl;
	cout<<"   sigma ="<<sigma<<endl;
	cout<<"Accepted optimization steps "<<accept_opt<<" of "<<nopt<<endl;
	cout << "Mean acceptance rate for Metropolis: "<< acc_tot/(nblock*nopt) <<" %"<<endl;

	rnd.SaveSeed();

	return 0;
}


void Initialize(void)
{
	//===========Setting Random Generator ================

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

	cout<<"Variational Monte Carlo simulation for the Ground State 1D particle with external potential:"<<endl;
	cout<<"     V(x)=x^4-2.5x^2"<<endl;
	cout<<"using the trial function:"<<endl;
	cout<<"     Psi_T(x)=exp(-(x-mu)^2/2sigma^2)+exp(-(x+mu)^2/2sigma^2)"<<endl<<endl;


	//Setting parameters

	nblock = 200;
	nstep = 1000;
	nopt = 500;

	//Starting position
	x = 0.;			
	delta_x = 2.7;

	//Parameters for trial function
	mu = 1.;	
	delta_m = 1.;

	sigma = 1.;
	delta_s = 1.;

	//Acceptance for metropolis
	accepted = 0;
	attempted = 0;
	acc_tot = 0;
	accept_opt = 0;

	cout<<"Number of blocks "<<nblock<<endl;
	cout<<"Number of steps per block "<<nstep<<endl;
	cout<<"Initial parameters: mu="<<mu<<" sigma="<<sigma<<endl;
	cout<<"Optimization steps: "<<nopt<<endl;
	cout<<"of lenght: d_mu="<<delta_m<<" d_sigma="<<delta_s<<endl;
	cout<<"Starting position x = "<<x<<" with step lenght dx = "<<delta_x<<endl<<endl;


	//Array for measurament
	iv = 0; //Epot
	ik = 1; //Ekin
	ie = 2; //Etot

	n_props = 3;

	print = 0;

}


void Move()	//Sampling with Metropolis
{	
	ofstream pos;

	double x_new = rnd.Rannyu(x-delta_x, x+delta_x);

	double p_old = pow(Psi(x),2);
	double p_new = pow(Psi(x_new),2);

	double p = p_new/p_old;	//Trail prob
	double r = rnd.Rannyu();

	if(p>r){
		x = x_new;
		accepted++;
	}

	attempted++;

	if(print==1){
		pos.open("Data/position.0",ios::app);
		pos<<x<<endl;
		pos.close();
	}
}


void Measure()
{
	double v = V(x);
	double k = K(x);
	
	walker[iv] = v;
	walker[ik] = k;
	walker[ie] = v + k;
}



void Averages(int iblk)
{
	ofstream Epot, Ekin, Etot;
	const int wd = 12;

   	acc_tot += 100*(double)accepted/attempted;

	stima_epot = blk_average[iv]/blk_norm;
	glob_av[iv] += stima_epot;
	glob_av2[iv] += stima_epot*stima_epot;
	err_epot = Error(glob_av[iv],glob_av2[iv],iblk);
	
	stima_ekin = blk_average[ik]/blk_norm;
	glob_av[ik] += stima_ekin;
	glob_av2[ik] += stima_ekin*stima_ekin;
	err_ekin = Error(glob_av[ik],glob_av2[ik],iblk);
	
	stima_etot = blk_average[ie]/blk_norm;
	glob_av[ie] += stima_etot;
	glob_av2[ie] += stima_etot*stima_etot;
	err_etot = Error(glob_av[ie],glob_av2[ie],iblk);
	
	if(print==1){	//Printing only when optimitazion has ended
		Epot.open("Data/output.epot.0",ios::app);
		Ekin.open("Data/output.ekin.0",ios::app);
		Etot.open("Data/output.etot.0",ios::app);
		
		Epot<<setw(wd)<<iblk<<setw(wd)<<stima_epot<<setw(wd)<<glob_av[iv]/(double)iblk<<setw(wd)<<err_epot<<endl;
		Ekin<<setw(wd)<<iblk<<setw(wd)<<stima_ekin<<setw(wd)<<glob_av[ik]/(double)iblk<<setw(wd)<<err_ekin<<endl;
		Etot<<setw(wd)<<iblk<<setw(wd)<<stima_etot<<setw(wd)<<glob_av[ie]/(double)iblk<<setw(wd)<<err_etot<<endl;
		
		Epot.close();
		Ekin.close();
		Etot.close();
	}
}


void Reset(int iblk){ //Reset block averages
   if(iblk == 1){
		for(int i=0; i<n_props; ++i){
			glob_av[i] = 0;
			glob_av2[i] = 0;
		}
	}

	for(int i=0; i<n_props; ++i){
		blk_average[i] = 0;
	}
	blk_norm = 0;
	accepted = 0;
	attempted = 0;
}

void Accumulate(void){ //Update block averages
	for(int i=0; i<n_props; ++i){
		blk_average[i] = blk_average[i] + walker[i];
	}
	blk_norm = blk_norm + 1.0;
}


double Error(double sum, double sum2, int iblk)
{
    if( iblk == 1 ) return 0.;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

double Psi(double x)
{
	double a = pow((x-mu)/sigma,2);
	double b = pow((x+mu)/sigma,2);
	return exp(-a/2)+exp(-b/2);
}

double K(double x)
{
	double a = pow((x-mu)/sigma,2);
	double b = pow((x+mu)/sigma,2);
	return -0.5*((a-1)*exp(-a/2)+(b-1)*exp(-b/2))/(sigma*sigma);
}

double V(double x)
{
	return pow(x,4)-2.5*pow(x,2);
}

