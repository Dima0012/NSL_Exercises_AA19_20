/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include<algorithm>
#include "MolDyn_NVE.h"

using namespace std;

int main(){
  Input();        //Inizialization
  if(stab) Stabilize();    //Equilibration if necessary
  else cout << "No stabilization perfomed" << endl;
  int nconf = 1;

  for(int istep=1; istep <= nstep; ++istep){
     Move();           //Move particles with Verlet algorithm
     if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
     if(istep%nmeasure == 0){
        Measure();     //Properties measurement

	      int j = istep/nmeasure;   //Preparing data for blocking

	      ave_temp[j] = stima_temp;
      	ave_epot[j] = stima_pot;
      	ave_ekin[j] = stima_kin;
      	ave_etot[j] = stima_etot;

        //ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        nconf += 1;
     }
  }

  Blocking();

  if(stab){
  ConfFinal();         //Write final configuration to config.0 to restart
  ConfOldFinal();      //Write old final configuration to old.0 to restart

  } else ConfLast();  //Write last configuration to old.final and config.final


  return 0;
 
}


void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf,ReadConfOld;
  double fs;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator
  
  ReadInput.open("input.dat"); //Read input

  ReadInput >> stab;

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> nblocks;
  ReadInput >> nmeasure;
  ReadInput >> iprint;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl;
  cout << "Number of blocks = " << nblocks << endl << endl;
  ReadInput.close();

  //Read initial configuration
  cout << "Read initial configuration from file config.0 " << endl << endl;
  ReadConf.open("config.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();

  ReadConfOld.open("old.0");

  if ( ReadConfOld.is_open() ) { //If file is present, read from file

    //Read initial old configuration
   cout << "Read initial old configuration from file old.0 " << endl << endl;

   for (int i=0; i<npart; ++i){
      ReadConfOld >> xold[i] >> yold[i] >> zold[i];
      xold[i] = xold[i] * box;
      yold[i] = yold[i] * box;
      zold[i] = zold[i] * box;
    }
    ReadConfOld.close();

  } else {  //If no file to read from is present

    cout << "File old.0 not present" << endl;

    //Prepare initial random velocities
     cout << "Prepare random velocities with center of mass velocity equal to zero " << endl;
    double sumv[3] = {0.0, 0.0, 0.0};
    for (int i=0; i<npart; ++i){
        vx[i] = rand()/double(RAND_MAX) - 0.5;
       vy[i] = rand()/double(RAND_MAX) - 0.5;
       vz[i] = rand()/double(RAND_MAX) - 0.5;

       sumv[0] += vx[i];
       sumv[1] += vy[i];
       sumv[2] += vz[i];
    }

   for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart; //CM velocity with mass = 1
   double sumv2 = 0.0;
   for (int i=0; i<npart; ++i){
     vx[i] = vx[i] - sumv[0];
     vy[i] = vy[i] - sumv[1];
     vz[i] = vz[i] - sumv[2];   //New velocity with no CM component

     sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
   }
   sumv2 /= (double)npart;

   fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor
   for (int i=0; i<npart; ++i){  
     vx[i] *= fs;
     vy[i] *= fs;
     vz[i] *= fs;

     xold[i] = Pbc(x[i] - vx[i] * delta);
     yold[i] = Pbc(y[i] - vy[i] * delta);  //Computing r(t-dt)
     zold[i] = Pbc(z[i] - vz[i] * delta);
   }

   ConfOldFinal(); //Saving this random config as the old one
   
   }

   return;
}

void Stabilize(void) {   //Algorithm to stabilize system

  Move();

  double kine = 0.;

  for (int i = 0; i < npart; i++)
  {
    kine += 0.5*(vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
  }

  double new_temp = 2./3. * kine/(double)npart;

  double f_scale = sqrt(temp/new_temp); //Salaing velocities to match new temp

  for (int i = 0; i < npart; i++)
  {
    vx[i] *= f_scale;
    vy[i] *= f_scale;
    vz[i] *= f_scale;

    xold[i] = Pbc( x[i] - delta*vx[i] );
    yold[i] = Pbc( y[i] - delta*vy[i] ); //Novel old spatial config
    zold[i] = Pbc( z[i] - delta*vz[i] ); 

  }

}



void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}


double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}

void Measure(){ //Properties measurement
  double v, t, vij;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp;

  Epot.open("Data/output_epot.dat",ios::app);
  Ekin.open("Data/output_ekin.dat",ios::app);
  Temp.open("Data/output_temp.dat",ios::app);
  Etot.open("Data/output_etot.dat",ios::app);

  v = 0.0; //reset observables
  t = 0.0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( xold[i] - xold[j] ); // here I use old configurations [old = r(t)]
     dy = Pbc( yold[i] - yold[j] ); // to be compatible with EKin which uses v(t)
     dz = Pbc( zold[i] - zold[j] ); // => EPot should be computed with r(t)

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);

//Potential energy
       v += vij;
     }
    }          
  }

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
    stima_pot = v/(double)npart; //Potential energy per particle
    stima_kin = t/(double)npart; //Kinetic energy per particle
    stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    stima_etot = (t+v)/(double)npart; //Total energy per particle

    Epot << stima_pot  << endl;
    Ekin << stima_kin  << endl;
    Temp << stima_temp << endl;
    Etot << stima_etot << endl;


    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();

    return;
}


void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;

  cout << endl << "Print final configuration to file config.0 " << endl << endl;
  WriteConf.open("config.0");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();
  return;
}

void ConfOldFinal(void){ //Write final old configuration
  ofstream WriteConf;

  cout << "Print final old configuration to file old.0 " << endl << endl;
  WriteConf.open("old.0");

  for (int i=0; i<npart; ++i){
    WriteConf << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
  }
  WriteConf.close();
  return;
}

void ConfLast(void){ //Write last configuration
  ofstream WriteConf1;
  ofstream WriteConf2;

  cout << endl << "Print last old configuration to file old.final " << endl  << endl;
  cout << "Print last configuration to file config.final " << endl << endl;

  WriteConf1.open("old.final");
  WriteConf2.open("config.final");


  for (int i=0; i<npart; ++i){
    WriteConf1 << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
    WriteConf2 << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf1.close();
  WriteConf2.close();
  return;
}


void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}

void DataPrint(double* data, double* error,  int n, std::string str) { //Print data blocking files

	string s = "Data/";
	s.append(str);

	ofstream out(s, ios::app);

	for(int i = 0; i<n; i++)
		out << data[i] << " " << error[i] << endl;
	
	out.close();	
}


void MC_Ave_Err_Block(int n_throws, int n_blocks, double* mean, double* error, double* data){ // Data blocking

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
		mean[i] = 0.;
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
		for(int j = 0; j<i+1; j++)
			mean[i] += ave[j];
		mean[i] /= (i+1);		//Valore medio
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

void Blocking(void) {

  MC_Ave_Err_Block(nstep/nmeasure, nblocks, mean, error, ave_temp);
  DataPrint(mean, error, nblocks, "ave_temp.out");

  MC_Ave_Err_Block(nstep/nmeasure, nblocks, mean, error, ave_epot);
  DataPrint(mean, error, nblocks, "ave_epot.out");

  MC_Ave_Err_Block(nstep/nmeasure, nblocks, mean, error, ave_ekin);
  DataPrint(mean, error, nblocks, "ave_ekin.out");

  MC_Ave_Err_Block(nstep/nmeasure, nblocks, mean, error, ave_etot);
  DataPrint(mean, error, nblocks, "ave_etot.out");

}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
