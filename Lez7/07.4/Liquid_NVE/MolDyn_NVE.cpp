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
#include <algorithm>
#include <iomanip>
#include "MolDyn_NVE.h"

using namespace std;

int main(){
  Input();        //Inizialization
  if(stab) Stabilize();    //Equilibration if necessary
  else cout << "No stabilization perfomed" << endl;
  int nconf = 1;

  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
     Move();           //Move particles with Verlet algorithm
     //if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
     if(istep%nmeasure == 0){
        Measure();     //Properties measurement
        Accumulate(); //Update block averages
        //ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        nconf += 1;
     }
    }
    Averages(iblk);
  }

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
  ReadInput >> nblk;
  ReadInput >> nmeasure;
  ReadInput >> iprint;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps in each block = " << nstep << endl;
  cout << "Number of blocks = " << nblk << endl << endl;
  ReadInput.close();

  //Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 4; //Number of observables

  //measurement of g(r)
  igofr = 4;
  nbins = 100;
  n_props = n_props + nbins;
  bin_size = (box/2.0)/(double)nbins;

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

  double f_scale = sqrt(temp/new_temp); //Scaling velocities to match new temp

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

  int bin;
  double v, t, vij;
  double dx, dy, dz, dr;

  v = 0.0; //reset observables
  t = 0.0;

  //reset the hystogram of g(r)
  for (int k=igofr; k<igofr+nbins; ++k) walker[k]=0.0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( xold[i] - xold[j] ); // here I use old configurations [old = r(t)]
     dy = Pbc( yold[i] - yold[j] ); // to be compatible with EKin which uses v(t)
     dz = Pbc( zold[i] - zold[j] ); // => EPot should be computed with r(t)

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

//update of the histogram of g(r)
      for (bin = 0; bin < nbins; bin++)
      {
        if ( dr >= bin_size*bin && dr < bin_size*(bin + 1) )
        {
          walker[igofr+bin]+=2.0;
        }
      }

     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);

//Potential energy
       v += vij;
     }
    }          
  }

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
    walker[iv] = v/(double)npart; //Potential energy per particle
    walker[ik] = t/(double)npart; //Kinetic energy per particle
    walker[it] = (2.0 / 3.0) * t/(double)npart; //Temperature
    walker[ie] = (t+v)/(double)npart; //Total energy per particle

    return;
}


void Averages(int iblk) //Print results for current block
{
    
   double r;
   ofstream Gofr, Gave, Epot, Ekin, Etot, Temp;
   const int wd=12;
    
    //cout << "Block number " << iblk << endl;
    //cout << "Acceptance rate " << 100 * accepted/attempted << "%" << endl << endl;
    
    Epot.open("Data/output_epot.dat",ios::app);
    Ekin.open("Data/output_ekin.dat",ios::app);
    Temp.open("Data/output_temp.dat",ios::app);
    Etot.open("Data/output_etot.dat",ios::app);
    Gofr.open("Data/output.gofr.0",ios::app);
    Gave.open("Data/output.gave.0",ios::app);
    
    stima_epot = blk_av[iv]/blk_norm/(double)npart; //Potential energy
    glob_av[iv] += stima_epot;
    glob_av2[iv] += stima_epot*stima_epot;
    err_epot=Error(glob_av[iv],glob_av2[iv],iblk);
    
    stima_ekin = blk_av[ik]/blk_norm/(double)npart; //Kinetic energy
    glob_av[ik] += stima_ekin;
    glob_av2[ik] += stima_ekin*stima_ekin;
    err_ekin=Error(glob_av[ik],glob_av2[ik],iblk);

    stima_etot = blk_av[ie]/blk_norm/(double)npart;//Total energy
    glob_av[ie] += stima_etot;
    glob_av2[ie] += stima_etot*stima_etot;
    err_etot=Error(glob_av[ie],glob_av2[ie],iblk);

    stima_temp = blk_av[it]/blk_norm/(double)npart; //Temperature
    glob_av[it] += stima_temp;
    glob_av2[it] += stima_temp*stima_temp;
    err_temp=Error(glob_av[it],glob_av2[it],iblk);
    

//Potential energy per particle
    Epot << setw(wd) << iblk <<  setw(wd) << stima_epot << setw(wd) << glob_av[iv]/(double)iblk << setw(wd) << err_epot << endl;
//Kinetic Energy
    Ekin << setw(wd) << iblk <<  setw(wd) << stima_ekin << setw(wd) << glob_av[ik]/(double)iblk << setw(wd) << err_ekin << endl;
//Total Energy
    Etot << setw(wd) << iblk <<  setw(wd) << stima_etot << setw(wd) << glob_av[ie]/(double)iblk << setw(wd) << err_etot << endl;
//Temperature
    Temp << setw(wd) << iblk <<  setw(wd) << stima_temp << setw(wd) << glob_av[it]/(double)iblk << setw(wd) << err_temp << endl;

//g(r) (normalized)
    r = bin_size;

    for (int ibin=0; ibin<nbins; ibin++)
    {
      double delta_V = 4.*M_1_PI/3. * (pow( r + bin_size, 3) - pow(r, 3) );
      double norm = rho*npart*delta_V;  //Normalization

      stima_gdir = blk_av[ibin+igofr]/blk_norm/norm;
      glob_av[ibin+igofr] += stima_gdir;
      glob_av2[ibin+igofr] +=  stima_gdir*stima_gdir;
      err_gdir = Error(glob_av[ibin+igofr], glob_av2[ibin+igofr],iblk);

      Gofr << setw(wd) << iblk << setw(wd) << ibin << setw(wd) << stima_gdir << setw(wd) << glob_av[ibin+igofr]/(double)iblk <<endl;
    
      if(iblk==nblk)  //Final average 
      {
        Gave << setw(wd) << ibin << setw(wd) << r << setw(wd) << glob_av[ibin+igofr]/(double)iblk << setw(wd) << err_gdir << endl;
      }

      r+=bin_size;

    }

    //cout << "----------------------------" << endl << endl;

    Epot.close();
    Ekin.close();
    Etot.close();
    Temp.close();
    Gofr.close();
    Gave.close();
}

void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
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

double Error(double sum, double sum2, int iblk)
{
    if( iblk == 1 ) return 0.0;
    else return sqrt(abs(sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
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
