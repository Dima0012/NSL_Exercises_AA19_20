/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
//parameters, observables
const int m_props=1000;
int n_props,nbins;
double bin_size;
int iv,ik,it,ie, igofr;
double stima_epot, stima_ekin, stima_etot, stima_temp, stima_gdir;
double walker[m_props];

// averages
double acc,att;
double blk_av[m_props],blk_norm,accepted,attempted;
double glob_av[m_props],glob_av2[m_props];
double err_ekin, err_epot, err_etot, err_temp, err_gdir;

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;

// simulation
int nstep, nblk, nmeasure;
int iprint, seed;
double delta;
bool stab;

//functions
void Input(void);
void Stabilize(void);
void Move(void);
void ConfFinal(void);
void ConfOldFinal(void);
void ConfLast(void);
void ConfXYZ(int);
void Accumulate(void);
void Averages(int iblk);
void Reset(int iblk);
void Measure(void);
double Force(int, int);
double Pbc(double);
double Error(double sum, double sum2, int iblk);

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
