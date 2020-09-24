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
const int m_props=4;
int n_props;
int iv,ik,it,ie;
double stima_pot, stima_kin, stima_etot, stima_temp;

// averages
double acc,att;

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;

// simulation
int nstep = 1E6, nblocks=1000, nmeasure=100;  //Parameters for arrays
int iprint, seed;
double delta;
bool stab;

// arrays for averages
double* ave_temp = new double[nstep/nmeasure];
double* ave_epot = new double[nstep/nmeasure];
double* ave_ekin = new double[nstep/nmeasure];
double* ave_etot = new double[nstep/nmeasure];

double* mean = new double[nblocks];
double* error = new double[nblocks];

//functions
void Input(void);
void Stabilize(void);
void Move(void);
void ConfFinal(void);
void ConfOldFinal(void);
void ConfLast(void);
void ConfXYZ(int);
void Measure(void);
double Force(int, int);
double Pbc(double);

void MC_Ave_Err_Block(int n_throws, int n_blocks, double* mean, double* error, double* data);
void DataPrint(double* data, double* error,  int n, std::string str);
void Blocking(void);
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
