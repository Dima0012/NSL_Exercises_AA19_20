#ifndef __population_h__
#define __population_h__

#include "individual.h"
#include "random.h"
#include <string>

using namespace std;

class population{
	public:
		//constructors
		population(){};
		population(double perm, double perm_group, double shift, double invert, double beta, double cool, int cool_step, vector<double> x, vector<double> y, Random rnd);

		//destructor
		~population(){};
		
		//functions
		double get_perm() const {return m_perm;}
		double get_perm_group() const {return m_perm_group;}
		double get_shift() const {return m_shift;}
		double get_invert() const {return m_invert;}
		int get_ncities() const {return m_ncities;}

		double get_beta() const {return m_beta;}
		double get_cool() const {return m_cool;}
		double get_step() const {return m_cool_step;}
		Random get_rnd() const {return m_rnd;}

		double get_accepted() const {return m_accepted;}
		double get_attempted() const {return m_attempted;}
		
		
		void first_gen();
		vector<int> random_path();
		void next_gen();
		void mutation();
		double Boltzmann(double lenght);
		
		//Printing data
		void data_print(const char* filename);
		void best_print(const char* filename, int print_test);
	
	protected:
		double m_perm, m_perm_group, m_shift, m_invert; //mutation rates
		int m_ncities, m_gen, wd;		//Population parameters	
		vector<double> m_x, m_y;		//cities coordinates
		Random m_rnd;

		//SA parameters
		double m_beta, m_cool; //Inverse T and decreasing unit (MUST BE < 1)
		int m_cool_step;	//Step before decreasing beta
		int m_accepted, m_attempted;
		individual m_path_old, m_path_new;
};

#endif