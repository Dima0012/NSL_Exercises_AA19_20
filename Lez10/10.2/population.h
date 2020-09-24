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
		population(double perm, double perm_group, double shift, double invert, double cross, int pop, vector<double> x, vector<double> y, Random rnd);

		//destructor
		~population(){};
		
		//functions
		double get_perm() const {return m_perm;}
		double get_perm_group() const {return m_perm_group;}
		double get_shift() const {return m_shift;}
		double get_invert() const {return m_invert;}
		double get_cross() const {return m_cross;}
		int get_ncities() const {return m_ncities;}
		int get_pop() const {return m_pop;}
		
		void Sort() {sort(m_population.begin(), m_population.end());}
		
		void first_gen();
		vector<int> random_path();
		void next_gen();
		void mutation();
		void crossover();
		int select();

		vector<int> get_path(int pos);
		void set_path(int pos, vector<int> path);
		
		//Printing data
		void data_print(string filename);
		void best_print(string filename);
	
	protected:
		double m_perm, m_perm_group, m_shift, m_invert, m_cross; //Mutation rates
		int m_pop, m_ncities, m_gen, wd;		//Population parameters
		vector<individual> m_population;	
		vector<double> m_x, m_y;		//cities coordinates
		Random m_rnd;
};

#endif