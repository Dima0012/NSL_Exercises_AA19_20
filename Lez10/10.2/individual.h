#ifndef __individual_h__
#define __individual_h__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
using namespace std;

class individual{
	public:
        //constructors
		individual(){};
		individual(vector<double> x, vector<double>y, vector<int> path);

        //destructor
		~individual(){};
		
		//functions
		void set_length();
		double get_length() const {return m_length;}
		double get_x(int pos) const {return m_x[pos];}
		double get_y(int pos) const {return m_y[pos];}
		int get_size() const {return m_cities;}
		vector<int> get_path() const {return m_path;}
		
		//check functions
		void print_path() const;
		void check();
		
		//mutations
		void permutation(int pos1, int pos2);
		void permutation_group(int start1, int start2, int length);
		void shift(int start, int length);
		void inversion(int start, int end);
	
	protected:
		int m_cities;	//Number of cities
		double m_length;	//Length of path
		vector<int> m_path;	//Path for TSP
		vector<double> m_x, m_y; //Cities coordinates
		int Pbc(int i);	//Periodic boundary conditions
};

bool operator<(const individual& path1, const individual& path2);

#endif