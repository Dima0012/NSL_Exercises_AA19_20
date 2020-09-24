#include"individual.h"
#include"population.h"
#include"random.h"
#include <vector>

#include<iostream>

using namespace std;

int main(int argc, char* argv[]) {

	//===========Setting Random Generator ================

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

	//Setting parameters

	int n_cities = 32, n_gen = 2000, n_pop = 200;
	double perm = 0.1, perm_group = 0.03, shift = 0.05, invs = 0.01, cross = 0.5; 	//Mutation rates
	double radius = 1., side = 1.;	//System dimensions

	cout << endl << endl;
	cout << "Traveling Salesman Problem (TSP) with a Genetic Algorithm (GA)" << endl << endl;
	cout << "Number of cities:\t" << n_cities << endl;
	cout << "Number of generations:\t" << n_gen << endl;
	cout << "Number of possible paths for GA:\t" << n_pop << endl << endl;

	cout << "Mutations rates:" << endl;
	cout << "Permutation:\t" << perm << endl;
	cout << "Group permutation:\t" << perm_group << endl;
	cout << "Shift:\t" << shift << endl;
	cout << "Inversion:\t" << invs << endl << endl;

	cout << "L1 is used as the cost function" << endl << endl;


	//Setting GA

	vector<double> x_sq, y_sq, x_c, y_c;

	for (int i = 0; i < n_cities; i++)
	{
		//Circumference
		double theta = rnd.Rannyu(0., 2.*M_PI);
		double x = radius*cos(theta);
		double y = radius*sin(theta);

		x_c.push_back(x);
		y_c.push_back(y);

		//Square
		x = rnd.Rannyu(-side/2., side/2.);
		y = rnd.Rannyu(-side/2., side/2.);

		x_sq.push_back(x);
		y_sq.push_back(y);
	}

	cout << "Performing GA for the following conditions:" << endl;
	cout << "Cities on a circumference of radius:\t" << radius << endl;
	cout << "Cties inside a square of side:\t" << side << endl << endl;

	population square_pop(perm, perm_group, shift, invs, cross, n_pop, x_sq, y_sq, rnd);
	population circle_pop(perm, perm_group, shift, invs, cross, n_pop, x_c, y_c, rnd);

	cout << "First generation ready, initiating GA ... " << endl << endl;

	//Starting GA

	for (int i = 0; i < n_gen; i++)
	{
		square_pop.next_gen();
		square_pop.data_print("Data/square_ave.out");

		circle_pop.next_gen();
		circle_pop.data_print("Data/circle_ave.out");

		if (i%100==0)
		{
			cout << "Generation:\t" << i << endl;
		}		
	}

	cout << endl << "GA completed, printing best path:" << endl << endl;

	square_pop.best_print("Data/square_best.out", 1);
	circle_pop.best_print("Data/circle_best.out", 0);	

	return 0;
}
