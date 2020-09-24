#include"individual.h"
#include"population.h"
#include"random.h"
#include <vector>

#include<iostream>
#include<string>

#include"mpi.h"

using namespace std;

int main(int argc, char* argv[]) {

	//Setting MPI enviorment

	MPI_Init(&argc, &argv);
	int size, rank; 
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	//===========Setting Random Generator ================
	//Parallelize

	Random rnd;
	int seed[4];
	int p[size][2];
	ifstream Primes("Primes");
	if (Primes.is_open())
	{
		for(int i=0; i<size; i++)
		{
		Primes >> p[i][0] >> p[i][1];
		}
	} else cerr << "PROBLEM: Unable to open Primes" << endl;
	Primes.close();

	ifstream input("seed.in");
	string property;
	if (input.is_open()){
		while ( !input.eof() ){
			input >> property;
			if( property == "RANDOMSEED" ){
				input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
				rnd.SetRandom(seed,p[rank][0],p[rank][1]);
			}
		}
	input.close();
	} else cerr << "PROBLEM: Unable to open seed.in" << endl;

	//=======================================================

	//Setting parameters

	int n_cities = 32, n_gen = 2000, n_pop = 200, n_mig = 100;
	double perm = 0.1, perm_group = 0.03, shift = 0.05, invs = 0.01, cross = 0.5; 	//Mutation rates
	double side = 1.;	//System dimensions

	if(rank==0)
	{
		cout << endl << endl;
		cout << "Traveling Salesman Problem (TSP) with a Genetic Algorithm (GA)" << endl << endl;
		cout << "Number of n_cities:\t" << n_cities << endl;
		cout << "Number of generations:\t" << n_gen << endl;
		cout << "Number of possible paths for GA:\t" << n_pop << endl << endl;

		cout << "Mutations rates:" << endl;
		cout << "Permutation:\t" << perm << endl;
		cout << "Group permutation:\t" << perm_group << endl;
		cout << "Shift:\t" << shift << endl;
		cout << "Inversion:\t" << invs << endl << endl;

		cout << "L2 is used as the cost function" << endl << endl;
	}

	//Setting GA

	vector<double> x_sq, y_sq;

	for (int i = 0; i < n_cities; i++)
	{
		//Square
		double x = rnd.Rannyu(-side/2., side/2.);
		double y = rnd.Rannyu(-side/2., side/2.);

		x_sq.push_back(x);
		y_sq.push_back(y);
	}

	if(rank==0)
	{
		cout << "Performing GA for the following conditions:" << endl;
		cout << "Cties inside a square of side:\t" << side << endl << endl;
	}

	population square_pop(perm, perm_group, shift, invs, cross, n_pop, x_sq, y_sq, rnd);

	if(rank==0) cout << "First generation ready, initiating GA ... " << endl << endl;

	//Starting GA

	for (int i = 0; i < n_gen; i++)
	{
		square_pop.next_gen();
		
		//Exchanging individuals
		if (i!=0 &&  i%n_mig==0)
		{
			MPI_Status stat1, stat2, stat3, stat4;
			MPI_Request req1, req2;
			int root = 0;
			int itag1 = 1, itag2 = 2, itag3 = 3, itag4 = 4;
			int target[3];

			if(rank==0)
			{
				target[0] = (int)rnd.Rannyu(1,4);
				do
				{
					target[1] = (int)rnd.Rannyu(1,4);
				} while(target[1]==0 || target[1]==target[0]);

				do
				{
					target[2] = (int)rnd.Rannyu(1,4);
				} while(target[2]==0 || target[2]==target[0] || target[2]==target[1]);
			}
			MPI_Bcast(target, 3, MPI_INT, root, MPI_COMM_WORLD);
			
			vector<int> path_exc = square_pop.get_path(0);
			vector<int> new_path(n_cities);

			int* path1 = new int[n_cities];
			int* path2 = new int[n_cities];

			for(int j=0; j<n_cities; j++)
			{
				path1[j] = path_exc[j];
			}
			MPI_Barrier(MPI_COMM_WORLD);
			
			if(rank==0)
			{
				MPI_Isend(&path1[0], n_cities, MPI_INTEGER, target[0], itag1, MPI_COMM_WORLD, &req1);
				MPI_Recv(&path2[0], n_cities, MPI_INTEGER, target[0], itag2, MPI_COMM_WORLD, &stat2);
			}
			else if(rank==target[0])
			{
				MPI_Send(&path1[0], n_cities, MPI_INTEGER, 0, itag2, MPI_COMM_WORLD);
				MPI_Recv(&path2[0], n_cities, MPI_INTEGER, 0, itag1, MPI_COMM_WORLD, &stat1);
			}
			else if(rank==target[1])
			{
				MPI_Isend(&path1[0], n_cities, MPI_INTEGER, target[2], itag3, MPI_COMM_WORLD, &req2);
				MPI_Recv(&path2[0], n_cities, MPI_INTEGER, target[2], itag4, MPI_COMM_WORLD, &stat4);
			}
			else if(rank==target[2])
			{
				MPI_Send(&path1[0], n_cities, MPI_INTEGER, target[1], itag4, MPI_COMM_WORLD);
				MPI_Recv(&path2[0], n_cities, MPI_INTEGER, target[1], itag3, MPI_COMM_WORLD, &stat3);
			}
			for(int j=0; j<n_cities; j++)
			{
				new_path[j] = path2[j];
			}
			
			square_pop.set_path(0, new_path);
			MPI_Barrier(MPI_COMM_WORLD);

			if(rank==0) cout<<" Generation "<<i<<endl;
		}
		square_pop.data_print("Data/sq_ave_"+to_string(rank)+".out");
	}

	if (rank==0)
	{
		cout << endl << "GA completed, printing best path" << endl << endl;
	}

	square_pop.best_print("Data/square_best_"+to_string(rank)+".out");	

	MPI_Finalize();
	return 0;
}
