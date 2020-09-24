#include "population.h"

population::population(double perm, double perm_group, double shift, double invert, double cross, int pop, vector<double> x, vector<double> y, Random rnd){
	//mutation rates
	m_perm = perm;
	m_perm_group = perm_group;
	m_shift = shift;
	m_invert = invert;
	m_cross = cross;
	
	m_rnd = rnd;
	m_pop = pop;
	m_ncities = x.size();
	m_x = x;
	m_y = y;
	m_gen = 0;
	wd = 12;

	first_gen();
	Sort();
}

//Generate first population
void population::first_gen()
{
	for(int i=0; i<m_pop; i++)
	{
		vector<int> path = random_path();
		individual random_ind(m_x, m_y, path);
		m_population.push_back(random_ind);
	}
}

//Create random path of cities
vector<int> population::random_path()
{
	vector<int> path;

	for(int i=0; i<m_ncities; i++)
	{
		path.push_back(i);
	}

	for(int i=0; i<100; i++)	//Cycle to randomize path
	{
		int pos1, pos2, p;

		pos1 = (int)m_rnd.Rannyu(1,m_ncities); 	//Start from 1 so first city always in first position
		pos2 = (int)m_rnd.Rannyu(1,m_ncities);
		p = path[pos1];

		path[pos1] = path[pos2];
		path[pos2] = p;
	}

	return path;
}


//Generate next population
void population::next_gen()
{
	Sort();
	vector<individual> next_pop;

	for(int i=0; i<m_pop; i++)
	{
		int pos = select();
		individual chosen_ind = m_population[pos];
		chosen_ind.set_length();

		next_pop.push_back(chosen_ind);
	}
	m_population = next_pop;

	//Crossover and mutation
	int k = (int)m_rnd.Rannyu(1,m_pop/5);
	for(int i=0; i<k; i++)
	{ 
		double cross = m_rnd.Rannyu();
		if(cross < m_cross)
		{
			crossover();
		}
	}
	mutation();
	m_gen++;
	Sort();
}

//Perform mutations
void population::mutation()
{
	//Permutation
	int k = (int)m_rnd.Rannyu(1,m_pop);
	for(int i=0; i<k; i++)
	{
		double perm = m_rnd.Rannyu();

		if(perm < m_perm)
		{
			int pos1, pos2, target;

			pos1 = (int)m_rnd.Rannyu(1,m_ncities);
			pos2 = (int)m_rnd.Rannyu(1,m_ncities);
			target = (int)m_rnd.Rannyu(0,m_pop);

			if(pos1==pos2) pos2+=1;
			m_population[target].permutation(pos1,pos2);
			m_population[target].set_length();
		}
	}
	
	//Group permutation
	for(int i=0; i<k; i++)
	{
		double perm_group = m_rnd.Rannyu();
		if(perm_group < m_perm_group)
		{
			int start1, start2, length, target;

			start1 = (int)m_rnd.Rannyu(1,m_ncities/2+1);
			start2 = (int)m_rnd.Rannyu(m_ncities/2+1,m_ncities);
			length = (int)m_rnd.Rannyu(1,m_ncities/2);
			target = (int)m_rnd.Rannyu(0,m_pop);
			
			m_population[target].permutation_group(start1,start2,length);
			m_population[target].set_length();
		}
	}

	//Shift
	for(int i=0; i<k; i++)
	{
		double shift = m_rnd.Rannyu();
		if(shift < m_shift)
		{
			int start, length, target;

			start = (int)m_rnd.Rannyu(1,m_ncities-1);
			length = (int)m_rnd.Rannyu(1,m_ncities/2);
			target = (int)m_rnd.Rannyu(0,m_pop);

			m_population[target].shift(start,length);
			m_population[target].set_length();
		}
	}
	
	//Inversion
	for(int i=0; i<k; i++)
	{
		double invert = m_rnd.Rannyu();
		if(invert < m_invert)
		{
			int start, end, target;

			start = (int)m_rnd.Rannyu(1,m_ncities/2);
			end = (int)m_rnd.Rannyu(m_ncities/2,m_ncities);
			target = (int)m_rnd.Rannyu(0,m_pop);

			m_population[target].inversion(start,end);
			m_population[target].set_length();
		}
	}
}

//Perform crossover
void population::crossover(){

	int pos1 = select();
	int pos2 = select();

	vector<int> parent1 = m_population[pos1].get_path();
	vector<int> parent2 = m_population[pos2].get_path();

	int cut = (int)m_rnd.Rannyu(1,m_ncities);

	vector<int> child1, child2;

	for(int i=0; i<cut; i++)
	{
		child1.push_back(parent1[i]);
		child2.push_back(parent2[i]);
	}

	for(int i=cut; i<m_ncities; i++)
	{
		for(int j=0; j<m_ncities; j++)
		{
			if(parent1[i]==parent2[j]) child1.push_back(parent2[j]);
			if(parent2[i]==parent1[j]) child2.push_back(parent1[j]);
		}
	}

	individual individual_1(m_x, m_y, child1);
	individual individual_2(m_x, m_y, child2);

	individual_1.set_length();
	individual_2.set_length();

	m_population[pos1] = individual_1;
	m_population[pos2] = individual_2;
}

//Select individual
int population::select()
{
	Sort();
	double r = m_rnd.Rannyu();
	return int(m_pop*pow(r,2));
}

//Print data for best half
void population::data_print(const char* filename)
{
	Sort();
	double ave = 0., std_dev = 0.;

	ofstream out(filename,ios::app);

	int k = m_pop/2;

	for(int i=0; i<k; i++)
	{
		ave+=m_population[i].get_length();
	}
	ave/=k;

	for(int i=0; i<k; i++)
	{
		std_dev+=pow( m_population[i].get_length()-ave, 2);
	}
	std_dev = sqrt(std_dev/(k-1));

	out<<setw(wd)<<m_gen<<setw(wd)<<ave<<setw(wd)<<std_dev<<endl;
	out.close();
}

//Print best path
void population::best_print(const char* filename, int print_test)
{
	Sort();

	if (print_test==1)
	{
		cout<<"Square Best ";
	}
	else
	{
		cout<<"Circumference Best ";
	}

	m_population[0].print_path();

	vector<int> best = m_population[0].get_path();
	ofstream out(filename,ios::app);

	for(int i=0; i<m_ncities; i++)
	{
		int pos = best[i];
		out<<setw(wd)<<i<<setw(wd)<<m_x[pos]<<setw(wd)<<m_y[pos]<<endl;
	}

	int pos = best[0];
	out<<setw(wd)<<m_ncities+1<<setw(wd)<<m_x[pos]<<setw(wd)<<m_y[pos]<<endl;
	out.close();
}