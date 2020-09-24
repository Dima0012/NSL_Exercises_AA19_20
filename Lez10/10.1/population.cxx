#include "population.h"

population::population(double perm, double perm_group, double shift, double invert, double beta, double cool, int cool_step, vector<double> x, vector<double> y, Random rnd){
	//mutation rates
	m_perm = perm;
	m_perm_group = perm_group;
	m_shift = shift;
	m_invert = invert;

	m_beta = beta;
	m_cool_step = cool_step;
	m_cool = cool;
	m_accepted = 0;
	m_attempted = 0;
	
	m_rnd = rnd;
	m_ncities = x.size();
	m_x = x;
	m_y = y;
	m_gen = 0;
	wd = 12;

	first_gen();
}

//Generate first population
void population::first_gen()
{
	vector<int> path = random_path();
	individual random_ind(m_x, m_y, path);

	m_path_old = random_ind;
	m_path_new = random_ind;

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
	m_path_new = m_path_old;
	mutation();

	//Metropolis sampling
	double p_old = Boltzmann(m_path_old.get_length());
	double p_new = Boltzmann(m_path_new.get_length());

	double p = p_new/p_old;

	if (m_rnd.Rannyu() < p )
	{
		m_path_old = m_path_new;
		m_accepted++;
	}
	m_attempted++;
	m_gen++;

	//Reduce beta every m_cool_step
	if(m_gen%m_cool_step == 0)
	{
		m_beta/=m_cool;
	}
}

//Perform mutations
void population::mutation()
{

	//Permutation	
	double perm = m_rnd.Rannyu();
	if(perm < m_perm)
	{
		int pos1, pos2;

		pos1 = (int)m_rnd.Rannyu(1,m_ncities);
		pos2 = (int)m_rnd.Rannyu(1,m_ncities);

		if(pos1==pos2) pos2+=1;
		m_path_new.permutation(pos1,pos2);
		m_path_new.set_length();
	}

	
	//Group permutation
	double perm_group = m_rnd.Rannyu();
	if(perm_group < m_perm_group)
	{
		int start1, start2, length;

		start1 = (int)m_rnd.Rannyu(1,m_ncities/2+1);
		start2 = (int)m_rnd.Rannyu(m_ncities/2+1,m_ncities);
		length = (int)m_rnd.Rannyu(1,m_ncities/2);
		
		m_path_new.permutation_group(start1,start2,length);
		m_path_new.set_length();
	}


	//Shift
	double shift = m_rnd.Rannyu();
	if(shift < m_shift)
	{
		int start, length;

		start = (int)m_rnd.Rannyu(1,m_ncities-1);
		length = (int)m_rnd.Rannyu(1,m_ncities/2);
		
		m_path_new.shift(start,length);
		m_path_new.set_length();
	}
	
	//Inversion
	double invert = m_rnd.Rannyu();
	if(invert < m_invert)
	{
		int start, end;

		start = (int)m_rnd.Rannyu(1,m_ncities/2);
		end = (int)m_rnd.Rannyu(m_ncities/2,m_ncities);

		m_path_new.inversion(start,end);
		m_path_new.set_length();
	}
	
}

//Compute Boltzamnn weight
double population::Boltzmann(double lenght)
{
	return exp(-m_beta*lenght);
}

//Print data
void population::data_print(const char* filename)
{
	ofstream out(filename,ios::app);
	out<<setw(wd)<<m_gen<<setw(wd)<<m_beta<<setw(wd)<<m_path_old.get_length()<<setw(wd)<<endl;
	out.close();
}

//Print best path
void population::best_print(const char* filename, int print_test)
{
	cout << "Metropolis acceptance rate: " << double(m_accepted)/m_attempted*100 << " %" << endl;

	if (print_test==1)
	{
		cout<<"Square Best ";
	}
	else
	{
		cout<<"Circumference Best ";
	}

	m_path_old.print_path();

	ofstream out(filename,ios::app);
	vector<int> best = m_path_old.get_path();

	for(int i=0; i<m_ncities; i++)
	{
		int pos = best[i];
		out<<setw(wd)<<i<<setw(wd)<<m_x[pos]<<setw(wd)<<m_y[pos]<<endl;
	}

	int pos = best[0];
	out<<setw(wd)<<m_ncities+1<<setw(wd)<<m_x[pos]<<setw(wd)<<m_y[pos]<<endl;
	out.close();

	cout << endl << endl;
}