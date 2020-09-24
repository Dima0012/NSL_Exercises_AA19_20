#include "individual.h"

individual::individual(vector<double> x, vector<double> y, vector<int> path){
	m_x = x;
	m_y = y;
	m_path = path;
	m_cities = m_path.size();
	set_length();
}

//Set lenght of the path for the individual
void individual::set_length()
{	
	double length = 0;
	for(int i=0; i<m_cities; i++)
	{
		int a = m_path[i];
		int b = m_path[Pbc(i+1)];
		double dx = m_x[a] - m_x[b];
		double dy = m_y[a] - m_y[b];
		length+=sqrt(dx*dx + dy*dy);
	}
	
	m_length = length;
}

//Print current path for the individual
void individual::print_path() const
{
	cout<<"Path: ";
	for(int i=0; i<m_cities; i++) cout<<m_path[i]<<" ";
	cout<<"Length "<<m_length<<endl;
}


//=========== Mutation functions ==============

//Exchange two cities
void individual::permutation(int pos1, int pos2)
{
	pos1 = Pbc(pos1);
	pos2 = Pbc(pos2);
	int p = m_path[pos1];

	m_path[pos1] = m_path[pos2];
	m_path[pos2] = p;

	check();
}

//Exchange a group of continuos cities
void individual::permutation_group(int start1, int start2, int length)
{
	for(int i=0; i<length; i++)
	{
		permutation(Pbc(start1+i),Pbc(start2+i));
	}

	check();
}

//Shift a group of continuos cities
void individual::shift(int start, int length)
{
	for(int i=0;i<length;i++)
	{
		permutation(Pbc(start+i),Pbc(start+length+i));
	}

	check();
}

//Invert the order of the cities
void individual::inversion(int start, int end)
{
	int length = floor((start-end)/2);
	for(int i=0; i<length; i++)
	{
		permutation(Pbc(start+i),Pbc(end-i));
	}
	check();
}

//===============================================

//Check if a city is visited more than once
void individual::check()
{
	int sum = 0;
	for(int i=0; i<m_cities; i++)
	{
		for(int j=i+1; j<m_cities; j++)
		{
			if(m_path[i]==m_path[j])
			{
				cout<<"City "<<m_path[i]<<" visited more than once"<<endl;
				sum++;
			}
		}
	}
	if(sum!=0) print_path();
}

//Periodic bounday conditions
int individual::Pbc(int i)
{
	if(i>=m_cities) i=i-m_cities;
	return i;
}

//Defining < for sorting
bool operator<(const individual& path1, const individual& path2)
{
	return path1.get_length() < path2.get_length();
}