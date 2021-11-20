/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <set>
#include "../include/Genetic_Algorithm.hpp"

using namespace std;


void Genetic_Algorithm::ReadInput(int rank){
	
	ifstream ReadInput;
	rnd.Initialize(rnd, rank);
	//Read input informations
	ReadInput.open("input.dat");
	
	ReadInput >> ncity;
	area = "square";
	ReadInput >> n_individuals;
	ReadInput >> n_generations;
	ReadInput >> nmigr;
	ReadInput >> fitness_type;
	ReadInput >> pm;
	ReadInput >> pc;
	
	if (rank == 0){
		cout << "The Traveling Salesman Problem                " << endl;
		cout << "Genetic Algorithm simulation                  " << endl << endl;
		cout << "Number of cities = " << ncity << endl;
		cout << "Area of Traveling = " << area << endl;
		cout << "Number of individuals (population) = " << n_individuals << endl;
		cout << "\nNumber of generations = " << n_generations << "\n";
		cout << "\nNumber of generations before each migration = " << nmigr << "\n";
		cout << "\nCost Function = " << fitness_type;
		cout << "\nProbability  of Mutations = " << pm << endl;
		cout << "Probability of Crossover = " << pc << endl;
	}
	
	ReadInput.close();
	
	return;
}
	
void Genetic_Algorithm::PopulationInput(void){
	pop.clear();
	
	//generate first generation:
	for (int i = 0; i < n_individuals; i++) {
		indiv.tour.clear();
		
		for (int j = 1; j <= ncity; j++) {
			indiv.tour.push_back(j);
		}
		
		ShuffleIndividuals();
		indiv.fitness = CostFunction(indiv.tour);
		
		pop.push_back(indiv);
		
	}
	
	//orders pop based on fitness
	sort(pop.begin(), pop.end());
	generation = 1;
	
	return;
}


void Genetic_Algorithm::circle_position(){
	double theta = rnd.Rannyu(0,2*pi);
	city.x = cos(theta);
	city.y = sin(theta);
}

void Genetic_Algorithm::random_position(){
	city.x = rnd.Rannyu(-1,1);
	city.y = rnd.Rannyu(-1,1);
}


void Genetic_Algorithm::ShuffleIndividuals(void){
	for (int i = ncity-1; i > 1; --i)
		swap(indiv.tour[i], indiv.tour[int(rnd.Rannyu(1,i+1))]);
}


double Genetic_Algorithm::CostFunction(vector<int> ind){
	
	if(fitness_type == "L1"){	
		double L1 = 0.;
		map<int, position>::iterator it = cities.find(ind.back());
		
		for (int i=0; i<ind.size(); i++){
			map<int, position>::iterator it2 = cities.find(ind[i]);
			L1 += sqrt(pow(it->second.x - it2->second.x,2) + pow(it->second.y - it2->second.y,2));
			it = it2;
		}
		
		return L1;
	}
	else if (fitness_type == "L2"){
		double L2 = 0.;
		map<int, position>::iterator it = cities.find(ind.back());
		
		for (int i=0; i<ind.size(); i++){
			map<int, position>::iterator it2 = cities.find(ind[i]);
			L2 += pow(it->second.x - it2->second.x,2) + pow(it->second.y - it2->second.y,2);
			it = it2;
		}
		
		return L2;
	}
	else{
		cerr << "Fitness type not valid" << endl;
		exit(1);
	}
}



bool Genetic_Algorithm::CheckIndividuals(void){
	
	bool check = true;
	int i = 0;
	
	
	while (check == true && i<pop.size()){
		
		set<int> s(pop[i].tour.begin(), pop[i].tour.end());
		check = (s.size() == pop[i].tour.size());
		
		i++;
	}
	
	if (check == false)
		cout << "\n\nProblem in individual number " << i << "\nA city is repeated\n\n";
	
	return check;
}

int Genetic_Algorithm::Selection(void){
	
	double p = 1.5;
	double x = rnd.Rannyu();
	double j = int(n_individuals*pow(x,p));
	
	return j;
}

void Genetic_Algorithm::NewPop(void){
	new_pop.clear();
	
	//pair of different individuals for cross over
	for (int i = 0; i < n_individuals; i=i+2) {
		int j = Selection();
		new_pop.push_back(pop[j]);
		int k = Selection();
		while (j == k)
			k = Selection();
		new_pop.push_back(pop[k]);
	}
	
	pop = new_pop;
}


//swap 1 city
void Genetic_Algorithm::Mutation1(void){
	
	int i = int(rnd.Rannyu(1,ncity));
	int k = int(rnd.Rannyu(1,ncity));
	
	while (k == i){
		k = int(rnd.Rannyu(1,ncity));
	}
	
	swap(pop[i_indiv].tour[i], pop[i_indiv].tour[k]);
	
	pop[i_indiv].fitness = CostFunction(pop[i_indiv].tour);

}

//swap 2 contiguous cities
void Genetic_Algorithm::Mutation2(void){
	
	int i = int(rnd.Rannyu(1,ncity));
	
	if (i == ncity - 1){
		int k = int(rnd.Rannyu(2,ncity-2));
		
		swap(pop[i_indiv].tour[i], pop[i_indiv].tour[k]);
		swap(pop[i_indiv].tour[1], pop[i_indiv].tour[k+1]);
	}
	else{
		int k = int(rnd.Rannyu(1,ncity));
		while (k >= i-1 && k <= i+1){
			k = int(rnd.Rannyu(1,ncity));
		}
		
		if (k == ncity - 1){
			swap(pop[i_indiv].tour[i], pop[i_indiv].tour[k]);
			swap(pop[i_indiv].tour[i+1], pop[i_indiv].tour[1]);
		}
		else
			swap_ranges(pop[i_indiv].tour.begin()+i, pop[i_indiv].tour.begin()+i+2, pop[i_indiv].tour.begin()+k);
	}
	
	pop[i_indiv].fitness = CostFunction(pop[i_indiv].tour);
}

//swap 3 contiguous cities
void Genetic_Algorithm::Mutation3(void){
	
	int i = int(rnd.Rannyu(1,ncity));
	
	if (i == ncity - 1){
		int k = int(rnd.Rannyu(3,ncity-3));
		
		swap(pop[i_indiv].tour[i], pop[i_indiv].tour[k]);
		swap(pop[i_indiv].tour[1], pop[i_indiv].tour[k+1]);
		swap(pop[i_indiv].tour[2], pop[i_indiv].tour[k+2]);
	}
	else if (i == ncity - 2){
		int k = int(rnd.Rannyu(2,ncity-4));
		
		swap(pop[i_indiv].tour[i], pop[i_indiv].tour[k]);
		swap(pop[i_indiv].tour[i+1], pop[i_indiv].tour[k+1]);
		swap(pop[i_indiv].tour[1], pop[i_indiv].tour[k+2]);
	}
	
	else{
		int k = int(rnd.Rannyu(1,ncity));
		while (k >= i-2 && k <= i+2){
			k = int(rnd.Rannyu(1,ncity));
		}
		
		if (k == ncity - 1){
			swap(pop[i_indiv].tour[i], pop[i_indiv].tour[k]);
			swap(pop[i_indiv].tour[i+1], pop[i_indiv].tour[1]);
			swap(pop[i_indiv].tour[i+2], pop[i_indiv].tour[2]);
		}
		else if (k == ncity - 2){
			swap(pop[i_indiv].tour[i], pop[i_indiv].tour[k]);
			swap(pop[i_indiv].tour[i+1], pop[i_indiv].tour[k+1]);
			swap(pop[i_indiv].tour[i+2], pop[i_indiv].tour[1]);
		}
		else
			swap_ranges(pop[i_indiv].tour.begin()+i, pop[i_indiv].tour.begin()+i+3, pop[i_indiv].tour.begin()+k);
	}
	
	pop[i_indiv].fitness = CostFunction(pop[i_indiv].tour);
}

void move(vector<int>& v, int oldIndex, int newIndex){
	if (oldIndex > newIndex)
		std::rotate(v.rend() - oldIndex - 1, v.rend() - oldIndex, v.rend() - newIndex);
	else        
		std::rotate(v.begin() + oldIndex, v.begin() + oldIndex + 1, v.begin() + newIndex + 1);
}

//shift 2 contiguous cities by 1
void Genetic_Algorithm::Mutation4(void){
	int i = int(rnd.Rannyu(1,ncity));
	
	if (i == ncity - 1){
		move(pop[i_indiv].tour,i,1);
		move(pop[i_indiv].tour,3,i);
	}
	else if (i == ncity - 2){
		move(pop[i_indiv].tour,i,1);
		swap(pop[i_indiv].tour[i+1], pop[i_indiv].tour[2]);
	}
	
	else{
		move(pop[i_indiv].tour,i,i+2);
		move(pop[i_indiv].tour,i,i+2);
	}
	
	
	pop[i_indiv].fitness = CostFunction(pop[i_indiv].tour);

}

//shift 2 contiguous cities by 2
void Genetic_Algorithm::Mutation5(void){
	int i = int(rnd.Rannyu(1,ncity));
	
	if (i == ncity - 1){
		swap(pop[i_indiv].tour[i],pop[i_indiv].tour[2]);
		swap(pop[i_indiv].tour[3], pop[i_indiv].tour[1]);
	}
	else if (i == ncity - 2){
		swap_ranges(pop[i_indiv].tour.begin()+i, pop[i_indiv].tour.begin()+i+2, pop[i_indiv].tour.begin()+1);
	}
	else if (i == ncity - 3){
		swap(pop[i_indiv].tour[i],pop[i_indiv].tour[i+2]);
		swap(pop[i_indiv].tour[i+1], pop[i_indiv].tour[1]);
	}
	
	else{
		swap_ranges(pop[i_indiv].tour.begin()+i, pop[i_indiv].tour.begin()+i+2, pop[i_indiv].tour.begin()+i+2);
	}
	
	pop[i_indiv].fitness = CostFunction(pop[i_indiv].tour);
	
}

//order inversion of 4 cities
void Genetic_Algorithm::Mutation6(void){
	int i = int(rnd.Rannyu(1,ncity));
	
	if (i == ncity - 1){
		swap(pop[i_indiv].tour[i],pop[i_indiv].tour[3]);
		swap(pop[i_indiv].tour[2], pop[i_indiv].tour[1]);
		
	}
	else if (i == ncity - 2){
		swap(pop[i_indiv].tour[i],pop[i_indiv].tour[2]);
		swap(pop[i_indiv].tour[i+1], pop[i_indiv].tour[1]);
	}
	else if (i == ncity - 3){
		swap(pop[i_indiv].tour[i],pop[i_indiv].tour[1]);
		swap(pop[i_indiv].tour[i+1], pop[i_indiv].tour[i+2]);
		
	}
	
	else{
		vector<int>::iterator it = pop[i_indiv].tour.begin();
		reverse(it+i, it+i+4);
	}
	
	
	pop[i_indiv].fitness = CostFunction(pop[i_indiv].tour);
	
}

void Genetic_Algorithm::Mutations(void){
	if (rnd.Rannyu()<0.1){
		int p = int(rnd.Rannyu(0,6));
		if (p == 0)
			Mutation1();
		if (p == 1)
			Mutation2();
		if (p == 2)
			Mutation3();
		if (p == 3)
			Mutation4();
		if (p == 4)
			Mutation5();
		if (p == 5)
			Mutation6();
	}
	
	return;
}

void removeDuplicates(vector<int> &v){
	auto end = v.end();
	for (auto it = v.begin(); it != end; ++it){
		//indicates the place in which there are elements to remove
		end = remove(it + 1, end, *it);
	}
	
	v.erase(end, v.end());
}

void Genetic_Algorithm::CrossOver(int j, int k){
	if (rnd.Rannyu()<=pc){
		
		vector<int> vec_j = pop[j].tour;
		int i = int(rnd.Rannyu(2,ncity));
		
		//insert in pop[j] the elements of pop[k]
		vector<int>::iterator it = pop[j].tour.insert(pop[j].tour.begin() + i, pop[k].tour.begin()+1, pop[k].tour.end());
		
		removeDuplicates(pop[j].tour);
		
		//insert in pop[k] the elements of vec_j
		it = pop[k].tour.insert(pop[k].tour.begin() + i, vec_j.begin()+1, vec_j.end());
		
		removeDuplicates(pop[k].tour);
		
		
		pop[j].fitness = CostFunction(pop[j].tour);
		pop[k].fitness = CostFunction(pop[k].tour);
	}
}

void Genetic_Algorithm::PrintPopulation(void){
	for(int i=0;i<pop.size();i++){
		cout << "\n\nIndividual " << i << "\n";
		for (auto j: pop[i].tour)
			cout << j << "  ";
	}
}

void Genetic_Algorithm::PrintSolution(int rank){
	ofstream out;
	
	string srank = to_string(rank);
	out.open("output/"+area+"/"+fitness_type+"/"+srank+"/solution.dat");
	for (auto j: pop[0].tour){
		map<int, position>::iterator it = cities.find(j);
		out << j << "  " << it->second.x << "    " << it->second.y << "\n";
	}

	map<int, position>::iterator it2 = cities.find(1);	
	out << 1 << "  " << it2->second.x << "    " << it2->second.y << "\n";		
	out.close();
	
	return;
}

double Genetic_Algorithm::FitnessAve(void){
	double fitness_average = 0;
	
	for (int i=0; i<pop.size()/2.; i++)
		fitness_average += pop[i].fitness;
	
	
	return fitness_average*2./pop.size();
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
