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
#include "../include/Simulated_Annealing.hpp"

using namespace std;


void Simulated_Annealing::Input(string ar){
	
	ifstream ReadInput;
	cout << "The Traveling Salesman Problem                " << endl;
	cout << "Simulated Annealing simulation                  " << endl << endl;
	//cout << "Variational Monte Carlo simulation            " << endl << endl;
	//cout << "Nearest neighbour interaction                 " << endl << endl;
	//cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
	cout << "The program uses k_B=1 and mu_B=1 units       " << endl;
	
	rnd.Initialize(rnd);
	
	//Read input informations
	ReadInput.open("input.dat");
	
	ReadInput >> ncity;
	cout << "Number of cities = " << ncity << endl;
	
	area = ar;
	
	ofstream out;
	out.open("output/"+area+"/"+fitness_type+"/cities.dat");
	
	// Cities on the map
	for (int i=1; i<=ncity; i++){
		if (area == "circle")
			circle_position();
		else if (area == "square")
			random_position();
		else
			cout << "\n\nArea not valid\n\n";
		cities.insert(pair<int,position>(i,city));
		out << i << "    " << city.x << "    " << city.y << endl;
		
	}
	
	ReadInput >> initial_t;
	cout << "Initial Temperature = " << initial_t << endl;
	
	ReadInput >> cooling_rate;
	cout << "Cooling Rate = " << cooling_rate << endl;
	
	ReadInput >> step_temperature;
	cout << "Temperature Step = " << step_temperature << endl;
	
	ReadInput >> MC_step;
	cout << "Monte Carlo Step = " << MC_step << endl;
	
	//generate first individual:
	indiv.tour.clear();
	
	for (int j = 1; j <= ncity; j++)
		indiv.tour.push_back(j);
	
	ShuffleIndividuals();
	indiv.fitness = CostFunction(indiv.tour);
	
	cout << endl << endl;
	ReadInput.close();
	
	accepted = 0;
	attempted = 0;
	
}

void Simulated_Annealing::circle_position(){
	double theta = rnd.Rannyu(0,2*pi);
	city.x = cos(theta);
	city.y = sin(theta);
}

void Simulated_Annealing::random_position(){
	city.x = rnd.Rannyu(-1,1);
	city.y = rnd.Rannyu(-1,1);
}

void Simulated_Annealing::ShuffleIndividuals(void){
	for (int i = ncity-1; i > 1; --i)
		swap(indiv.tour[i], indiv.tour[int(rnd.Rannyu(1,i+1))]);
}


double Simulated_Annealing::CostFunction(vector<int> ind){
	
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

//swap 1 city
void Simulated_Annealing::Mutation1(void){

	int i = int(rnd.Rannyu(1,ncity));
	int k = int(rnd.Rannyu(1,ncity));
	
	while (k == i){
		k = int(rnd.Rannyu(1,ncity));
	}
	
	swap(indiv.tour[i], indiv.tour[k]);
	
	indiv.fitness = CostFunction(indiv.tour);
	
}

//swap 2 contiguous cities
void Simulated_Annealing::Mutation2(void){
	
	int i = int(rnd.Rannyu(1,ncity));
	
	if (i == ncity - 1){
		int k = int(rnd.Rannyu(2,ncity-2));
		
		swap(indiv.tour[i], indiv.tour[k]);
		swap(indiv.tour[1], indiv.tour[k+1]);
	}
	else{
		int k = int(rnd.Rannyu(1,ncity));
		while (k >= i-1 && k <= i+1){
			k = int(rnd.Rannyu(1,ncity));
		}
		
		if (k == ncity - 1){
			swap(indiv.tour[i], indiv.tour[k]);
			swap(indiv.tour[i+1], indiv.tour[1]);
		}
		else
			swap_ranges(indiv.tour.begin()+i, indiv.tour.begin()+i+2, indiv.tour.begin()+k);
	}
	
	indiv.fitness = CostFunction(indiv.tour);

}

//swap 3 contiguous cities
void Simulated_Annealing::Mutation3(void){
		
	int i = int(rnd.Rannyu(1,ncity));
	
	if (i == ncity - 1){
		int k = int(rnd.Rannyu(3,ncity-3));
		
		swap(indiv.tour[i], indiv.tour[k]);
		swap(indiv.tour[1], indiv.tour[k+1]);
		swap(indiv.tour[2], indiv.tour[k+2]);
	}
	else if (i == ncity - 2){
		int k = int(rnd.Rannyu(2,ncity-4));
		
		swap(indiv.tour[i], indiv.tour[k]);
		swap(indiv.tour[i+1], indiv.tour[k+1]);
		swap(indiv.tour[1], indiv.tour[k+2]);
	}
	
	else{
		int k = int(rnd.Rannyu(1,ncity));
		while (k >= i-2 && k <= i+2){
			k = int(rnd.Rannyu(1,ncity));
		}
		
		if (k == ncity - 1){
			swap(indiv.tour[i], indiv.tour[k]);
			swap(indiv.tour[i+1], indiv.tour[1]);
			swap(indiv.tour[i+2], indiv.tour[2]);
		}
		else if (k == ncity - 2){
			swap(indiv.tour[i], indiv.tour[k]);
			swap(indiv.tour[i+1], indiv.tour[k+1]);
			swap(indiv.tour[i+2], indiv.tour[1]);
		}
		else
			swap_ranges(indiv.tour.begin()+i, indiv.tour.begin()+i+3, indiv.tour.begin()+k);
	}
	
	indiv.fitness = CostFunction(indiv.tour);
	
}

void move(vector<int>& v, int oldIndex, int newIndex){
	if (oldIndex > newIndex)
		std::rotate(v.rend() - oldIndex - 1, v.rend() - oldIndex, v.rend() - newIndex);
	else        
		std::rotate(v.begin() + oldIndex, v.begin() + oldIndex + 1, v.begin() + newIndex + 1);
}

//shift 2 contiguous cities by 1
void Simulated_Annealing::Mutation4(void){
		
	int i = int(rnd.Rannyu(1,ncity));
	
	if (i == ncity - 1){
		move(indiv.tour,i,1);
		move(indiv.tour,3,i);
	}
	else if (i == ncity - 2){
		move(indiv.tour,i,1);
		swap(indiv.tour[i+1], indiv.tour[2]);
	}
	
	else{
		move(indiv.tour,i,i+2);
		move(indiv.tour,i,i+2);
	}
	
	
	indiv.fitness = CostFunction(indiv.tour);
	
}

//shift 2 contiguous cities by 2
void Simulated_Annealing::Mutation5(void){
		
	int i = int(rnd.Rannyu(1,ncity));
	
	if (i == ncity - 1){
		swap(indiv.tour[i],indiv.tour[2]);
		swap(indiv.tour[3], indiv.tour[1]);
	}
	else if (i == ncity - 2){
		swap_ranges(indiv.tour.begin()+i, indiv.tour.begin()+i+2, indiv.tour.begin()+1);
	}
	else if (i == ncity - 3){
		swap(indiv.tour[i],indiv.tour[i+2]);
		swap(indiv.tour[i+1], indiv.tour[1]);
	}
	
	else{
		swap_ranges(indiv.tour.begin()+i, indiv.tour.begin()+i+2, indiv.tour.begin()+i+2);
	}
	
	indiv.fitness = CostFunction(indiv.tour);
	
}

//order inversion of 4 cities
void Simulated_Annealing::Mutation6(void){

	int i = int(rnd.Rannyu(1,ncity));
	
	if (i == ncity - 1){
		swap(indiv.tour[i],indiv.tour[3]);
		swap(indiv.tour[2], indiv.tour[1]);
		
	}
	else if (i == ncity - 2){
		swap(indiv.tour[i],indiv.tour[2]);
		swap(indiv.tour[i+1], indiv.tour[1]);
	}
	else if (i == ncity - 3){
		swap(indiv.tour[i],indiv.tour[1]);
		swap(indiv.tour[i+1], indiv.tour[i+2]);
		
	}
	
	else{
		vector<int>::iterator it = indiv.tour.begin();
		reverse(it+i, it+i+4);
	}
	
	
	indiv.fitness = CostFunction(indiv.tour);
	
}

void Simulated_Annealing::PrintSolution(void){
	ofstream out;
	out.open("output/"+area+"/"+fitness_type+"/solution.dat");
	for (auto j: indiv.tour){
		map<int, position>::iterator it = cities.find(j);
		out << j << "  " << it->second.x << "    " << it->second.y << "\n";
	}

	map<int, position>::iterator it2 = cities.find(1);	
	out << 1 << "  " << it2->second.x << "    " << it2->second.y << "\n";		
	out.close();
	
	return;
}

void Simulated_Annealing::Mutations(void){
	
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
	
	return;
}

void Simulated_Annealing::Move(void){
	individual old_indiv;
	
	old_indiv = indiv;
	Mutations();
	
	alpha = min(1.,exp(beta*(old_indiv.fitness-indiv.fitness)));
	
	if (rnd.Rannyu() <= alpha)
		accepted++;
	else
		indiv = old_indiv;
	
	attempted++;
}

void Simulated_Annealing::Reset(void){ //Reset acceptance
	attempted = 0;
	accepted = 0;
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
