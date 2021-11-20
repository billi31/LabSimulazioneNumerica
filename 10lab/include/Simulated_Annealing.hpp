/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __ISING__
#define __ISING__

#include <vector>
#include <map>
#include <math.h>
#include "random.hpp"


const int m_props=1000;
const double pi=3.1415927;


struct position {
	double x;
	double y;
};

struct individual {
	std::vector<int> tour;
	double fitness;
	
	//redefine < operator for the struct individual: in this way I can order the tours by confronting L2
	bool operator < (const individual& indiv) const{
		return (fitness < indiv.fitness);
	}
	
};

void one_shift(std::vector<int>& v, int oldIndex, int newIndex);
void removeDuplicates(std::vector<int> &v);

class Simulated_Annealing{
public:
	
	
	void Input(std::string ar);
	void circle_position();
	void random_position();
	void ShuffleIndividuals(void);
	double CostFunction(std::vector<int> indiv);
	void Mutations(void);
	void Mutation1(void);
	void Mutation2(void);
	void Mutation3(void);
	void Mutation4(void);
	void Mutation5(void);
	void Mutation6(void);
	void PrintSolution(void);
	void Move(void);
	void Reset(void);
	
	std::string dir_output;
	
	Random rnd;
	
	std::string area;
	std::string fitness_type;
	
	int ncity;
	double t;
	double beta;
	double initial_t;
	double cooling_rate;
	int step_temperature;
	int MC_step;
	
	int generation;
	int n_individuals;	
	int n_generations;
	int i_indiv;
	int i_gener;
	
	
	//pop contains the population of individuals with each evaluation of the cost function
	std::vector<individual> pop;
	std::vector<individual> new_pop;
	std::map<int, position> cities;
	individual indiv;
	position city;
	
	bool check_ind;
	
	double alpha;
	int accepted, attempted;
	
};



#endif

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
