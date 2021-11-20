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

//Random numbers
#include "random.hpp"

const int m_props=10;
const double pi=3.1415927;


class VMC{
public:
	
	//Ising();
	//~Ising();
	
	
	void Input(void);
	void Reset(void);
	void Accumulate(void);
	void Averages(std::string);
	void Move(void);
	void Measure(void);
	void Optimization(void);
	double SquaredDistribution(double x);
	double Error(double, double);
	double KineticEnergy(void);	
	double PotentialEnergy(void);
	void Instantaneus(void);
	
	std::string dir_output;
	
	//parameters, observables
	int n_props,ie;
	double walker[m_props];


	int seed[4];
	Random rnd;

	// averages
	double blk_av[m_props],blk_norm,accepted,attempted;
	double glob_av[m_props],glob_av2[m_props];
	double stima_e,err_e;
		
	double x, edge;
	double mu, sigma;
	double mu_opt, sigma_opt, e_gs;
	double initial_mu, initial_sigma;
	double delta_mu, delta_sigma;
	int i_mu, i_sigma;
	int grid_mu, grid_sigma;
	
	// simulation
	int nstep, nblk, iblk;
	
	double alpha;
	
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
