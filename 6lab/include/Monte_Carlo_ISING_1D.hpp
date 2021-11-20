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
#include "../../random/random.hpp"

const int m_props=10;
const int m_spin=50;


class Ising{
public:
	
	//Ising();
	//~Ising();
	
	
	void Input(std::string, int);
	void Reset(void);
	void Accumulate(void);
	void Averages(void);
	void Move(void);
	void Equilibrate(void);
	void ConfFinal(void);
	void Measure(void);
	double Boltzmann(int, int);
	double Error(double, double);
	int Pbc(int);
	
	
	std::string dir_output;
	
	//parameters, observables
	int n_props,iu,ic,im,ix,ig;
	double nbins;
	double walker[m_props];


	int seed[4];
	Random rnd;

	// averages
	double blk_av[m_props],blk_norm,iblk,accepted,attempted;
	double glob_av[m_props],glob_av2[m_props];
	double stima_u,stima_c,stima_m,stima_x,stima_g;
	double err_u,err_c,err_m,err_x,err_g;
	

	//configuration
	double s[m_spin];
	
	// thermodynamical state
	int nspin;
	double beta,temp,J,h;
	
	// simulation
	int nstep, nblk, metro, nequil;
	
	
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
