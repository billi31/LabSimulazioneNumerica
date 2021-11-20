/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __NVT__
#define __NVT__

//Random numbers
#include "random.hpp"

//pigreco
const double pi=3.1415927;
const int m_part=108;
const int m_props=1000;


class FluidLJ{

public:
	int seed[4];
	Random rnd;
	
	//parameters, observables
	int n_props, iv, iw, igofr;
	double vtail,ptail,bin_size,nbins,sd;
	double walker[m_props];
	
	std::string state, mode;
	
	
	// averages
	double blk_av[m_props],blk_norm,accepted,attempted;
	double glob_av[m_props],glob_av2[m_props];
	double stima_pot,stima_pres,err_pot,err_press,err_gdir;
	
	//configuration
	double x[m_part],y[m_part],z[m_part];
	
	// thermodynamical state
	int npart;
	double beta,temp,vol,rho,box,rcut;
	
	// simulation
	int nstep, nblk, iblk, nequil;
	double delta;
	
	
	
	//functions
	void Input(void);
	void Reset(void);
	void Accumulate(void);
	void Averages(void);
	void Move(void);
	void Equilibrate(void);
	void ConfFinal(void);
	void ConfXYZ(int);
	void Measure(void);
	double Boltzmann(double, double, double, int);
	double Error(double, double);
	double Pbc(double);
	double Acceptance(void);
	void MeasureUP(void);
	void InstantaneusUP(void);
	
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
