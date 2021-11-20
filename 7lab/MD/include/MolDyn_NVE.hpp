/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef MolDyn_NVE_hpp
#define MolDyn_NVE_hpp


//#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <iomanip>
#include "random.hpp"

const int m_part=108;
const int m_props=1000;
const double pi=3.1415927;



class MolDynamics{
public:
	
	MolDynamics();
	~MolDynamics();

	int nstep, iprint, iblk, nequil, nequil_step;
	
	
	
	//functions
	void Input(void);
	void Rescale(void);
	void Move(void);	
	void Accumulate(void);
	void Averages(void);
	void Equilibrate(void);
	void ConfOld(void);
	void ConfFinal(void);
	void ConfXYZ(int);
	void Measure(void);
	void Reset(void);
	void PrintMeasure(int);
	double Force(int, int);
	double Pbc(double);
	double Error(double, double);
	
	
	
	std::string state, mode;
	int nblk;
	Random rnd;
	
private:
	
	
	//parameters, observables
	int n_props;
	int iv, ik, it, ie, ip, igofr;
	double stima_pot, stima_kin, stima_etot, stima_temp, stima_pres;
	double err_pot, err_kin, err_etot, err_temp, err_pres, err_gdir;
	double walker[m_props];
	double bin_size, nbins;
	

	
	double blk_av[m_props],blk_norm;
	double glob_av[m_props],glob_av2[m_props];
	
	
	//configuration
	double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
	double vx[m_part],vy[m_part],vz[m_part];
	
	// thermodynamical state
	int npart;
	double temp,vol,rho,box,rcut;
	
	// simulation
	double delta;

};





#endif /* MolDyn_NVE_hpp */
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
