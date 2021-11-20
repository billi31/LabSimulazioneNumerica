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
const int nblock=100;
const int m_props=5;


class MolDynamics{
public:
	
	MolDynamics();
	~MolDynamics();

	int nstep, iprint, iblk, nequil, nequil_step;
	
	
	
	//functions
	void Input(const char *input_file, const char *config0);
	void Input(const char *input_file, const char *config0, const char *old0);
	void Rescale(void);
	void Reset(void);
	void Accumulate(void);
	void Averages(void);
	void Move(void);
	void Equilibrate(const char *input_file);
	void ConfOld(void);
	void ConfFinal(void);
	void ConfXYZ(int);
	void Measure(void);
	void PrintMeasure();
	double Force(int, int);
	double Pbc(double);
	double Error(double, double);
	
private:
	
	
	//parameters, observables
	int n_props;
	int iv,ik,ie,it,ip;
	double stima_epot, stima_ekin, stima_etot, stima_temp, stima_pres;
	double err_v,err_k,err_e,err_t,err_p;

	
	double inst_measure[m_props];
	
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
