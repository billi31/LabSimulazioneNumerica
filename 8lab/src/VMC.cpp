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
#include "../include/VMC.hpp"

using namespace std;


void VMC::Input(void){
	
	ifstream ReadInput;
	cout << "1D Single Quantum Particle                    " << endl;
	cout << "Variational Monte Carlo simulation            " << endl << endl;
	//cout << "Nearest neighbour interaction                 " << endl << endl;
	//cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
	//cout << "The program uses k_B=1 and mu_B=1 units       " << endl;
	
	rnd.Initialize(rnd);
	
	//Read input informations
	ReadInput.open("input.dat");
	
	ReadInput >> mu;
	cout << "Mu = " << mu << endl;
	
	ReadInput >> sigma;
	cout << "Sigma = " << sigma << endl;
	
	ReadInput >> edge;
	cout << "Edge = " << edge << endl << endl;

	ReadInput >> nblk;	
	cout << "Number of blocks = " << nblk << endl;

	ReadInput >> nstep;
	cout << "Number of steps in one block = " << nstep << endl << endl;

	ReadInput >> initial_mu;	
	cout << "Initial mu for grid search = " << initial_mu << endl;
	
	ReadInput >> initial_sigma;
	cout << "Initial sigma for grid search = " << initial_sigma << endl << endl;
	
	ReadInput >> delta_mu;	
	cout << "Mu step for grid search = " << delta_mu << endl;
	
	ReadInput >> delta_sigma;
	cout << "Sigma step for grid search = " << delta_sigma << endl << endl;
	
	ReadInput >> grid_mu;	
	cout << "Number of evaluation of mu for grid search = " << grid_mu << endl;
	
	ReadInput >> grid_sigma;
	cout << "Number of evaluation of sigma for grid search = " << grid_sigma << endl;
	cout << endl << endl;
	ReadInput.close();
	
	x = 0.;
	
	//Prepare arrays for measurements
	ie = 0; //Total energy
	n_props = 1; //Number of observables
	
}

void VMC::Move(void){
	double xnew;
	double distribution2_old, distribution2_new;
	
	xnew = x + edge * rnd.Rannyu(-0.5,0.5);
	
	distribution2_old = SquaredDistribution(x);
	distribution2_new = SquaredDistribution(xnew);
	
	alpha = min(1.,distribution2_new/distribution2_old);
	
	if (rnd.Rannyu() <= alpha) {
		x = xnew;
		accepted++;
	}
	
	attempted++;
	
}

double VMC::SquaredDistribution(double x){
	double exp_plus, exp_minus;
	
	exp_plus = exp(-pow(x+mu,2)/(2*pow(sigma,2)));
	exp_minus = exp(-pow(x-mu,2)/(2*pow(sigma,2)));
	
	double psi2 = pow(exp_plus + exp_minus,2);
	return psi2;
}

void VMC::Measure(){
	double t = 0.0, u = 0.0;
	
	//Potential and kinetic energy
	u = PotentialEnergy();
	t = KineticEnergy();
	
	walker[ie] = u + t;
	
}

double VMC::KineticEnergy(void){
	double exp_plus, exp_minus;
	double t;
	
	exp_plus = exp(-pow(x+mu,2)/(2.*pow(sigma,2)));
	exp_minus = exp(-pow(x-mu,2)/(2.*pow(sigma,2)));
	
	t = (1./(2.*pow(sigma,2))) * (1. - (pow(x,2) + pow(mu,2))/pow(sigma,2) - (2*x*mu/pow(sigma,2))*(exp_plus - exp_minus)/(exp_minus + exp_plus));
	
	return t;
}
double VMC::PotentialEnergy(void){
	double u;
	
	u = pow(x,4) - pow(x,2)*2.5;
	
	return u;
}


void VMC::Instantaneus(void){
	ofstream Out;
	int wd = 14;
	
	stima_e = walker[ie]; //Total energy
	
	Out.open("output/instantaneus.psi2_etot.0",ios::app);
	
	//Position x and instantaneus total energy
	Out << setw(wd) << x << setw(wd) << stima_e << endl;
	
	Out.close();
}


void VMC::Reset(void){ //Reset block averages
	
	if(iblk == 1){
		for(int i=0; i<n_props; ++i){
			glob_av[i] = 0;
			glob_av2[i] = 0;
		}
	}
	
	for(int i=0; i<n_props; ++i)
		blk_av[i] = 0;
	
	blk_norm = 0;
	attempted = 0;
	accepted = 0;
}


void VMC::Accumulate(void){ //Update block averages
	
	for(int i=0; i<n_props; ++i)
		blk_av[i] += walker[i];
	
	blk_norm = blk_norm + 1.0;
}


void VMC::Averages(string mode){ //Print results for current block
	
	ofstream Ene;
	
	const int wd = 14;
	
	if (mode == "grid"){
		if (iblk == nblk){
			
			cout << "Mu = " << mu << "\tSigma = " << sigma << endl;
			cout << "Metropolis' acceptance rate " << accepted/attempted << endl << endl;
			
			cout << "----------------------------" << endl << endl;
		}
		Ene.open("output/output.etot." + to_string(i_mu) + to_string(i_sigma) + ".0",ios::app);
	}
	
	else{
		cout << "Block number " << iblk << endl;
		
		cout << "Metropolis' acceptance rate " << accepted/attempted << endl << endl;
		
		cout << "----------------------------" << endl << endl;

		Ene.open("output/output.etot.0",ios::app);
	}
	
	stima_e = blk_av[ie]/blk_norm; //Energy
	glob_av[ie]  += stima_e;
	glob_av2[ie] += stima_e*stima_e;
	err_e = Error(glob_av[ie],glob_av2[ie]);
	Ene << setw(wd) << iblk << setw(wd) << stima_e << setw(wd) << glob_av[ie]/(double)iblk << setw(wd) << err_e << endl;
	Ene.close();
	
	
}


void VMC::Optimization(void){
	int nmove_opt = 1e5;
	e_gs = 100.;
	double e;
	
	int wd = 14;
	ofstream Etot;
	
	Etot.open("output/grid.etot.0",ios::app);
	
	for (int i=0; i<grid_mu; i++){
		mu = initial_mu + i*delta_mu; 
		for (int j=0; j<grid_sigma; j++){
			sigma = initial_sigma + j*delta_sigma;
			
			double t = 0.;
			double u = 0.;
			e = 0.;
			
			for (int k=0; k<nmove_opt; k++){
				Move();
				t += KineticEnergy();
				u += PotentialEnergy();
			}
			e = t + u;
			e /= nmove_opt;
			
			cout << "Mu = " << mu << endl;
			cout << "Sigma = " << sigma << endl;
			cout << "Energy = " << e << endl << endl;
			
			cout << "----------------------------" << endl << endl;
			
			
			//Potential energy per particle
			Etot << setw(wd) << mu << setw(wd) << sigma << setw(wd) << e << endl;
			
			
			
			if (e < e_gs){
				e_gs = e;
				mu_opt = mu;
				sigma_opt = sigma;
			}
			
		}
	}
	Etot.close();
	
	return;
}

double VMC::Error(double sum, double sum2){
	
	if (iblk==1)
		return 0.0;
	else
		return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
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
