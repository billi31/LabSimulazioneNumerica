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
#include "../include/Monte_Carlo_ISING_1D.hpp"

using namespace std;


void Ising::Input(string mode, int itemp){
	
	ifstream ReadInput;
	cout << "Classic 1D Ising model             " << endl;
	cout << "Monte Carlo simulation             " << endl << endl;
	cout << "Nearest neighbour interaction      " << endl << endl;
	cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
	cout << "The program uses k_B=1 and mu_B=1 units " << endl;
	
	rnd.Initialize(rnd);
	
	//Read input informations
	ReadInput.open("input.dat");
	
	ReadInput >> nspin;
	cout << "Number of spins = " << nspin << endl;
	
	ReadInput >> J;
	cout << "Exchange interaction = " << J << endl;
	
	ReadInput >> h;
	cout << "External field = " << h << endl << endl;
	
	ReadInput >> metro; // if=1 Metropolis else Gibbs
	
	ReadInput >> nblk;
	
	ReadInput >> nstep;
	
	ReadInput >> nequil;
	
	//for the directories
	std::ostringstream streamObj;
	streamObj << std::setprecision(1);
	streamObj << h;
	string str_h = streamObj.str();
	
	if(metro==1){
		cout << "The program perform Metropolis moves" << endl;
		dir_output = "metro/output" + str_h + "/";
	}
	else {
		cout << "The program perform Gibbs moves" << endl;
		dir_output = "gibbs/output" + str_h + "/";
	}
	
	cout << "Number of blocks = " << nblk << endl;
	cout << "Number of steps in one block = " << nstep << endl;
	cout << "Number of equilibration steps = " << nequil << endl << endl;
	ReadInput.close();
	
	attempted = 0;
	accepted = 0;
	
	
	//Prepare arrays for measurements
	iu = 0; //Energy
	ic = 1; //Heat capacity
	im = 2; //Magnetization
	ix = 3; //Magnetic susceptibility
	
	n_props = 4; //Number of observables
	
	if (mode == "start" && itemp == 0){
		//initial configuration
		for (int i=0; i<nspin; ++i){
			if (rnd.Rannyu() >= 0.5)
				s[i] = 1;
			else
				s[i] = -1;
		}
		
		cout << "Print initial configuration to file config.0 " << endl << endl;
		ofstream WriteConf;
		WriteConf.open("config.0");
		for (int i=0; i<nspin; ++i)
			WriteConf << s[i] << endl;
		WriteConf.close();
	}
	else {
		ifstream ReadConf;
		ReadConf.open("config.final");
		
		for (int i=0; i<nspin; ++i)
			ReadConf >> s[i];
	}

	
	//Evaluate energy etc. of the initial configuration
	Measure();
	
	//Print initial values for the potential energy and virial
	cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
}

void Ising::Move(void){
	int o;
	double p, energy_old, energy_new;
	double energy_up, energy_down;
	
	
	for(int i=0; i<nspin; ++i){
		
		//Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
		o = (int)(rnd.Rannyu(0,nspin));
		
		if(metro==1){ //Metropolis
			
			//sto cercando di flippare s[o]
			energy_old = Boltzmann(s[o],o);
			energy_new = Boltzmann(-s[o],o);
			
			alpha = min(1.,exp(-beta*energy_new + beta*energy_old));
			
			if (rnd.Rannyu() <= alpha) {
				s[o] = -s[o];
				accepted++;
			}
			
			attempted++;
			
		}
		
		else{ //Gibbs sampling
			energy_up = Boltzmann(1,o);
			energy_down = Boltzmann(-1,o);
			
			//probability for +1
			p = 1./(1. + exp(-beta*energy_down+beta*energy_up));
			
			if (rnd.Rannyu() <= p)
				s[o] = 1;
			else
				s[o] = -1;
			
		}
	}
}

//valore sm dello spin numero ip
double Ising::Boltzmann(int sm, int ip){
	double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
	
	return ene;
}

void Ising::Measure(){
	//int bin;
	double u = 0.0, m = 0.0;
		
	//cycle over spins of only one Ising
	for (int i=0; i<nspin; ++i){
		//internal energy and heat capacity
		u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
		
		//magnetization and magnetic susceptibility
		m += s[i];
		
	}
	
	//u farà l'energia media, c farà l'energia quadratica media
	//m farà la magnetizzazione media, x la magnetizzazione quadratica media
	walker[iu] = u;
	walker[ic] = u*u;
	walker[im] = m;
	walker[ix] = m*m;
}

void Ising::Equilibrate(void){
	cout << "\nStarting the equilibration...\n\n";
	for (int i=0; i<nequil; i++)
		Move();
	
	return;
}


void Ising::Reset(void){ //Reset block averages
	
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


void Ising::Accumulate(void){ //Update block averages
	
	for(int i=0; i<n_props; ++i)
		blk_av[i] += walker[i];
	
	blk_norm = blk_norm + 1.0;
}


void Ising::Averages(void){ //Print results for current block
	
	ofstream Ene, Heat, Mag, Chi;

	const int wd = 12;
	
	cout << "Block number " << iblk << endl;
	if (metro == 1)
		cout << "Metropolis' acceptance rate " << accepted/attempted << endl << endl;
	
	//for the directories
	std::ostringstream streamObj;
	streamObj << std::setprecision(3);
	streamObj << temp;
	string str_temp = streamObj.str();
	
	string dir = dir_output + "out" + str_temp + "/";
	
	
	Ene.open(dir+"output.ene.0",ios::app);
	stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy
	glob_av[iu]  += stima_u;
	glob_av2[iu] += stima_u*stima_u;
	err_u = Error(glob_av[iu],glob_av2[iu]);
	Ene << setw(wd) << iblk << setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
	Ene.close();
	
	
	Heat.open(dir+"output.heat.0",ios::app);
	stima_c = beta*beta*(blk_av[ic]/blk_norm/(double)nspin-stima_u*stima_u*nspin); //Heat capacity
	glob_av[ic]  += stima_c;
	glob_av2[ic] += stima_c*stima_c;
	err_c = Error(glob_av[ic],glob_av2[ic]);
	Heat << setw(wd) << iblk << setw(wd) << stima_c << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
	Heat.close();
	
	
	Mag.open(dir+"output.mag.0",ios::app);
	if (h == 0)
		stima_m = 0;
	else
		stima_m = blk_av[im]/blk_norm/(double)nspin; //Magnetization
	
	glob_av[im]  += stima_m;
	glob_av2[im] += stima_m*stima_m;
	err_m = Error(glob_av[im],glob_av2[im]);
	Mag << setw(wd) << iblk << setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
	Mag.close();
	
	
	Chi.open(dir+"output.chi.0",ios::app);
	stima_x = beta*(blk_av[ix]/blk_norm/(double)nspin-stima_m*stima_m*nspin); //Magnetic susceptibility
	glob_av[ix]  += stima_x;
	glob_av2[ix] += stima_x*stima_x;
	err_x = Error(glob_av[ix],glob_av2[ix]);
	Chi << setw(wd) << iblk << setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
	Chi.close();
	
	
	cout << "----------------------------" << endl << endl;
}


void Ising::ConfFinal(void){
	ofstream WriteConf;
	const int wd = 12;

	ofstream Ene_f, Heat_f, Mag_f, Chi_f;
	
	Ene_f.open(dir_output + "output.ene.temp.0",ios::app);
	Ene_f << setw(wd) << temp << setw(wd) << glob_av[iu]/(double)nblk << setw(wd) << err_u << endl;
	
	Ene_f.close();

	Heat_f.open(dir_output + "output.heat.temp.0",ios::app);
	Heat_f << setw(wd) << temp << setw(wd) << glob_av[ic]/(double)nblk << setw(wd) << err_c << endl;
	Heat_f.close();

	Mag_f.open(dir_output + "output.mag.temp.0",ios::app);
	Mag_f << setw(wd) << temp << setw(wd) << glob_av[im]/(double)nblk << setw(wd) << err_m << endl;
	Mag_f.close();

	Chi_f.open(dir_output + "output.chi.temp.0",ios::app);
	Chi_f << setw(wd) << temp << setw(wd) << glob_av[ix]/(double)nblk << setw(wd) << err_x << endl;
	
	Chi_f.close();
	
	
	cout << "Print final configuration to file config.final " << endl << endl;
	WriteConf.open("config.final");
	
	for (int i=0; i<nspin; ++i)
		WriteConf << s[i] << endl;
	
	WriteConf.close();
	
	rnd.SaveSeed();
}

double Ising::Error(double sum, double sum2){
	if (iblk==1)
		return 0.0;
	else
		return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

int Ising::Pbc(int i){  //Algorithm for periodic boundary conditions

	if (i >= nspin)
		i = i - nspin;
	else if (i < 0)
		i = i + nspin;
	
	return i;
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
