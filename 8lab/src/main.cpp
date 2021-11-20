//
//  main.cpp
//  6lab
//
//  Created by Administrator on 30/08/21.
//

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cmath>        // rint, pow
#include "../include/VMC.hpp"
#include "boost/program_options.hpp"



namespace po = boost::program_options;

int main(int argc, const char** argv){
	
	std::string mode;
	
	// Declare the supported options.
	po::options_description desc("Allowed options");
	desc.add_options()
	("help", "produce help message")
	("insta", "print 5e3 values of instantaneus energy and psi^2 distribution for mu and sigma in input.dat")
	("fixed", "evaluate energy blocking average for mu and sigma in input.dat")
	("grid", "search for mu and sigma that minimize the energy, starting from values in input.dat: perform two different grids")
	;
	
	po::variables_map vm;
	//takes the options from the command line
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);
	
	if (vm.count("help")) {
		std::cout << desc << "\n";
		return 1;
	}
	
	
	if (argc == 1) {
		std::cout << desc << "\n";
		return 1;
	}
	
	VMC Vmc;
	Vmc.Input(); //Inizialization
	
	
	//instantaneus for autocorrelation and for histogram
	if (vm.count("insta")) {
		int ntrials = 5e5;
		std::cout << "\n\nEvaluation of |psi|^2 \n";
		std::cout << "Mu = " << Vmc.mu << "\tSigma = " << Vmc.sigma << "\n";
		std::cout << "Number of Evaluations = " << ntrials << "\n";
		
		for(int istep=1; istep <= ntrials; ++istep){
			Vmc.Move();
			Vmc.Measure();
			Vmc.Instantaneus();		
		}
		
		std::cout << "\nAcceptance = " << Vmc.accepted/Vmc.attempted << "\n";
	}
	
	if (vm.count("fixed")) {
		
		for(Vmc.iblk=1; Vmc.iblk <= Vmc.nblk; ++Vmc.iblk){ //Simulation
			Vmc.Reset();   //Reset block averages
			
			for(int istep=1; istep <= Vmc.nstep; ++istep){
				Vmc.Move();
				Vmc.Measure();
				Vmc.Accumulate(); //Update block averages
			}
			Vmc.Averages("");   //Print results for current block
		}
	}
	
	
	
	if (vm.count("grid")) {
		
		//optimization of mu and sigma
		Vmc.Optimization();
		
		std::cout << "\n\nOptimal mu = " << Vmc.mu_opt << "\n";
		std::cout << "Optimal sigma = " << Vmc.sigma_opt << "\n";
		std::cout << "Ground State Energy = " << Vmc.e_gs << "\n\n";
		
		//second optimization: with block average
		for (Vmc.i_mu = 0; Vmc.i_mu < 6; Vmc.i_mu ++){
			Vmc.mu = Vmc.mu_opt + (Vmc.i_mu-3)*Vmc.delta_mu/10.;
			for (Vmc.i_sigma = 0; Vmc.i_sigma < 6; Vmc.i_sigma ++){
				Vmc.sigma = Vmc.sigma_opt + (Vmc.i_sigma-3)*Vmc.delta_sigma/10.;
				
				for(Vmc.iblk=1; Vmc.iblk <= Vmc.nblk; ++Vmc.iblk){ //Simulation
					Vmc.Reset();   //Reset block averages
					
					for(int istep=1; istep <= Vmc.nstep; ++istep){
						Vmc.Move();
						Vmc.Measure();
						Vmc.Accumulate(); //Update block averages
					}
					Vmc.Averages("grid");   //Print results for current block
				}
			}
		}
	}
	
	Vmc.rnd.SaveSeed();
	
	return 0;
}
