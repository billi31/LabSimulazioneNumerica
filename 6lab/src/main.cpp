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
#include "../include/Monte_Carlo_ISING_1D.hpp"
#include "boost/program_options.hpp"



namespace po = boost::program_options;



int main(int argc, const char** argv){
	
	
	std::string mode;
	
	// Declare the supported options.
	po::options_description desc("Allowed options");
	desc.add_options()
	("help", "produce help message")
	("mode", po::value<std::string>(), "start from random or old configuration (MUST be start or repeat)")
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
	
	if (vm.count("mode")) {
		if (vm["mode"].as<std::string>() == "start" || vm["mode"].as<std::string>() == "repeat")
			std::cout << "Mode to be used is " << vm["mode"].as<std::string>() << ".\n\n";
		else{
			std::cerr << "Mode not valid\n\n";
			exit(1);
		}
		mode = vm["mode"].as<std::string>();
	}
	else {
		std::cout << desc << "\n";
		return 1;
	}
	
	Ising IS;
	
	for (int i = 0; i<31 ; i++){
		IS.temp = 0.5 + i*0.05;
		IS.beta = 1.0/IS.temp;
		
		IS.Input(mode, i); //Inizialization
		
		std::cout << "Temperature = " << IS.temp << std::endl;
		
		IS.Equilibrate();
		
		for(IS.iblk=1; IS.iblk <= IS.nblk; ++IS.iblk){ //Simulation
			IS.Reset();   //Reset block averages
			
			for(int istep=1; istep <= IS.nstep; ++istep){
				IS.Move();
				IS.Measure();
				IS.Accumulate(); //Update block averages
			}
			
			IS.Averages();   //Print results for current block
		}
		IS.ConfFinal(); //Write final configuration
	}
	
	IS.rnd.SaveSeed();
	
	return 0;
}
