//
//  main.cpp
//  7lab
//
//  Created by Administrator on 03/09/21.
//

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cmath>        // rint, pow
#include "../include/Monte_Carlo_NVT.hpp"
#include "boost/program_options.hpp"


namespace po = boost::program_options;


int main(int argc, const char** argv){
	
	FluidLJ LJ;
	
	
	// Declare the supported options.
	po::options_description desc("Allowed options");
	desc.add_options()
	("help", "produce help message")
	("state", po::value<std::string>(), "state of fluid (MUST be gas, liquid or solid)")
	("insta", "evaluate 5e5 instantaneus values of potential energy and pressure")
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
	
	if (vm.count("state")) {
		if (vm["state"].as<std::string>() == "gas" || vm["state"].as<std::string>() == "liquid" || vm["state"].as<std::string>() == "solid")
			std::cout << "State of the fluid is " << vm["state"].as<std::string>() << ".\n\n";
		else{
			std::cerr << "State not valid\n\n";
			exit(1);
		}
		LJ.state = vm["state"].as<std::string>();
	}
	else {
		std::cout << desc << "\n";
		return 1;
	}	
	
	
	LJ.Input();			//Inizialization
	LJ.Equilibrate();	//Equilibrate the system
	LJ.Reset();			//Reset acceptance
	
	if (vm.count("insta")) {
		int ntrials = 5e5;
		std::cout << "Starting the simulation...\nMeasuring instantaneus values of Potential Energy per particle and Pressure\n\n";
		for(int istep=1; istep <= ntrials; ++istep){
			LJ.Move();
			LJ.MeasureUP();
			LJ.InstantaneusUP();		
		}
		LJ.ConfFinal(); //Write final configuration
		
		std::cout << "Acceptance of " << LJ.Acceptance() << std::endl << std::endl;
	}
	
	
	//int nconf = 1;
	
	for(LJ.iblk=1; LJ.iblk <= LJ.nblk; ++LJ.iblk){ //Simulation
		LJ.Reset();   //Reset block averages and acceptance
		
		for(int istep=1; istep <= LJ.nstep; ++istep){
			LJ.Move();
			LJ.Measure();
			LJ.Accumulate(); //Update block averages
			
			/*if(istep%10 == 0){
				LJ.ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
				nconf += 1;
			}*/
		}
		LJ.Averages();   //Print results for current block
	}
	
	LJ.ConfFinal(); //Write final configuration

	LJ.rnd.SaveSeed();
	
	
	return 0;
}
