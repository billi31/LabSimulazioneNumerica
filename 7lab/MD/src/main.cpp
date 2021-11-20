//#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include "../include/MolDyn_NVE.hpp"
#include "boost/program_options.hpp"


namespace po = boost::program_options;


int main(int argc, const char** argv){
	
	MolDynamics MD;
	
	
	// Declare the supported options.
	po::options_description desc("Allowed options");
	desc.add_options()
	("help", "produce help message")
	("state", po::value<std::string>(), "state of fluid (MUST be gas, liquid or solid)")
	("mode", po::value<std::string>(), "start from fcc or old configuration (MUST be start or repeat)")
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
		MD.state = vm["state"].as<std::string>();
	}
	else {
		std::cout << desc << "\n";
		return 1;
	}
	
	if (vm.count("mode")) {
		if (vm["mode"].as<std::string>() == "start" || vm["mode"].as<std::string>() == "repeat")
		std::cout << "Mode is " << vm["mode"].as<std::string>() << ".\n\n";
		else{
			std::cerr << "Mode not valid\n\n";
			exit(1);
		}
		MD.mode = vm["mode"].as<std::string>();
	}
	else {
		std::cout << desc << "\n";
		return 1;
	}

	std::cout << "\nSimulation Molecular Dynamics \n\n";
	
	MD.Input(); //Inizialization
	
	
	//int nconf = 1;
	int totstep = 1;

	for(MD.iblk=1; MD.iblk <= MD.nblk; ++MD.iblk){ //Simulation
		MD.Reset();   //Reset block averages
		
		for(int istep=1; istep <= MD.nstep; ++istep){
			MD.Move();
			
			MD.Measure();
			MD.Accumulate(); //Update block averages
			
			//if(totstep%MD.iprint == 0) 
				//std::cout << "Number of time-steps: " << totstep << std::endl;
			
			if(totstep%10 == 0){
				MD.PrintMeasure(0);
				//MD.ConfXYZ(nconf);					//Write actual configuration in XYZ format
				//nconf += 1;
			}
			totstep++;
		}
		MD.Averages();   //Print results for current block
	}
	
	MD.ConfOld();
	MD.ConfFinal(); //Write final configuration
	
	MD.rnd.SaveSeed();
	
	return 0;
}
