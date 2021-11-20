//
//  main.cpp
//  10lab
//
//  Created by Administrator on 30/08/21.
//

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cmath>        // rint, pow
#include "../include/Simulated_Annealing.hpp"
#include "boost/program_options.hpp"



namespace po = boost::program_options;


int main(int argc, const char** argv){
	
	
	std::string area;
	// Declare the supported options.
	po::options_description desc("Allowed options");
	desc.add_options()
	("help", "produce help message")
	("area",po::value<std::string>(), "shape of the area in which the cities will be (MUST be circle or square)")
	("fitness",po::value<std::string>(), "cost function to be used (MUST be L1 or L2)")
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
	
	
	if (vm.count("area")) {
		if (vm["area"].as<std::string>() == "circle" || vm["area"].as<std::string>() == "square")
			std::cout << "Area to be used is " << vm["area"].as<std::string>() << ".\n\n";
		else{
			std::cerr << "Area not valid\n\n";
			exit(1);
		}
		area = vm["area"].as<std::string>();
	}
	else {
		std::cout << desc << "\n";
		return 1;
	}
	
	
	
	Simulated_Annealing TSP;
	
	
	if (vm.count("fitness")) {
		if (vm["fitness"].as<std::string>() == "L1" || vm["fitness"].as<std::string>() == "L2")
			std::cout << "Cost function to be used is " << vm["fitness"].as<std::string>() << ".\n\n";
		else{
			std::cerr << "Cost function not valid\n\n";
			exit(1);
		}
		TSP.fitness_type = vm["fitness"].as<std::string>();
	}
	else {
		std::cout << desc << "\n";
		return 1;
	}
	
	TSP.Input(area); //Inizialization
	
	std::ofstream out;
	out.open("output/"+area+"/"+TSP.fitness_type+"/best.fitness.dat");
	
	TSP.t = TSP.initial_t;
	
	for (int i=0; i<TSP.step_temperature; i++){
		TSP.t = TSP.t * (1 - TSP.cooling_rate);
		TSP.beta = 1./TSP.t;
		std::cout << "\nTemperature = " << TSP.t;
		
		TSP.Reset();
		
		for (int istep=0; istep<TSP.MC_step; istep++)
			TSP.Move();
		
		out << TSP.t << "    " << TSP.indiv.fitness << "\n";
		std::cout << "\nAcceptance = " << (double)TSP.accepted/TSP.attempted;		
		std::cout << "\n\n---------------------------------\n";
		
	}
	TSP.PrintSolution();
	out.close();
	
	
	TSP.rnd.SaveSeed();
	
	return 0;
}
