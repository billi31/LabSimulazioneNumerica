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
#include "../include/Genetic_Algorithm.hpp"
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
	("grid", "if present returns a brief analysis on # of generations and # of individuals")
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
	
	
	
	Genetic_Algorithm TSP;
	
	
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
	
	
	
	
	if (vm.count("grid")) {
		std::ofstream out_grid;
		out_grid.open("output/"+area+"/"+TSP.fitness_type+"/grid.dat");
		
		
		for (int i_individuals = 0; i_individuals < 20; i_individuals ++){
			if (area == "circle")
				TSP.n_individuals = i_individuals*10 + 50;
			else if (area == "square")
				TSP.n_individuals = i_individuals*10 + 100;
			
			for (TSP.i_gener = 0; TSP.i_gener < 10; TSP.i_gener ++){
				if (area == "circle")
					TSP.n_generations = TSP.i_gener*50 + 200;
				else if (area == "square")
					TSP.n_generations = TSP.i_gener*50 + 400;
				
				TSP.pop.clear();
				TSP.pop.resize(TSP.n_individuals);
				
				std::cout << "Number of individuals = " << TSP.n_individuals << "\n";
				std::cout << "Number of generations = " << TSP.n_generations << "\n";
				
				TSP.Input(area, "grid"); //Inizialization
				
				
				for (int i=0; i<TSP.n_generations ; i++){
					TSP.NewPop();
					
					for (TSP.i_indiv=0; TSP.i_indiv<TSP.n_individuals; TSP.i_indiv=TSP.i_indiv+2)
						TSP.CrossOver(TSP.i_indiv,TSP.i_indiv+1);
					
					
					for (TSP.i_indiv=0; TSP.i_indiv<TSP.n_individuals; TSP.i_indiv++)
						TSP.Mutations();
					
					
					std::sort(TSP.pop.begin(), TSP.pop.end());
					
					
					TSP.generation++;
				}
				out_grid << TSP.n_individuals << "     " << TSP.n_generations << "     " << TSP.pop[0].fitness << "\n";
			}
		}
		out_grid.close();
	}
	
	
	TSP.Input(area, ""); //Inizialization
	
	
	
	std::ofstream out;
	out.open("output/"+area+"/"+TSP.fitness_type+"/best.fitness.dat");
	
	std::ofstream out_ave;
	out_ave.open("output/"+area+"/"+TSP.fitness_type+"/ave.fitness.dat");
	
	
	for (int i=0; i<TSP.n_generations ; i++){
		std::cout << "Generazione " << TSP.generation << "\n";
		TSP.NewPop();
		//TSP.PrintPopulation();
		
		for (TSP.i_indiv=0; TSP.i_indiv<TSP.n_individuals; TSP.i_indiv=TSP.i_indiv+2){
			TSP.CrossOver(TSP.i_indiv,TSP.i_indiv+1);
			//TSP.PrintPopulation();
			//TSP.CheckIndividuals();
		}
		
		for (TSP.i_indiv=0; TSP.i_indiv<TSP.n_individuals; TSP.i_indiv++){
			TSP.Mutations();
			//TSP.CheckIndividuals();
		}
		
		std::sort(TSP.pop.begin(), TSP.pop.end());
		
		
		TSP.generation++;
		
		out << TSP.pop[0].fitness << "\n";
		out_ave << TSP.FitnessAve() << "\n";
		
		std::cout << "\n---------------------------------\n";
		
	}
	TSP.PrintSolution();
	out.close();
	
	
	TSP.rnd.SaveSeed();
	
	return 0;
}
