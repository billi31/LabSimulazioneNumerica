//
//  main.cpp
//  5lab
//
//  Created by Administrator on 26/07/21.
//
//

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cmath>        // rint, pow
#include "../include/Hydrogen.hpp"
#include "boost/program_options.hpp"

using namespace std;

namespace po = boost::program_options;


inline bool file_exists(const char* name);



int main(int argc, const char** argv){
	
	
	Hydrogen hyd;
	
	std::string type = "null";
	
	
	// Declare the supported options.
	po::options_description desc("Allowed options");
	desc.add_options()
	("help", "produce help message")
	("type", po::value<string>(), "evaluate radius or visualise wave function (MUST be radius or graph)")
	("state", po::value<string>(), "set Hydrogen state (can be 100 or 210)")
	("mode", po::value<string>(), "set probability distribution of Metropolis (can be unif or gauss)")
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
		if (vm["state"].as<string>() == "100" || vm["state"].as<string>() == "210")
	 		std::cout << "Hydrogen' state to be analysed is " << vm["state"].as<string>() << ".\n\n";
		else{
	 		std::cerr << "Hydrogen's state not valid\n\n";
			exit(1);
		}
		hyd.state = vm["state"].as<string>();
	}
	
	if (vm.count("mode")) {
		if (vm["mode"].as<string>() == "unif" || vm["mode"].as<string>() == "gauss")
	 		std::cout << "Metropolis' probability distribution function to be used is " << vm["mode"].as<string>() << ".\n\n";
		else{
	 		std::cerr << "Metropolis' probability distribution function not valid\n\n";
			exit(1);
		}
		hyd.mode = vm["mode"].as<string>();
	}
	
	if (vm.count("type")) {
		if (vm["type"].as<string>() == "radius" || vm["type"].as<string>() == "graph")
			std::cout << "Program option is " << vm["type"].as<string>() << ".\n\n";
		else{
			std::cerr << "Program option not valid\n\n";
			exit(1);
		}
		type = vm["type"].as<string>();
	}
	else {
		std::cout << desc << "\n";
		return 1;
	}
	
	
	ofstream clearfile;
	
	if (type == "graph"){
		int nstep_graph = nstep/10;
		if (hyd.state == "null"){
			if (hyd.mode == "null"){
				hyd.Initialize(0.,0.,1.,"unif","100");
				
				clearfile.open("output/graph.100unif.dat");
				clearfile.close();
				
				//generates nstep moves of the algorithm
				for (int i=0; i<nstep_graph; i++){
					hyd.Move();
					hyd.PrintPos();
				}
				
				std::cout << "Acceptance of " << hyd.Acceptance() << endl << endl;
				
				
				//gauss 100
				hyd.Initialize(0.,0.,1.,"gauss","100");
				
				clearfile.open("output/graph.100gauss.dat");
				clearfile.close();
				
				//generates nstep moves of the algorithm
				for (int i=0; i<nstep_graph; i++){
					hyd.Move();
					hyd.PrintPos();
				}
				
				std::cout << "Acceptance of " << hyd.Acceptance() << endl << endl;
				
				
				//unif 210
				hyd.Initialize(0.,0.,1.,"unif","210");
				
				clearfile.open("output/graph.210unif.dat");
				clearfile.close();
				
				//generates nstep moves of the algorithm
				for (int i=0; i<nstep_graph; i++){
					hyd.Move();
					hyd.PrintPos();
				}
				
				std::cout << "Acceptance of " << hyd.Acceptance() << endl << endl;
				
				
				//gauss 210
				hyd.Initialize(0.,0.,1.,"gauss","210");
				
				clearfile.open("output/graph.210gauss.dat");
				clearfile.close();
				
				//generates nstep moves of the algorithm
				for (int i=0; i<nstep_graph; i++){
					hyd.Move();
					hyd.PrintPos();
				}
				
				std::cout << "Acceptance of " << hyd.Acceptance() << endl << endl;
				
			}
			
			else {
				
				//mode 100
				hyd.Initialize(0.,0.,1.,hyd.mode,"100");
				
				clearfile.open("output/graph.100" + hyd.mode + ".dat");
				clearfile.close();
				
				//generates nstep moves of the algorithm
				for (int i=0; i<nstep_graph; i++){
					hyd.Move();
					hyd.PrintPos();
				}
				
				std::cout << "Acceptance of " << hyd.Acceptance() << endl << endl;
				
				
				//mode 210
				hyd.Initialize(0.,0.,1.,hyd.mode,"210");
				
				clearfile.open("output/graph.210" + hyd.mode + ".dat");
				clearfile.close();
				
				//generates nstep moves of the algorithm
				for (int i=0; i<nstep_graph; i++){
					hyd.Move();
					hyd.PrintPos();
				}
				
				std::cout << "Acceptance of " << hyd.Acceptance() << endl << endl;
				
			}
			
		}
		else {
			if (hyd.mode == "null"){
				//unif state
				hyd.Initialize(0.,0.,1.,"unif",hyd.state);
				
				clearfile.open("output/graph." + hyd.state + "unif.dat");
				clearfile.close();
				
				//generates nstep moves of the algorithm
				for (int i=0; i<nstep_graph; i++){
					hyd.Move();
					hyd.PrintPos();
				}
				
				std::cout << "Acceptance of " << hyd.Acceptance() << endl << endl;
				
				
				//gauss state
				hyd.Initialize(0.,0.,1.,"gauss",hyd.state);
				
				clearfile.open("output/graph." + hyd.state + "gauss.dat");
				clearfile.close();
				
				//generates nstep moves of the algorithm
				for (int i=0; i<nstep_graph; i++){
					hyd.Move();
					hyd.PrintPos();
				}
				
				std::cout << "Acceptance of " << hyd.Acceptance() << endl << endl;
				
			}
			
			else {
				
				hyd.Initialize(0.,0.,1.,hyd.mode,hyd.state);
				
				clearfile.open("output/graph." + hyd.state + hyd.mode + ".dat");
				clearfile.close();
				
				//generates nstep moves of the algorithm
				for (int i=0; i<nstep_graph; i++){
					hyd.Move();
					hyd.PrintPos();
				}
				
				std::cout << "Acceptance of " << hyd.Acceptance() << endl << endl;
				
				
			}
		}
	}
	
	
	if (type == "radius"){
		
		int L = int(nstep/nblock);					//number of evaluations in each block
		
		if (hyd.state == "null"){
			if (hyd.mode == "null"){
				int istep = 1;

				hyd.Initialize(0.,0.,1.,"unif","100");
				
				//inizio un ciclo sul numero di blocchi
				for (hyd.iblk=1; hyd.iblk<=nblock; hyd.iblk++){		
					hyd.Reset();

					//ciclo sugli L passi in ogni blocco
					for (int j=0; j<L; j++){
						hyd.Move();
						hyd.Accumulate();
						
						if(istep%100 == 0)
							hyd.PrintInstantaneusRadius();
						istep++;
						
					}
					
					hyd.Averages();

				}
				
				istep = 1;
				//gauss 100
				hyd.Initialize(0.,0.,1.,"gauss","100");
				
				//inizio un ciclo sul numero di blocchi
				for (hyd.iblk=1; hyd.iblk<=nblock; hyd.iblk++){		
					hyd.Reset();

					//ciclo sugli L passi in ogni blocco
					for (int j=0; j<L; j++){
						hyd.Move();
						hyd.Accumulate();
						
						if(istep%100 == 0)
							hyd.PrintInstantaneusRadius();
						istep++;
						
					}
					
					hyd.Averages();
				}
				
				istep = 1;
				
				//unif 210
				hyd.Initialize(0.,0.,1.,"unif","210");
				
				//inizio un ciclo sul numero di blocchi
				for (hyd.iblk=1; hyd.iblk<=nblock; hyd.iblk++){		
					hyd.Reset();

					//ciclo sugli L passi in ogni blocco
					for (int j=0; j<L; j++){
						hyd.Move();
						hyd.Accumulate();
						
						if(istep%100 == 0)
							hyd.PrintInstantaneusRadius();

						istep++;
						
					}
					
					hyd.Averages();

				}
				
				istep = 1;
				
				//gauss 210
				hyd.Initialize(0.,0.,1.,"gauss","210");
				
				//inizio un ciclo sul numero di blocchi
				for (hyd.iblk=1; hyd.iblk<=nblock; hyd.iblk++){		
					hyd.Reset();

					//ciclo sugli L passi in ogni blocco
					for (int j=0; j<L; j++){
						hyd.Move();
						hyd.Accumulate();
						
						if(istep%100 == 0)
							hyd.PrintInstantaneusRadius();

						istep++;
						
					}
					
					hyd.Averages();
					
				}
				
			}
			else {
				int istep = 1;
				
				//mode 100
				hyd.Initialize(0.,0.,1.,hyd.mode,"100");
				
				//inizio un ciclo sul numero di blocchi
				for (hyd.iblk=1; hyd.iblk<=nblock; hyd.iblk++){		
					hyd.Reset();
					
					//ciclo sugli L passi in ogni blocco
					for (int j=0; j<L; j++){
						hyd.Move();
						hyd.Accumulate();
						
						if(istep%100 == 0)
							hyd.PrintInstantaneusRadius();
						
						istep++;
						
					}
					
					hyd.Averages();
					
					
				}
				
				istep = 1;
				
				//mode 210
				hyd.Initialize(0.,0.,1.,hyd.mode,"210");
				
				
				//inizio un ciclo sul numero di blocchi
				for (hyd.iblk=1; hyd.iblk<=nblock; hyd.iblk++){		
					hyd.Reset();

					//ciclo sugli L passi in ogni blocco
					for (int j=0; j<L; j++){
						hyd.Move();
						hyd.Accumulate();
						
						if(istep%100 == 0)
							hyd.PrintInstantaneusRadius();
						
						istep++;
						
					}
					
					hyd.Averages();
					
				}
				
			}
		}
		else {
			if (hyd.mode == "null"){
				int istep = 1;

				//unif state
				hyd.Initialize(0.,0.,1.,"unif",hyd.state);
				
				//inizio un ciclo sul numero di blocchi
				for (hyd.iblk=1; hyd.iblk<=nblock; hyd.iblk++){		
					hyd.Reset();

					//ciclo sugli L passi in ogni blocco
					for (int j=0; j<L; j++){
						hyd.Move();
						hyd.Accumulate();
						
						if(istep%100 == 0)
							hyd.PrintInstantaneusRadius();

						istep++;
						
					}
					
					hyd.Averages();
				
				}
				
				istep = 1;
				
				//gauss state
				hyd.Initialize(0.,0.,1.,"gauss",hyd.state);
				
				//inizio un ciclo sul numero di blocchi
				for (hyd.iblk=1; hyd.iblk<=nblock; hyd.iblk++){		
					hyd.Reset();

					//ciclo sugli L passi in ogni blocco
					for (int j=0; j<L; j++){
						hyd.Move();
						hyd.Accumulate();
						
						if(istep%100 == 0)
							hyd.PrintInstantaneusRadius();
						
						istep++;
						
					}
					
					hyd.Averages();
					
				}
			}
			
			else {
				int istep = 1;

				//gauss state
				hyd.Initialize(0.,0.,1.,hyd.mode,hyd.state);
				
				//inizio un ciclo sul numero di blocchi
				for (hyd.iblk=1; hyd.iblk<=nblock; hyd.iblk++){		
					hyd.Reset();

					//ciclo sugli L passi in ogni blocco
					for (int j=0; j<L; j++){
						hyd.Move();
						hyd.Accumulate();
						
						if(istep%100 == 0)
							hyd.PrintInstantaneusRadius();

						istep++;
						
					}
					
					hyd.Averages();
					
				}
			}
		}
	}
	
	
	
	return 0;
}



//check if the file exists
inline bool file_exists(const char* name){
	ifstream f(name);
	return f.good();
}
