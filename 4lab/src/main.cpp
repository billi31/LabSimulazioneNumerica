//#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include "../include/MolDyn_NVE.hpp"
#include "../include/random.hpp"

using namespace std;

inline bool file_exists(const char* name);

int main(int argc, const char** argv){
	
	if (argc!=2){
		cout << "\nSimulation Molecular Dynamics, specify the desired mode:\n\n";
		cout << "	1. \"start\" if you want to start from the spatial configuration in file \"config.0\";\n\n";
		cout << "	2. \"repeat\" if you want to start from the two spatial configuration in file \"old.0\" and \"config.0\".\n\n";
		exit(1);
	}
	
	cout << "\nSimulation Molecular Dynamics \n\n";
	
	
	MolDynamics MD;
	
	
	if (string(argv[1])=="start"){
		
		if (!file_exists("input.dat")){
			cerr << "Unable to open \"input.dat\"\n"; 
			exit(1);
		}
		if (!file_exists("config.0")){
			cerr << "Unable to open \"config.0\"\n"; 
			exit(1);
		}
		cout << "Simulation of Molecular Dynamics: the setup is in file \"input.dat\" and the spatial configuration is in file \"config.0\"\n\n";
		
		MD.Input("input.dat","config.0");				//Inizialization
		
		
	}
	
	
	else if (string(argv[1])=="repeat"){
		
		if (!file_exists("input.dat")){
			cerr << "Unable to open \"input.dat\"\n"; 
			exit(1);
		}
		if (!file_exists("old.0")){
			cerr << "Unable to open \"old.0\", it is necessary to run the code in \"restart\" mode\n"; 
			exit(1);
		}
		if (!file_exists("config.0")){
			cerr << "Unable to open \"config.0\", it is necessary to run the code in \"restart\" mode\n"; 
			exit(1);
		}
		cout << "Simulation of Molecular Dynamics: the setup is in file \"input.dat\" and the spatial configurations are in file \"old.0\" and in file \"config.0\"\n\n";
		
		MD.Input("input.dat","config.0","old.0");		//Inizialization
		
	}
	
	else {
		cout << "\nSimulation Molecular Dynamics, specify the desired mode:\n\n";
		cout << "	1. \"start\" if you want to start from the spatial configuration in file \"config.0\";\n\n";
		cout <<	"	2. \"repeat\" if you want to start from the two spatial configuration in file \"old.0\" and \"config.0\".\n\n";
		exit(1);
	}
	
	int L = int(MD.nstep/nblock);					//number of evaluations in each block
	
	//int nconf = 1;
	//int istep = 1;
	
	//inizio un ciclo sul numero di blocchi
	for (MD.iblk=1; MD.iblk<=nblock; MD.iblk++){		
		MD.Reset();
		
		//ciclo sugli L passi in ogni blocco
		for (int j=0; j<L; j++){
			MD.Move();
			MD.Measure();
			MD.Accumulate();
			
			/*
			if(istep%MD.iprint == 0) 
				cout << "Number of time-steps: " << istep << endl;
			*/
			
			//if(istep%10 == 0){
				//MD.PrintMeasure();
				//MD.ConfXYZ(nconf);					//Write actual configuration in XYZ format
				//nconf += 1;
			//}
			
			//istep++;
			
		}
		
		MD.Averages();
		
	}
	MD.ConfOld();
	MD.ConfFinal();								//Write final configuration to restart
	
	return 0;
}


//check if the file exists
inline bool file_exists(const char* name){
	ifstream f(name);
	return f.good();
}
