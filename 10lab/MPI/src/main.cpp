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
#include <algorithm>
#include "../include/Genetic_Algorithm.hpp"
#include <mpi.h>
#include <string>     // std::to_string


// mpirun -np 4 ./Esercizio10_2.exe    comando da eseguire


int main(int argc, char* argv[]){
	int size;
	int rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Status stat[4];
	MPI_Request req, req2;
	
	Genetic_Algorithm TSP;
	
	TSP.ReadInput(rank); //Inizialization

	
	double *x_cities = new double[TSP.ncity];
	double *y_cities = new double[TSP.ncity];
	
	std::ofstream out;
	std::ofstream out_best;

	if (rank == 0){
		out_best.open("output/square/"+TSP.fitness_type+"/best.fitness.rank.dat");

		out.open("output/square/"+TSP.fitness_type+"/cities.dat");

		// Cities on the map
		for (int i=1; i<=TSP.ncity; i++){
			TSP.random_position();
			
			out << i << "    " << TSP.city.x << "    " << TSP.city.y << std::endl;
			x_cities[i-1] = TSP.city.x;
			y_cities[i-1] = TSP.city.y;
		}
		out.close();
	}
	
	//std::stringstream ss;
	//ss << rank;
	std::string srank = std::to_string(rank);	
	
	//send to everyone x_cities and y_cities of rank 0
	MPI_Bcast(x_cities,TSP.ncity,MPI_DOUBLE,0, MPI_COMM_WORLD);		
	MPI_Bcast(y_cities,TSP.ncity,MPI_DOUBLE,0, MPI_COMM_WORLD);
	
	
	for (int i=1; i<=TSP.ncity; i++){
		TSP.city.x = x_cities[i-1];
		TSP.city.y = y_cities[i-1];
		TSP.cities.insert(std::pair<int,position>(i,TSP.city));
	}
	//every process now have the same list of cities
	
	TSP.PopulationInput();
	
	out.open("output/square/"+TSP.fitness_type+"/"+srank+"/best.fitness.dat");
	
	std::ofstream out_ave;
	out_ave.open("output/square/"+TSP.fitness_type+"/"+srank+"/ave.fitness.dat");
	
	double costs[4];
	double best_fitness;

	for (int i=0; i<TSP.n_generations ; i++){
		TSP.NewPop();
		
		for (TSP.i_indiv=0; TSP.i_indiv<TSP.n_individuals; TSP.i_indiv=TSP.i_indiv+2)
			TSP.CrossOver(TSP.i_indiv,TSP.i_indiv+1);
		
		for (TSP.i_indiv=0; TSP.i_indiv<TSP.n_individuals; TSP.i_indiv++)
			TSP.Mutations();
		
		std::sort(TSP.pop.begin(), TSP.pop.end());
		

		//scambio informazioni ogni nmigr generazioni casualmente tra i 4 rank
		if (TSP.generation%TSP.nmigr == 0){
			MPI_Status stat1, stat2, stat3, stat4,stat5, stat6, stat7, stat8;
			MPI_Request req1,req2,req3,req4;
			int itag1=1; int itag2=2; int itag3=3; int itag4=4; int itag5=5; int itag6=6; int itag7=7; int itag8=8;
			int scambio[3];
			for (int i=0; i<3; i++)
				scambio[i] = 0;
			
			if (rank == 0){
				scambio[0] = int(TSP.rnd.Rannyu(1,4));

				if (scambio[0] == 1){
					scambio[1] = 2;
					scambio[2] = 3;
				}
				else if (scambio[0] == 2){
					scambio[1] = 1;
					scambio[2] = 3;
				}
				else {
					scambio[1] = 1;
					scambio[2] = 2;
				}
				
			}
			//faccio sì che abbiano gli stessi valori in scambio
			MPI_Bcast(scambio,4,MPI_INTEGER,0, MPI_COMM_WORLD);
			
			//copio il miglior percorso
			int *best_tour = new int[TSP.ncity];
			int *best_tour2 = new int[TSP.ncity];
			
			for (int i=0; i<TSP.ncity; i++){
				best_tour[i] = TSP.pop[0].tour[i];
				best_tour2[i] = TSP.pop[1].tour[i];
			}
			
			if (rank == 0){		//invia a scambio[0] e riceve da scambio[0]
				//manda a scambio[0] i due percorsi migliori (vettore di 32 int)
				MPI_Isend(best_tour,TSP.ncity,MPI_INTEGER,scambio[0],itag1,MPI_COMM_WORLD,&req1);
				MPI_Isend(best_tour2,TSP.ncity,MPI_INTEGER,scambio[0],itag5,MPI_COMM_WORLD,&req3);
				MPI_Recv(best_tour,TSP.ncity,MPI_INTEGER,scambio[0],itag2, MPI_COMM_WORLD,&stat2);
				MPI_Recv(best_tour2,TSP.ncity,MPI_INTEGER,scambio[0],itag6, MPI_COMM_WORLD,&stat6);
			}
			if (rank == scambio[0]){		//invia a rank 0 e riceve da rank 0
				//manda a rank 0 i due percorsi migliori (vettore di 32 int)
				MPI_Send(best_tour,TSP.ncity,MPI_INTEGER,0,itag2,MPI_COMM_WORLD);
				MPI_Send(best_tour2,TSP.ncity,MPI_INTEGER,0,itag6,MPI_COMM_WORLD);
				MPI_Recv(best_tour,TSP.ncity,MPI_INTEGER,0,itag1, MPI_COMM_WORLD,&stat1);
				MPI_Recv(best_tour2,TSP.ncity,MPI_INTEGER,0,itag5, MPI_COMM_WORLD,&stat5);
			}
			
			if (rank == scambio[1]){		//invia a scambio[2] e riceve da scambio[2]
				//manda a scambio[2] i due percorsi migliori (vettore di 32 int)
				MPI_Isend(best_tour,TSP.ncity,MPI_INTEGER,scambio[2],itag3,MPI_COMM_WORLD,&req2);
				MPI_Isend(best_tour2,TSP.ncity,MPI_INTEGER,scambio[2],itag7,MPI_COMM_WORLD,&req4);
				MPI_Recv(best_tour,TSP.ncity,MPI_INTEGER,scambio[2],itag4, MPI_COMM_WORLD,&stat4);
				MPI_Recv(best_tour2,TSP.ncity,MPI_INTEGER,scambio[2],itag8, MPI_COMM_WORLD,&stat8);
			}
			if (rank == scambio[2]){		//invia a scambio[1] e riceve da scambio[1]
				//manda a scambio[1] i due percorsi migliori (vettore di 32 int)
				MPI_Send(best_tour,TSP.ncity,MPI_INTEGER,scambio[1],itag4,MPI_COMM_WORLD);
				MPI_Send(best_tour2,TSP.ncity,MPI_INTEGER,scambio[1],itag8,MPI_COMM_WORLD);
				MPI_Recv(best_tour,TSP.ncity,MPI_INTEGER,scambio[1],itag3, MPI_COMM_WORLD,&stat3);
				MPI_Recv(best_tour2,TSP.ncity,MPI_INTEGER,scambio[1],itag7, MPI_COMM_WORLD,&stat7);
			}
			
			for (int i=0; i<TSP.ncity; i++){
				TSP.pop[0].tour[i] = best_tour[i];
				TSP.pop[1].tour[i] = best_tour2[i];
			}
		}
		if (rank == 0){
			std::cout << "Generazione " << TSP.generation << "\n";
			std::cout << "\n---------------------------------\n";
		}
		
		TSP.pop[0].fitness = TSP.CostFunction(TSP.pop[0].tour);
		TSP.pop[1].fitness = TSP.CostFunction(TSP.pop[1].tour);
		std::sort(TSP.pop.begin(), TSP.pop.end());
		
		
		double fitness;
		fitness = TSP.pop[0].fitness;
		
		for (int i=0;i<4; i++)
			costs[i] = 0;
		
		MPI_Gather(&fitness,1,MPI_DOUBLE,costs,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		if(rank==0){			
			best_fitness = *std::min_element(std::begin(costs), std::end(costs));
			out_best << best_fitness << "\n";
		}
		
		out << TSP.pop[0].fitness << "\n";
		out_ave << TSP.FitnessAve() << "\n";
		
		TSP.generation++;
	}
	
	//stampo le soluzioni di ogni rank
	TSP.PrintSolution(rank);
	out.close();
	
	if(rank==0){
		int best_rank;
		best_rank = std::distance(std::begin(costs), std::min_element(std::begin(costs), std::end(costs)));
		out_best.close();		
		
		std::string sbest_rank = std::to_string(best_rank);
		std::cout << "\n\n\nBEST PATH FOUND IS IN FILE " << "output/square/"+TSP.fitness_type+"/"+sbest_rank+"/solution.dat" << std::endl;
	}
	
	
	TSP.rnd.SaveSeed();
	
	MPI_Finalize();
	
	return 0;
}
