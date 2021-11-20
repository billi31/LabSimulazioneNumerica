//
//  3.cpp
//  lab1
//
//  Created by Administrator on 23/03/21.
//
//

#include "../include/lib.hpp"
#include "../include/random.hpp"


using namespace std;

int main (int argc, char *argv[]){
	
	//definisco una variabile Random che utilizzerò per generare i numeri casuali
	Random rnd;
	rnd.Initialize(rnd);
	
	
	//definisco il numero di aghi che lancerò, il numero di blocchi ed il numero di punti in ogni blocco
	int M = 10000;					//numero totale di valutazioni di pi greco
	int N = 100;					//numero di blocchi
	int L = int(M/N);				//numero di valutazioni in ogni blocco
	
	int N_tot = 50000;				//numero di lanci per ogni valutazione
	
	
	//definisco i valori di d ed l (ovvero la separazione tra le righe e la lunghezza dell'ago)
	double d = 3.;
	double l = 2.5;	 
	
	double *ave = new double[N]();
	double *err_prog = new double[N]();
	
	
	//apro un ciclo sul numero di blocchi
	for (int i=0; i<N; i++){
		
		double sum = 0;
		
		//Apro un ciclo sul numero di valutazioni di pi greco
		for (int j=0; j<L; j++){
			
			int N_theta = 0;
			int N_tocco = 0;
			
			//faccio girare il programma fino a che gli angoli theta generati siano N_tot
			while(N_theta < N_tot){
				double x_j;
				double y_j;
				double sum_2;
				
				//genero due punti casuali x_j e y_j tra 0 ed 1: se la loro somma in quadratura è inferiore a 1, ricavo il valore del seno di theta
				//in caso contrario, rigenero i due punti casuali
				x_j = rnd.Rannyu();
				y_j = rnd.Rannyu();
				
				//calcolo la somma in quadratura
				sum_2 = pow(x_j,2) + pow(y_j,2);
				if(sum_2 <= 1){
					//genero l'ascissa del centro dell'ago lanciato, che avrà un valore compreso tra 0 e d/2
					double y;
					y = rnd.Rannyu(0,d*0.5);
					
					double sintheta;
					sintheta = y_j/sqrt(sum_2);
					N_theta++;
					
					//se il seno dell'angolo sarà maggiore di d/2-y, allora l'ago tocca la riga
					if (sintheta*l*0.5 > d*0.5-y)
						N_tocco++;
					
				}
			}
			
			sum += 1. / N_tocco;
			
		}
		
		ave[i] = (2.*l*sum*N_tot) / (L*d);	
	}
	
	block_ave("pi.txt", N, ave, err_prog);
	
	
	delete[] ave;
	delete[] err_prog;	
	
	rnd.SaveSeed();
	
	return 0;
}

	
	
	
