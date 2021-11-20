//
//  1.cpp
//  lab1
//
//  Created by Administrator on 18/03/21.
//
//

#include "../include/lib.hpp"
#include "../include/random.hpp"


using namespace std;

int main (int argc, char *argv[]){
	
	int M = 1000000;					//numero di throw
	int N = 100;						//numero di blocchi
	int L = int(M/N);					//numero di throw in ogni blocco
	
	
	//definisco una variabile Random che utilizzerò per generare i numeri casuali
	Random rnd;
	rnd.Initialize(rnd);
	
	double *ave = new double[N]();
	double *err_prog = new double[N]();
	
	
	//inizio un ciclo sul numero di blocchi
	for (int i=0; i<N; i++){
		double sum = 0;
		
		//genero L numeri random per ogni blocco
		for (int j=0; j<L; j++){
			double appo = rnd.Rannyu();
			sum += appo;
		}
		ave[i] = sum/L;
	}
	
	block_ave("mean.txt", N, ave, err_prog);

	
	delete[] ave;
	delete[] err_prog;	
	
	
	
	
	//rieseguo gli stessi passaggi, ma valutando la varianza dei dati
	double *var = new double[N]();
	double *varerr_prog = new double[N]();
	
	//inizio un ciclo sul numero di blocchi
	//i valori in var saranno le medie degli scarti quadratici per gli L numeri generati
	
	for (int i=0; i<N; i++){
		double sum = 0;
		
		//genero L numeri random per ogni blocco
		for (int j=0; j<L; j++){
			double appo = rnd.Rannyu();
			sum += pow((appo-0.5),2);
		}
		var[i] = sum/L;
	}	
	
	block_ave("var.txt", N, var, varerr_prog);

	delete[] var;
	delete[] varerr_prog;
	
	
	
	
	//a questo punto divido l'intervallo [0,1) in M=100 parti e considero n=10000 punti casuali.
	//per le n/M=100 iterazioni calcolo il chi quadro e ne mostro l'andamento.
	
	M = 100;				//numero di intervalli
	int n_chi = 100;			//numero di ripetizioni di calcolo di chi quadro
	
	int n = 10000;			//numero di punti che considero di volta in volta
	L = int(n/M);			//numero di punti che mi aspetto di ottenere in ogni intervallo
	
	
	
	//inizializzo il vettore dei chi quadri a 0
	double *chi_sq = new double[n_chi] ();
	
	//inizio un ciclo per il chi quadro da calcolare
	for (int k=0; k<n_chi; k++){
		
		//creo un array che conterrà il numero di punti nell'intervallo i-esimo per ogni ripetizione k-esima
		int *n_conta = new int[M]();
		
		//ciclo sugli n numeri casuali del file
		for (int j=0; j<n; j++){
			double appo = rnd.Rannyu();
			
			//associo ad ogni numero random in [0,1) un int tra 0 e M, che indicherà l'intervallo di appartenenza
			int i = int(appo*M);
			
			n_conta[i]++;
		}
		for (int i=0; i<M; i++)
			chi_sq[k] += pow((n_conta[i]-L),2);
		
		delete[] n_conta;
		
		chi_sq[k] /= L;
		
	}
	
	//scrivo i risultati nel file "output/chi_sq.txt"
	outfile_vett("chi_sq.txt",n_chi,chi_sq);
	
	
	delete[] chi_sq;
	
	rnd.SaveSeed();
	
	return 0;
}
