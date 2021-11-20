//
//  main.cpp
//  3lab
//
//  Created by Administrator on 12/04/21.
//
//



#include "../include/lib.hpp"
#include "../include/random.hpp"


using namespace std;

int main (int argc, char *argv[]){
	
	int M = 100000;						//numero di throw (numero effettivo di price calcolati)
	int N = 100;						//numero di blocchi
	int L = int(M/N);					//numero di throw in ogni blocco
	
	
	//definisco una variabile Random che utilizzerò per generare i numeri casuali
	Random rnd;
	rnd.Initialize(rnd);
	
	
	//definisco i parametri da utilizzare
	double S0 = 100.;
	double K = 100.;
	double T = 1.;
	double r = 0.1;
	double sigma = 0.25;
	double t = 0;
	
	
	
	
	double *ave_C = new double[N]();
	double *err_prog_C = new double[N]();
	
	double *ave_P = new double[N]();
	double *err_prog_P = new double[N]();
		
	
	//riempio ave e av2 con le medie degli L valori di CALL per il caso diretto
	//riempio ave_P e av2_P con le medie degli L valori di PUT per il caso diretto
	
	//apro un ciclo sul numero di blocchi
	for (int i=0; i<N; i++){
		
		double sum_C = 0;
		double sum_P = 0;
		
		//apro un ciclo sul numero di valutazioni del prezzo di CALL
		for (int j=0; j<L; j++){
			
			double ST;
			double CT;
			double PT;
			
			//valuto il prezzo al tempo T
			ST = price(T, S0, r, sigma, rnd);
			
			//valuto il prezzo di CALL-option
			CT = exp(-r*T) * max(0.,ST-K);
			
			//valuto il prezzo di PUT-option
			PT = exp(-r*T) * max(0.,K-ST);
			
			//aggiorno sum_C con la somma dei valori di prezzo di CALL-option appena calcolati
			sum_C += CT;
			//aggiorno sum_P con la somma dei valori di prezzo di PUT-option appena calcolati
			sum_P += PT;
			
		}
		//valuto la media degli L valori di prezzo di CALL_option calcolati per l'i-esimo blocco e ne calcolo il quadrato
		ave_C[i] = sum_C/L;	
		
		//valuto la media degli L valori di prezzo di PUT_option calcolati per l'i-esimo blocco e ne calcolo il quadrato
		ave_P[i] = sum_P/L;
		
	}
	
	block_ave("CALL.txt", N, ave_C, err_prog_C);
	block_ave("PUT.txt", N, ave_P, err_prog_P);
	
	
	delete[] ave_C;
	delete[] err_prog_C;	
	delete[] ave_P;
	delete[] err_prog_P;
	
	double *ave_CD = new double[N]();
	double *err_prog_CD = new double[N]();
	
	double *ave_PD = new double[N]();
	double *err_prog_PD = new double[N]();
	
	
	//definisco il numero di passi per il calcolo di S(t_k)
	int passi = 100;
	double delta_t = (double) T/passi;
		
	//riempio ave_CD con le medie degli L valori di CALL per il caso discreto
	//riempio ave_PD con le medie degli L valori di PUT per il caso discreto

	//apro un ciclo sul numero di blocchi
	for (int i=0; i<N; i++){
		
		double sum_CD = 0;
		double sum_PD = 0;
		
		//apro un ciclo sul numero di valutazioni del prezzo di CALL
		for (int j=0; j<L; j++){
			
			double ST = S0;
			double CT;
			double PT;
			
			//apro un ciclo sul valore del prezzo S(t_{k+1})
			for (int k=0; k<passi; k++){
				//aggiorno il prezzo al tempo t_{k+1}
				ST = price_discr(t + k*delta_t, t + (k+1)*delta_t, ST, r, sigma, rnd);
			}
			
			//valuto il prezzo di CALL-option
			CT = exp(-r*T) * max(0.,ST-K);
			//valuto il prezzo di PUT-option
			PT = exp(-r*T) * max(0.,K-ST);
			
			//aggiorno sum con la somma dei valori di prezzo di CALL-option appena calcolati
			sum_CD += CT;
			//aggiorno sum_P con la somma dei valori di prezzo di PUT-option appena calcolati
			sum_PD += PT;
			
		}
		//valuto la media degli L valori di prezzo di CALL_option calcolati per l'i-esimo blocco
		ave_CD[i] = sum_CD/L;		
		//valuto la media degli L valori di prezzo di PUT_option calcolati per l'i-esimo blocco
		ave_PD[i] = sum_PD/L;		
	}
	
	block_ave("CALL_d.txt", N, ave_CD, err_prog_CD);
	block_ave("PUT_d.txt", N, ave_PD, err_prog_PD);
	
	
	delete[] ave_CD;
	delete[] err_prog_CD;
	delete[] ave_PD;
	delete[] err_prog_PD;
	
	rnd.SaveSeed();
	
	return 0;
}

