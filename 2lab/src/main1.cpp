//
//  main1.cpp
//  lab2
//
//  Created by Administrator on 26/03/21.
//
//


#include "../include/lib.hpp"
#include "../include/random.hpp"


using namespace std;

int main (int argc, char *argv[]){
	
	int M = 10000;						//numero di throw (numero effettivo di integrali calcolati)
	int N = 100;						//numero di blocchi
	int L = int(M/N);					//numero di throw in ogni blocco
	int R = 50000;						//numero di variabili random da utilizzare nel calcolo dell'integrale
	
	
	//definisco una variabile Random che utilizzerò per generare i numeri casuali
	Random rnd;
	rnd.Initialize(rnd);
	
	
	double *ave = new double[N]();
	double *err_prog = new double[N]();
	
	
	//riempio ave con le medie degli L integrali valutati ognuno con R numeri random, per N volte
	
	//apro un ciclo sul numero di blocchi
	for (int i=0; i<N; i++){
		
		double sum = 0;
		
		//apro un ciclo sul numero di valutazioni dell'integrale
		for (int j=0; j<L; j++){
			
			double sum_int = 0;
			
			//apro un ciclo sul numero di punti da utilizzare per il calcolo dell'integrale
			for (int k=0; k<R; k++){
				double gx;
				
				//valuto l'integranda in un numero distribuito uniformemente tra 0 e 1
				gx = eval_cos(rnd.Rannyu());
				
				//sommo i valori appena trovati
				sum_int += gx;
				
			}
			
			//valuto la media dell'integranda valutata a partire da R numeri distribuiti uniformemente tra 0 e 1
			sum_int /= R;
			
			//aggiorno sum con la somma dei valori degli integrali appena valutati 
			sum += sum_int;
			
		}
		//valuto la media degli L integrali calcolati per l'i-esimo blocco e ne calcolo il quadrato
		ave[i] = sum/L;		
	}
	
	
	//svolgo la media a blocchi e scrivo su file i risultati	
	block_ave("int_unif.txt", N, ave, err_prog);
	
	
	delete[] ave;
	delete[] err_prog;
	
	
	
	//affronto lo stesso problema, ma ora utilizzo l'importance sampling
	double *ave_IS = new double[N]();
	double *err_prog_IS = new double[N]();
	
	//voglio che d(x)dx sia molto simile a cos(x*pi/2)*pi/2, a questo punto, decido di generare dei numeri casuali
	//distribuiti secondo p(x)=CN*(1-x)*pi/2 dove CN è la costante di normalizzazione, pari a 4/pi
	//la retta (1-x)pi/2 è infatti la retta le cui intersezioni con gli assi coincidono con quelle di cos(x*pi/2)*pi/2
	
	
	
	//per fare questo, devo calcolare l'inverso della cumulata di p(x), che porta a scrivere x=1-sqrt(1-y),
	//dove y è un numero distribuito uniformemente tra 0 e 1
	
	//in questo modo, la valutazione dell'integrale sarà data dalla media della funzione pi*cos(x*pi/2)/(4(1-x))
	//dove x è distribuita come sopra
	
	//riempio ave_IS con le medie degli L integrali valutati ognuno con R numeri random, per N volte
	
	//apro un ciclo sul numero di blocchi
	for (int i=0; i<N; i++){
		
		double sum = 0;
		
		//apro un ciclo sul numero di valutazioni dell'integrale
		for (int j=0; j<L; j++){
			
			double sum_int = 0;
						
			//apro un ciclo sul numero di punti da utilizzare per il calcolo dell'integrale (R)
			for (int k=0; k<R; k++){
				double gx;
				double x = rnd.cos_IS();
				
				//valuto l'integranda in un numero distribuito secondo la retta 2(1-x)
				gx = eval_cos_IS(x);
				
				//sommo i valori appena trovati
				sum_int += gx;
				
			}
			
			//valuto la media dell'integranda valutata a partire da R numeri distribuiti secondo la retta 2(1-x)
			sum_int /= R;
			
			//aggiorno sum con la somma dei valori degli integrali appena valutati 
			sum += sum_int;
			
		}
		
		//valuto la media degli L integrali calcolati per l'i-esimo blocco e ne calcolo il quadrato
		ave_IS[i] = sum/L;
		
	}
	
	
	//svolgo la media a blocchi e scrivo su file i risultati	
	block_ave("int_cos_IS.txt", N, ave_IS, err_prog_IS);
	
	
	delete[] ave_IS;
	delete[] err_prog_IS;
	
	rnd.SaveSeed();
	
	
	return 0;
}

