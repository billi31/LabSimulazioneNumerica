//
//  main2.cpp
//  lab2
//
//  Created by Administrator on 26/03/21.
//
//


#include "../include/lib.hpp"
#include "../include/random.hpp"


using namespace std;

int main (int argc, char *argv[]){
	
	
	//definisco una variabile Random che utilizzerò per generare i numeri casuali
	Random rnd;
	rnd.Initialize(rnd);
	
	
	
	
	int n = 100;				//numero massimo di passi
	int M = 100000;				//numero di random walk
	int N = 100;				//numero di blocchi
	int L = int(M/N);			//numero di random walk per ogni blocco
	
	
	
	
	
	
	
	
	double *passi_d = new double[n+1];				//numero del passo (da 0 a 100 compresi)
	double *sum_prog_d = new double[n+1]();			//contiene la somma progressiva delle medie di r_i^2 valutata su tutti gli N blocchi
	double *su2_prog_d = new double[n+1]();
	double *err_prog_d = new double[n+1]();			//contiene l'errore per il passo i-esimo valutato considerando tutti gli N blocchi
	
	
	
	//apro un ciclo sul numero di blocchi
	for (int j=0; j<N; j++){
		
		double *ave_d = new double[n+1]();			//contiene le medie per il passo i-esimo relativo al blocco j-esimo
		double *av2_d = new double[n+1]();
		
		//apro un ciclo sul numero di RW all'interno del blocco j-esimo
		for (int k=0; k<L; k++){
			double rx = 0;
			double ry = 0;
			double rz = 0;
			
			//apro un ciclo sul numero di passi da effettuare per il RW k-esimo
			for (int i=0; i<=n; i++){
				ave_d[i] += pow(rx,2) + pow(ry,2) + pow(rz,2);
				//effettuo a posteriori il passo, in modo tale da avere per ave[0] una distanza quadratica media nulla
				//in ave[n] avremo la distanza quadratica media al 100-esimo passo
				rnd.Walk_d(rx,ry,rz);
			}
			
		}
		
		//ave_d[i] contiene la media della distanza quadratica all'i-esimo passo per gli L RW del j-esimo blocco
		for (int i=0; i<=n; i++){
			ave_d[i] /= L;
			av2_d[i] = pow(ave_d[i],2);
			sum_prog_d[i] += ave_d[i];
			su2_prog_d[i] += av2_d[i];
		}
		
		delete[] ave_d;
		delete[] av2_d;
		
		
	}
	
	//sum_prog_d[i] contiene la media delle N medie della distanza quadratica valutate per il passo i-esimo
	for (int i=0; i<=n; i++){
		sum_prog_d[i] /= N;
		su2_prog_d[i] /= N;
		
		//calcolo l'errore sul passo i-esimo, considerando gli N blocchi di RW
		err_prog_d[i] = sqrt((su2_prog_d[i] - pow(sum_prog_d[i],2))/N);
		
		passi_d[i] = (int) i;
	}
	
	//scrivo su file i risultati necessari: il vettore passi, il vettore delle medie progressive ed il vettore degli errori
	outfile_3vett("RW_d.txt",n+1,passi_d,sum_prog_d,err_prog_d);
	
	delete[] sum_prog_d;
	delete[] su2_prog_d;
	delete[] passi_d;
	delete[] err_prog_d;
	
	
	
	
	double *passi_c = new double[n+1];				//numero del passo (da 0 a 100 compresi)
	double *sum_prog_c = new double[n+1]();			//contiene la somma progressiva delle medie di r_i^2 valutata su tutti gli N blocchi
	double *su2_prog_c = new double[n+1]();
	double *err_prog_c = new double[n+1]();			//contiene l'errore per il passo i-esimo valutato considerando tutti gli N blocchi
	
	
	
	//apro un ciclo sul numero di blocchi
	for (int j=0; j<N; j++){
		
		double *ave_c = new double[n+1]();			//contiene le medie per il passo i-esimo relativo al blocco j-esimo
		double *av2_c = new double[n+1]();
		
		//apro un ciclo sul numero di RW all'interno del blocco j-esimo
		for (int k=0; k<L; k++){
			double rx = 0;
			double ry = 0;
			double rz = 0;
			
			//apro un ciclo sul numero di passi da effettuare per il RW k-esimo
			for (int i=0; i<=n; i++){
				ave_c[i] += pow(rx,2) + pow(ry,2) + pow(rz,2);
				//effettuo a posteriori il passo, in modo tale da avere per ave[0] una distanza quadratica media nulla
				//in ave[n] avremo la distanza quadratica media al 100-esimo passo
				rnd.Walk_d(rx,ry,rz);
			}
			
		}
		
		//ave[i] contiene la media della distanza quadratica all'i-esimo passo per gli L RW del j-esimo blocco
		for (int i=0; i<=n; i++){
			ave_c[i] /= L;
			av2_c[i] = pow(ave_c[i],2);
			sum_prog_c[i] += ave_c[i];
			su2_prog_c[i] += av2_c[i];
		}
		
		
		delete[] ave_c;
		delete[] av2_c;
		
		
	}
	
	//sum_prog_c[i] contiene la media delle N medie della distanza quadratica valutate per il passo i-esimo
	for (int i=0; i<=n; i++){
		sum_prog_c[i] /= N;
		su2_prog_c[i] /= N;
		
		//calcolo l'errore sul passo i-esimo, considerando gli N blocchi di RW
		err_prog_c[i] = sqrt((su2_prog_c[i] - pow(sum_prog_c[i],2))/N);
		
		passi_c[i] = (int) i;
	}
	
	
	
	//scrivo su file i risultati necessari: il vettore passi, il vettore delle medie progressive ed il vettore degli errori
	outfile_3vett("RW_c.txt",n+1,passi_c,sum_prog_c,err_prog_c);
	
	delete[] sum_prog_c;
	delete[] su2_prog_c;
	delete[] passi_c;
	delete[] err_prog_c;
	
	rnd.SaveSeed();
	
	return 0;
}
