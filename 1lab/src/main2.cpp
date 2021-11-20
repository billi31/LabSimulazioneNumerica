//
//  2.cpp
//  lab1
//
//  Created by Administrator on 25/03/21.
//
//

#include "../include/lib.hpp"
#include "../include/random.hpp"

using namespace std;


int main (int argc, char *argv[]){
	
	int L = 10000;						//numero di realizzazioni
	
	string output_dir = "output/";

	
	//definisco una variabile Random che utilizzerò per generare i numeri casuali
	Random rnd;
	rnd.Initialize(rnd);
	
	//definisco il numero di punti che si sommano di volta in volta
	int N_sum[4] = {1,2,10,100};
	
	
	
	
	
	//DISTRIBUZIONE UNIFORME
	
	//creo una matrice che conterrà tutti gli L valori della somma di (1,2,10,100) termini
	//per la distribuzione uniforme
	double (*sN_u)[4] = new double[L][4]();
	
	//creo un ciclo per i diversi valori di N_sum considerati
	for (int i=0; i<4; i++){
		
		for (int k=0; k<L; k++){
			//genero N_sum[i] numeri casuali
			for (int j=0; j<N_sum[i]; j++){
				double appo = rnd.Rannyu();
				sN_u[k][i] += appo;
			}
			
			sN_u[k][i] /= N_sum[i];
			
		}
		
	}
	
	ofstream file_out;
	file_out.open(output_dir + "uniform.txt");
	
	if (file_out.is_open()){
		for (int k=0; k<L; k++){
			//esporto i valori delle L somme di (1,2,10,100) numeri distribuiti uniformemente in un file di testo
			file_out << sN_u[k][0] << "    " << sN_u[k][1] << "    " << sN_u[k][2] << "    " << sN_u[k][3] << endl;
		}
		
		file_out.close();
		file_out.clear();
		
	}
	else cout << "Unable to open output file" << endl;
	
	
	delete [] sN_u;
	
	
	
	
	
	//DISTRIBUZIONE ESPONENZIALE
	
	//creo una matrice che conterrà tutti gli L valori della somma di (1,2,10,100) termini
	//per la distribuzione esponenziale (lambda=1) 
	double (*sN_exp)[4] = new double[L][4]();
	
	//creo un ciclo per i diversi valori di N_sum considerati
	for (int i=0; i<4; i++){
		
		for (int k=0; k<L; k++){
			//genero N_sum[i] numeri casuali
			for (int j=0; j<N_sum[i]; j++){
				double appo = rnd.exponential(1.);
				sN_exp[k][i] += appo;
			}
			
			sN_exp[k][i] /= N_sum[i];
			
		}
		
	}
	
	file_out.open(output_dir + "exponential.txt");
	
	if (file_out.is_open()){
		for (int k=0; k<L; k++){
			//esporto i valori delle L somme di (1,2,10,100) numeri distribuiti esponenzialmente in un file di testo
			file_out << sN_exp[k][0] << "    " << sN_exp[k][1] << "    " << sN_exp[k][2] << "    " << sN_exp[k][3] << endl;
		}
		
		file_out.close();
		file_out.clear();
		
	}
	else cout << "Unable to open output file" << endl;
	
	
	delete [] sN_exp;
	
	
	
	
	
	
	//DISTRIBUZIONE LORENTZIANA
	
	//creo una matrice che conterrà tutti gli L valori della somma di (1,2,10,100) termini
	//per la distribuzione di Lorentz (mu=0, Gamma=1)
	double (*sN_CL)[4] = new double[L][4]();
	
	//creo un ciclo per i diversi valori di N_sum considerati
	for (int i=0; i<4; i++){
		
		for (int k=0; k<L; k++){
			//genero N_sum[i] numeri casuali
			for (int j=0; j<N_sum[i]; j++){
				double appo = rnd.CauchyLorentz(0,1.);
				sN_CL[k][i] += appo;
			}
			
			sN_CL[k][i] /= N_sum[i];
			
		}
		
	}
	
	
	
	file_out.open(output_dir + "lorentzian.txt");
	
	if (file_out.is_open()){
		for (int k=0; k<L; k++){
			//esporto i valori delle L somme di (1,2,10,100) numeri distribuiti secondo una distribuzione di Lorentz in un file di testo
			file_out << sN_CL[k][0] << "    " << sN_CL[k][1] << "    " << sN_CL[k][2] << "    " << sN_CL[k][3] << endl;
		}
		
		file_out.close();
		file_out.clear();
		
	}
	else cout << "Unable to open output file" << endl;
	
	
	delete [] sN_CL;
	
	rnd.SaveSeed();
	
	return 0;
}
