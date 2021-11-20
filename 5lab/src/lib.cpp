//
//  lib.cpp
//  4lab
//
//  Created by Administrator on 08/07/21.
//
//

#include "../include/lib.hpp"

using namespace std;



// errore statistico
double error(double * AV, double * AV2, int n){
	if(n==0){
		return 0;
	}
	else{
		return sqrt((AV2[n] - pow(AV[n],2))/n);
	}
}

//effettuo una somma progressiva
void sumprog(int N, double *sum_prog, double *vett){
	for (int i=0; i<N; i++){
		for (int j=0; j<i+1; j++){
			sum_prog[i] += vett[j];
		}
		
		sum_prog[i] /= (i+1);
	}
}

void block_ave(string fileName, int N, int L, double *vett){
	double *vett2 = new double[N]();
	double *sum_prog = new double[N]();
	double *su2_prog = new double[N]();
	double *err_prog = new double[N]();
	
	for (int i=0; i<N; i++){
		vett[i] /= L;
		vett2[i] = pow(vett[i],2);
	}
	
	sumprog(N,sum_prog,vett);
	sumprog(N,su2_prog,vett2);
	
	
	//valuto l'errore sui valori appena calcolati
	for (int i=0; i<N; i++)
		err_prog[i] = error(sum_prog,su2_prog,i);
	
	ofstream file_out;
	
	file_out.open(fileName);
	
	if (file_out.is_open()){
		for (int i=0; i<N; i++){	
			//carico sul file i tre vettori 
			file_out << (int) (i+1)*L << "    " << sum_prog[i] << "    " << err_prog[i] << endl;
		}
		
		file_out.close();
		file_out.clear();
	}
	else cout << "Unable to open output file" << endl;
	
	delete[] vett2;
	delete[] sum_prog;
	delete[] su2_prog;
	delete[] err_prog;

}
