//
//  lib.cpp
//  lab2
//
//  Created by Administrator on 26/03/21.
//
//

#include "../include/lib.hpp"

using namespace std;

string output_dir = "output/";


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

void block_ave(string fileName, int N, double *vett, double *err_prog){
	double *vett2 = new double[N]();
	double *sum_prog = new double[N]();
	double *su2_prog = new double[N]();
	
	
	for (int i=0; i<N; i++){
		vett2[i]=pow(vett[i],2);
	}
	
	sumprog(N,sum_prog,vett);
	sumprog(N,su2_prog,vett2);
	
	
	//valuto l'errore sui valori appena calcolati
	for (int i=0; i<N; i++)
		err_prog[i] = error(sum_prog,su2_prog,i);
	
	
	ofstream file_out;
	
	file_out.open(output_dir+fileName);
	
	if (file_out.is_open()){
		for (int i=0; i<N; i++){	
			//carico sul file i tre vettori 
			file_out << (int) (i+1) << "    " << sum_prog[i] << "    " << err_prog[i] << endl;
		}
		
		file_out.close();
		file_out.clear();
	}
	else cout << "Unable to open output file" << endl;
	
	delete[] vett2;
	delete[] sum_prog;
	delete[] su2_prog;
}


//valuto la funzione pi/2*cos(x*pi/2)
double eval_cos(double x){
	double gx;
	gx = M_PI_2*cos(M_PI_2*x);
	
	return gx;
}


//valuto la funzione ricavata per IS di cos(x*pi/2)*pi/2
//che sarà pi*cos(pi/2*x)/(4*(1-x)), dove x è un numero random
//estratto dalla disribuzione di probabilità 2(1-x)
double eval_cos_IS(double x){
	double gx;
	gx = M_PI_4*cos(M_PI_2*x)/(1-x);
	
	return gx;
}



//carico su file un vettore double
void outfile_vett(string fileName, int N, double *vett){
	ofstream file_out;
	
	file_out.open(output_dir+fileName);
	
	
	if (file_out.is_open())
	{
		for (int i=0; i<N; i++){	
			//carico sul file il vettore
			file_out << vett[i] << endl;
		}
		
		file_out.close();
		file_out.clear();
	}
	else cout << "Unable to open output file" << endl;
}

//carico su file tre vettori di double 
void outfile_3vett(string fileName, int N, double *vett1, double *vett2, double *vett3){
	ofstream file_out;
	
	file_out.open(output_dir+fileName);
	
	
	if (file_out.is_open())
	{
		for (int i=0; i<N; i++){	
			//carico sul file i tre vettori 
			file_out << vett1[i] << "    " << vett2[i] << "    " << vett3[i] << endl;
		}
		
		file_out.close();
		file_out.clear();
	}
	else cout << "Unable to open output file" << endl;
}

