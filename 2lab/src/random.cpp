//
//  random.cpp
//  lab2
//
//  Created by Administrator on 26/03/21.
//
//

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "../include/random.hpp"

using namespace std;

Random :: Random(){}

Random :: ~Random(){}

void Random :: SaveSeed(){
	ofstream WriteSeed;
	WriteSeed.open("src/seed.out");
	if (WriteSeed.is_open()){
		WriteSeed << l1 << " " << l2 << " " << l3 << " " << l4 << endl;;
	} else cerr << "PROBLEM: Unable to open random.out" << endl;
	WriteSeed.close();
	return;
}


void Random :: Initialize(Random & rnd){
	int seed[4];
	int p1, p2;
	ifstream Primes("src/Primes");
	if (Primes.is_open()){
		Primes >> p1 >> p2 ;
	} else cerr << "PROBLEM: Unable to open Primes" << endl;
	Primes.close();
	
	ifstream input("src/seed.in");
	string property;
	if (input.is_open()){
		while ( !input.eof() ){
			input >> property;
			if( property == "RANDOMSEED" ){
				input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
				rnd.SetRandom(seed,p1,p2);
			}
		}
		input.close();
	} else cerr << "PROBLEM: Unable to open seed.in" << endl;
	
	return;
}



//genera un numero casuale distribuito secondo la distribuzione 2*(1-x)
double Random :: cos_IS(){
	double y = Rannyu();
	double x = 1-sqrt(1-y);
	return x;
}


//realizza un passo di random walk nello spazio discreto (a=1)
void Random :: Walk_d(double & rx, double & ry, double & rz){
	double appo = Rannyu(0,3);
	double appo2 = Rannyu();
	if (appo < 1){
		if (appo2 < 0.5)
			rx++;
		else
			rx--;
	}
	else if (appo >= 1 & appo <2){
		if (appo2 < 0.5)
			ry++;
		else
			ry--;
	}
	else{
		if (appo2 < 0.5)
			rz++;
		else
			rz--;
	}
	
	return;
}


//realizza un passo di random walk nello spazio continuo (a=1)
void Random :: Walk_c(double & rx, double & ry, double & rz){
	double theta = Rannyu(0,M_PI);
	double phi = Rannyu(0,2*M_PI);
	
	rx += sin(theta)*cos(phi);
	ry += sin(theta)*sin(phi);
	rz += cos(theta);
	
	return;
}



double Random :: Rannyu(double min, double max){
	return min+(max-min)*Rannyu();
}

double Random :: Rannyu(void){
	const double twom12=0.000244140625;
	int i1,i2,i3,i4;
	double r;
	
	i1 = l1*m4 + l2*m3 + l3*m2 + l4*m1 + n1;
	i2 = l2*m4 + l3*m3 + l4*m2 + n2;
	i3 = l3*m4 + l4*m3 + n3;
	i4 = l4*m4 + n4;
	l4 = i4%4096;
	i3 = i3 + i4/4096;
	l3 = i3%4096;
	i2 = i2 + i3/4096;
	l2 = i2%4096;
	l1 = (i1 + i2/4096)%4096;
	r=twom12*(l1+twom12*(l2+twom12*(l3+twom12*(l4))));
	
	return r;
}

void Random :: SetRandom(int * s, int p1, int p2){
	m1 = 502;
	m2 = 1521;
	m3 = 4071;
	m4 = 2107;
	l1 = s[0];
	l2 = s[1];
	l3 = s[2];
	l4 = s[3];
	n1 = 0;
	n2 = 0;
	n3 = p1;
	n4 = p2;
	
	return;
}
