//
//  lib.hpp
//  lab2
//
//  Created by Administrator on 26/03/21.
//
//

#ifndef lib_hpp
#define lib_hpp

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include "../include/random.hpp"


//errore statistico

double error(double *, double *, int);

//valuta la media progressiva di un vettore vett di dimensione N

void sumprog(int N, double *sum_prog, double *vett);


void block_ave(std::string fileName, int N, double *vett, double *err_prog);

//valuto la cumulativa di una gaussiana di media nulla e varianza pari a 1

double N(double);

// valuta l'evoluzione del prezzo al tempo t

double price(double t, double S0, double mu, double sigma, Random &rnd);

//evoluzione discreta del prezzo da un tempo t1 ad un tempo t2

double price_discr(double t1, double t2, double S1, double mu, double sigma, Random &rnd);

//carica su file 1 vettore double

void outfile_vett(std::string, int, double *);

//carica su file 3 vettori double

void outfile_3vett(std::string, int, double *, double *, double *);


#endif /* lib_hpp */
