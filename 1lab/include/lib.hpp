//
//  lib.hpp
//  lab1
//
//  Created by Administrator on 18/03/21.
//
//

#ifndef lib_hpp
#define lib_hpp

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>


//errore statistico

double error(double *, double *, int);

//valuta la media progressiva di un vettore vett di dimensione N

void sumprog(int N, double *sum_prog, double *vett);

void block_ave(std::string fileName, int N, double *vett, double *err_prog);


//carica su file 1 vettore double

void outfile_vett(std::string, int, double *);

//carica su file 3 vettori double

void outfile_3vett(std::string, int, double *, double *, double *);


#endif /* lib_hpp */
