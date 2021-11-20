//
//  lib.hpp
//  4lab
//
//  Created by Administrator on 08/07/21.
//
//

#ifndef lib_hpp
#define lib_hpp

#include <stdio.h>
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow


double error(double * AV, double * AV2, int n);

void sumprog(int N, double *sum_prog, double *vett);

void block_ave(std::string fileName, int N, int L, double *vett);


#endif /* lib_hpp */
