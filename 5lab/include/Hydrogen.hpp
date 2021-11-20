//
//  Hydrogen.hpp
//  5lab
//
//  Created by Administrator on 27/07/21.
//
//

#ifndef Hydrogen_hpp
#define Hydrogen_hpp

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cmath>        // rint, pow
#include <iomanip>
#include "random.hpp"


//const int ntrial = 1E6;
const int nblock = 100;
const int nstep = 1E6;
const int nequil = 150;

class Hydrogen{
public:
	
	Hydrogen();
	~Hydrogen();
	
	void Initialize(double x, double y, double z, std::string mod, std::string stat);
	void SetEdge(double edg){edge = edg; return;};
	void SetPos(double x, double y, double z){pos[0] = x; pos[1] = y; pos[2] = z; return;};
	void Radius(void);
	void PrintPos(void);
	void Reset(void);
	void Alpha(double xnew, double ynew, double znew);
	double Acceptance(void);
	void Move();
	void Equilibration(void);
	void Accumulate(void);
	void Averages(void);
	void Measure(void);
	void PrintInstantaneusRadius(void);
	double Error(double sum, double sum2);
	
	
	double r[nblock];
	
	int iblk;
	int accepted,attempted;	
	
	double radius;
	
	std::string mode;
	std::string state;
	
	Random rnd;
	
private:
	double edge;
	double pos[3];
	double alpha;
	
	double stima_r;
	double err_r;
	
	double blk_av,blk_norm;
	double glob_av,glob_av2;
	
};


#endif /* Hydrogen_hpp */
