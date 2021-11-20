//
//  Hydrogen.cpp
//  5lab
//
//  Created by Administrator on 03/08/21.
//

#include "../include/Hydrogen.hpp"

using namespace std;

string output_directory = "output/";

Hydrogen::Hydrogen(){
	state = "null";
	mode = "null";
	rnd.Initialize(rnd);
}

Hydrogen::~Hydrogen(){
}

void Hydrogen::Reset(void){ //Reset block averages
	
	if(iblk == 1){
		glob_av = 0;
		glob_av2 = 0;
	}
	
	blk_av = 0;
	blk_norm = 0;
	
}


void Hydrogen::Initialize(double x, double y, double z, string mod, string stat){
	state = stat;
	mode = mod;
	SetPos(x,y,z);
	
	attempted = 0;
	accepted = 0;
	
	if (mod == "unif"){
		if (stat == "100")
			SetEdge(2.43);
		else if (stat == "210")
			SetEdge(5.95);
		else{
			cerr << "Invalid state. Must be 100 or 210.\n\n";
			exit(1);
		}
	}
	else if (mod == "gauss"){
		if (stat == "100")
			SetEdge(1.5);
		else if (stat == "210")
			SetEdge(3.74);
		else{
			cerr << "Invalid state. Must be 100 or 210.\n\n";
			exit(1);
		}
	}
	else{
		cerr << "Invalid mode. Must be unif or gauss.\n\n";
		exit(1);
	}
	
	cout << "\nSimulation Hydrogen Wave Function\n\n";
	cout << "state = " << state << "\nmode = " << mode << "\n\n";
	
	Equilibration();
	
	return;
}

void Hydrogen::Radius(void){
	radius = sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
	
	return;
}


void Hydrogen::PrintPos(void){
	ofstream out;
	out.open(output_directory + "graph." + state + mode + ".dat",ios::app);
	
	out << pos[0] << "   " <<  pos[1] << "   " << pos[2] << endl;
	
	out.close();
	
	return;
}


void Hydrogen::Alpha(double xnew, double ynew, double znew){
	double rnew = sqrt(xnew*xnew + ynew*ynew + znew*znew);
	
	if (state == "100")
		alpha = min(1.,exp(-2.*rnew + 2.*radius));
	
	else if (state == "210")
		alpha = min(1.,exp(-rnew + radius)*(znew*znew)/(pos[2]*pos[2]));
	
	else{
		cerr << "Invalid state. Must be 100 or 210.\n\n";
		exit(1);
	}
	
	return;
}
	
	
double Hydrogen::Acceptance(void){
	double acc = (double) accepted/attempted;
	
	return acc;
}

void Hydrogen::Move(){
	double xnew, ynew, znew;
	
	if (mode == "unif"){
		xnew = pos[0] + edge * rnd.Rannyu(-0.5,0.5); 
		ynew = pos[1] + edge * rnd.Rannyu(-0.5,0.5); 
		znew = pos[2] + edge * rnd.Rannyu(-0.5,0.5);
	}
	else if (mode == "gauss"){
		xnew = pos[0] + edge * rnd.Gauss(0,0.5); 
		ynew = pos[1] + edge * rnd.Gauss(0,0.5); 
		znew = pos[2] + edge * rnd.Gauss(0,0.5);
	}
	else {
		cerr << "Invalid mode. Must be unif or gauss.\n\n";
		exit(1);
	}
	
	Radius();
	Alpha(xnew,ynew,znew);
	
	if (rnd.Rannyu() <= alpha) {
		pos[0] = xnew;
		pos[1] = ynew;
		pos[2] = znew;
		accepted++;
	}
	
	attempted++;
	
	return;
}

void Hydrogen::Equilibration(void){
	for (int i=0; i<nequil; i++)
		Move();
	
	return;
}
	

void Hydrogen::Accumulate(void){ //Update block averages
	
	blk_av += radius;
	blk_norm = blk_norm + 1.0;
}


void Hydrogen::Averages(void){ //Print results for current block
	
	ofstream R;
	
	const int wd = 12;
	cout << "Block number " << iblk << endl;
	
	R.open(output_directory + "output.r" + state + mode + ".0",ios::app);
	stima_r = blk_av/blk_norm; //Potential Energy
	glob_av  += stima_r;
	glob_av2 += stima_r*stima_r;
	err_r = Error(glob_av,glob_av2);
	R << setw(wd) << iblk <<  setw(wd) << stima_r << setw(wd) << glob_av/(double)iblk << setw(wd) << err_r << endl;
	R.close();
	
	cout << "----------------------------" << endl << endl;
}

void Hydrogen::PrintInstantaneusRadius(void){
	ofstream R;
	
	R.open(output_directory + "insta.r" + state + mode + ".dat",ios::app);
	
	R << radius << endl;
	
	R.close();
}


double Hydrogen::Error(double sum, double sum2){
	if (iblk==1)
		return 0.0;
	else
		return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}


