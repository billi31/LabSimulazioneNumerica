//
//  MolDyn_NVE.cpp
//  4lab
//
//  Created by Administrator on 22/04/21.
//
//

#include "../include/MolDyn_NVE.hpp"

using namespace std;


MolDynamics::MolDynamics(){
}


MolDynamics::~MolDynamics(){
}


void MolDynamics::Input(const char *input_file, const char *config0){ //Prepare all stuff for the simulation
	ifstream ReadInput,ReadConf;
	
	Random rnd;
	rnd.Initialize(rnd);	
	
	cout << "Classic Lennard-Jones fluid        " << endl;
	cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
	cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
	cout << "The program uses Lennard-Jones units " << endl;
	
	
	ReadInput.open(input_file); //Read input
	if(ReadInput.fail()){
		cerr << "Unable to open " << input_file << endl ;
		exit(1);
	}
	cout << "Reading setup from file " << input_file << endl << endl;
	
	ReadInput >> temp;
	
	ReadInput >> npart;
	cout << "Number of particles = " << npart << endl;
	
	ReadInput >> rho;
	cout << "Density of particles = " << rho << endl;
	vol = (double)npart/rho;
	cout << "Volume of the simulation box = " << vol << endl;
	box = pow(vol,1.0/3.0);
	cout << "Edge of the simulation box = " << box << endl;
	
	ReadInput >> rcut;
	ReadInput >> delta;
	ReadInput >> nstep;
	ReadInput >> iprint;
	ReadInput >> nequil;
	ReadInput >> nequil_step;
	
	cout << "The program integrates Newton equations with the Verlet method " << endl;
	cout << "Time step = " << delta << endl;
	cout << "Number of steps = " << nstep << endl << endl;
	ReadInput.close();
	
	//Prepare array for measurements
	iv = 0; //Potential energy
	ik = 1; //Kinetic energy
	ie = 2; //Total energy
	it = 3; //Temperature
	ip = 4;	//Pressure
	n_props = 5; //Number of observables
	
	
	//Read initial configuration
	ReadConf.open(config0);
	if(ReadConf.fail()){
		cerr << "Unable to open " << config0 << endl ;
		exit(1);
	}
	cout << "Read initial configuration from file " << config0 << endl << endl;
	
	for (int i=0; i<npart; ++i){
		ReadConf >> x[i] >> y[i] >> z[i];
		x[i] = x[i] * box;
		y[i] = y[i] * box;
		z[i] = z[i] * box;
	}
	ReadConf.close();
	
	//Prepare initial velocities
	cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
	double sumv[3] = {0.0, 0.0, 0.0};
	for (int i=0; i<npart; ++i){
		vx[i] = rnd.Rannyu(-0.5,0.5);			//we now have random velocities (-0.5;0.5]
		vy[i] = rnd.Rannyu(-0.5,0.5);
		vz[i] = rnd.Rannyu(-0.5,0.5);
		
		sumv[0] += vx[i];
		sumv[1] += vy[i];
		sumv[2] += vz[i];
	}
	for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
	double sumv2 = 0.0, fs;
	for (int i=0; i<npart; ++i){				//velocities are centered in 0
		vx[i] = vx[i] - sumv[0];
		vy[i] = vy[i] - sumv[1];
		vz[i] = vz[i] - sumv[2];
		
		sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
	}
	sumv2 /= (double)npart;
	
	
	fs = sqrt(3 * temp / sumv2);				//fs = velocity scale factor
	
	cout << "Rescale velocities to match the desired temperature\n";
	cout << "Scaling factor = " << fs << "\n\n";
	
	for (int i=0; i<npart; ++i){				//we now have velocities centered in 0 in an interval that matches the temperature requested
		vx[i] *= fs;
		vy[i] *= fs;
		vz[i] *= fs;
		
		xold[i] = Pbc(x[i] - vx[i] * delta);	//evaluation of positions at time t
		yold[i] = Pbc(y[i] - vy[i] * delta);
		zold[i] = Pbc(z[i] - vz[i] * delta);
	}
	
	Equilibrate(input_file);
	
	rnd.SaveSeed();
	
	return;
}

void MolDynamics::Input(const char *input_file, const char *config0, const char *old0){ //Prepare all stuff for the simulation
	ifstream ReadInput,ReadConf;
	
	cout << "Classic Lennard-Jones fluid        " << endl;
	cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
	cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
	cout << "The program uses Lennard-Jones units " << endl;
	
	//seed = 1;    //Set seed for random numbers
	//srand(seed); //Initialize random number generator
	
	ReadInput.open(input_file); //Read input
	if(ReadInput.fail()){
		cerr << "Unable to open " << input_file << endl ;
		exit(1);
	}
	cout << "Reading setup from file " << input_file << endl << endl;
	
	
	ReadInput >> temp;
	
	ReadInput >> npart;
	cout << "Number of particles = " << npart << endl;
	
	ReadInput >> rho;
	cout << "Density of particles = " << rho << endl;
	vol = (double)npart/rho;
	cout << "Volume of the simulation box = " << vol << endl;
	box = pow(vol,1.0/3.0);
	cout << "Edge of the simulation box = " << box << endl;
	
	ReadInput >> rcut;
	ReadInput >> delta;
	ReadInput >> nstep;
	ReadInput >> iprint;
	ReadInput >> nequil;
	ReadInput >> nequil_step;
	
	cout << "The program integrates Newton equations with the Verlet method " << endl;
	cout << "Time step = " << delta << endl;
	cout << "Number of steps = " << nstep << endl << endl;
	ReadInput.close();
	
	//Prepare array for measurements
	iv = 0; //Potential energy
	ik = 1; //Kinetic energy
	ie = 2; //Total energy
	it = 3; //Temperature
	ip = 4; //Pressure
	n_props = 5; //Number of observables
	
	//Read initial configuration
	ReadConf.open(old0);
	if(ReadConf.fail()){
		cerr << "Unable to open " << old0 << endl ;
		exit(1);
	}
	cout << "Read initial old configuration from file " << old0 << endl << endl;
	
	for (int i=0; i<npart; ++i){
		ReadConf >> xold[i] >> yold[i] >> zold[i];
		xold[i] = xold[i] * box;
		yold[i] = yold[i] * box;
		zold[i] = zold[i] * box;
	}
	ReadConf.close();
	
	
	//Read initial configuration
	ReadConf.open(config0);
	if(ReadConf.fail()){
		cerr << "Unable to open " << config0 << endl ;
		exit(1);
	}
	cout << "Read initial configuration from file " << config0 << endl << endl;
	
	for (int i=0; i<npart; ++i){
		ReadConf >> x[i] >> y[i] >> z[i];
		x[i] = x[i] * box;
		y[i] = y[i] * box;
		z[i] = z[i] * box;
	}
	ReadConf.close();
	
	Equilibrate(input_file);
	
	return;
}

void MolDynamics::Reset(void){ //Reset block averages
	
	if(iblk == 1){
		for(int i=0; i<n_props; ++i){
			glob_av[i] = 0;
			glob_av2[i] = 0;
		}
	}
	
	for(int i=0; i<n_props; ++i)
		blk_av[i] = 0;
	
	blk_norm = 0;
}


void MolDynamics::Accumulate(void){ //Update block averages
	
	for(int i=0; i<n_props; ++i)
		blk_av[i] += inst_measure[i];
	
	blk_norm = blk_norm + 1.0;
}


void MolDynamics::Averages(void){ //Print results for current block
	
	ofstream Epot, Ekin, Etot, Temp, Pres;
	
	const int wd = 12;
	cout << "Block number " << iblk << endl;
	
	Epot.open("output.epot.0",ios::app);
	stima_epot = blk_av[iv]/blk_norm/(double)npart; //Potential Energy
	glob_av[iv]  += stima_epot;
	glob_av2[iv] += stima_epot*stima_epot;
	err_v = Error(glob_av[iv],glob_av2[iv]);
	Epot << setw(wd) << iblk <<  setw(wd) << stima_epot << setw(wd) << glob_av[iv]/(double)iblk << setw(wd) << err_v << endl;
	Epot.close();
	
	Ekin.open("output.ekin.0",ios::app);
	stima_ekin = 0.5*blk_av[ik]/blk_norm/(double)npart; //Kinetic Energy
	glob_av[ik]  += stima_ekin;
	glob_av2[ik] += stima_ekin*stima_ekin;
	err_k = Error(glob_av[ik],glob_av2[ik]);
	Ekin << setw(wd) << iblk <<  setw(wd) << stima_ekin << setw(wd) << glob_av[ik]/(double)iblk << setw(wd) << err_k << endl;
	Ekin.close();
	
	Etot.open("output.etot.0",ios::app);
	stima_etot = blk_av[ie]/blk_norm/(double)npart; //Total Energy
	glob_av[ie]  += stima_etot;
	glob_av2[ie] += stima_etot*stima_etot;
	err_e = Error(glob_av[ie],glob_av2[ie]);
	Etot << setw(wd) << iblk <<  setw(wd) << stima_etot << setw(wd) << glob_av[ie]/(double)iblk << setw(wd) << err_e << endl;
	Etot.close();
	
	Temp.open("output.temp.0",ios::app);
	stima_temp = blk_av[it]/3./blk_norm/(double)npart; //Temperature
	glob_av[it]  += stima_temp;
	glob_av2[it] += stima_temp*stima_temp;
	err_t = Error(glob_av[it],glob_av2[it]);
	Temp << setw(wd) << iblk <<  setw(wd) << stima_temp << setw(wd) << glob_av[it]/(double)iblk << setw(wd) << err_t << endl;
	Temp.close();
	
	Pres.open("output.pres.0",ios::app);
	stima_pres = blk_av[ip]/3./vol/blk_norm; //Pressure
	glob_av[ip]  += stima_pres;
	glob_av2[ip] += stima_pres*stima_pres;
	err_p = Error(glob_av[ip],glob_av2[ip]);
	Pres << setw(wd) << iblk <<  setw(wd) << stima_pres << setw(wd) << glob_av[ip]/(double)iblk << setw(wd) << err_p << endl;
	Pres.close();
	
	
	cout << "----------------------------" << endl << endl;
}


void MolDynamics::Rescale(void){ //Rescaling velocities
	
	double sumv2 = 0.0;
	double xnew[m_part], ynew[m_part], znew[m_part], fx[m_part], fy[m_part], fz[m_part];
	
	for(int i=0; i<npart; ++i){ //Force acting on particle i
		fx[i] = Force(i,0);
		fy[i] = Force(i,1);
		fz[i] = Force(i,2);
	}
	
	
	
	//Making the first step of Verlet, in order to obtain the positions at time t+dt and the velocities at time t
	for(int i=0; i<npart; ++i){
		
		xnew[i] = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );		//positions at time t+dt
		ynew[i] = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
		znew[i] = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );
		
		vx[i] = Pbc(xnew[i] - x[i])/delta;						//velocities at time t
		vy[i] = Pbc(ynew[i] - y[i])/delta;
		vz[i] = Pbc(znew[i] - z[i])/delta;
		
		
		sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
		
	}
	
	sumv2 /= (double)npart;
	
	//double old_temp;
	//old_temp = sumv2 / 3.0;			//Temperature at time t before rescaling
	
	
	double fs;
	fs = sqrt(3 * temp / sumv2);
	
	
	cout << "Rescale velocities to match the desired temperature\n";
	cout << "Scaling factor = " << fs << "\n\n";
	
	for (int i=0; i<npart; ++i){
		vx[i] *= fs;								//corrected velocities at time t
		vy[i] *= fs;
		vz[i] *= fs;
		
		x[i] = Pbc(xnew[i] - vx[i] * delta);		
		y[i] = Pbc(ynew[i] - vy[i] * delta);
		z[i] = Pbc(znew[i] - vz[i] * delta);
		
		xold[i] = x[i];								//new initial positions at time t and t+dt
		yold[i] = y[i];
		zold[i] = z[i];
		
		x[i] = xnew[i];
		y[i] = ynew[i];
		z[i] = znew[i];
		
	}
	
	return;
}

void MolDynamics::Move(void){ //Move particles with Verlet algorithm
	double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];
	
	for(int i=0; i<npart; ++i){ //Force acting on particle i
		fx[i] = Force(i,0);
		fy[i] = Force(i,1);
		fz[i] = Force(i,2);
	}
	
	for(int i=0; i<npart; ++i){ //Verlet integration scheme
		
		xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );		//positions at time t+dt
		ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
		znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );
		
		vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);						//velocities at time t
		vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
		vz[i] = Pbc(znew - zold[i])/(2.0 * delta);
		
		xold[i] = x[i];													//positions at time t
		yold[i] = y[i];
		zold[i] = z[i];
		
		x[i] = xnew;
		y[i] = ynew;
		z[i] = znew;
	}
	return;
}

double MolDynamics::Force(int ipar, int idir){ //Compute forces as -Grad_ipar V(r)
	double f=0.0;
	double dvec[3], dr;
	
	for (int i=0; i<npart; ++i){
		if(i != ipar){
			dvec[0] = Pbc( x[ipar] - x[i] );  // distance ipar-i in pbc
			dvec[1] = Pbc( y[ipar] - y[i] );
			dvec[2] = Pbc( z[ipar] - z[i] );
			
			dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
			dr = sqrt(dr);
			
			if(dr < rcut){
				f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_iparV(r)
			}
		}
	}
	
	return f;
}

void MolDynamics::Equilibrate(const char* input_file){ //Equilibrate the system
	
	//Possibility to rescale the velocity
	string ans;
	do {
		cout << "Would you like to equilibrate the system in order to match the temperature in file " << input_file << "? ( y / n )\n";
		cin >> ans;
	} while (ans != "y" && ans != "n");
	
	if (ans == "y"){
		for (int i=0; i<nequil; i++){
			Rescale();
			for (int j=0; j<nequil_step; j++){
				Move();
				Measure();
				PrintMeasure();
			}
			
		}
	}
	
	return;
}

void MolDynamics::Measure(void){ //Properties measurement
	double v, t, vij, p;
	double dx, dy, dz, dr;
	
	
	v = 0.0; //reset observables
	t = 0.0;
	p = 0.0;
	
	//cycle over pairs of particles
	for (int i=0; i<npart-1; ++i){
		for (int j=i+1; j<npart; ++j){
			
			dx = Pbc( xold[i] - xold[j] ); // here I use old configurations [old = r(t)]
			dy = Pbc( yold[i] - yold[j] ); // to be compatible with EKin which uses v(t)
			dz = Pbc( zold[i] - zold[j] ); // => EPot should be computed with r(t)
			
			dr = dx*dx + dy*dy + dz*dz;
			dr = sqrt(dr);
			
			if(dr < rcut){
				vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
				p += 1./pow(dr,12) - 0.5/pow(dr,6);
				
				//Potential energy
				v += vij;
			}
		}
	}
	
	//Kinetic energy
	for (int i=0; i<npart; ++i) t += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
	
	inst_measure[iv] = v;
	inst_measure[ik] = t;
	inst_measure[ie] = v + 0.5*t;
	inst_measure[it] = t;
	inst_measure[ip] = t + 48.*p;
	
	return;
}


void MolDynamics::PrintMeasure(void){
	ofstream Epot, Ekin, Etot, Temp, Pres;
	
	Epot.open("insta.epot.dat",ios::app);
	Ekin.open("insta.ekin.dat",ios::app);
	Etot.open("insta.etot.dat",ios::app);
	Temp.open("insta.temp.dat",ios::app);
	Pres.open("insta.pres.dat",ios::app);
	
	stima_epot = inst_measure[iv]/(double)npart;
	stima_ekin = 0.5*inst_measure[ik]/(double)npart;
	stima_etot = inst_measure[ie]/(double)npart;
	stima_temp = inst_measure[it]/3./(double)npart;
	stima_pres = inst_measure[ip]/3./vol;
	
	Epot << stima_epot << endl;
	Ekin << stima_ekin << endl;
	Etot << stima_etot << endl;
	Temp << stima_temp << endl;
	Pres << stima_pres << endl;
		
	Epot.close();
	Ekin.close();
	Etot.close();
	Temp.close();
	Pres.close();
}


void MolDynamics::ConfOld(void){ //Write old configuration
	ofstream WriteConf;
	
	cout << "Print old configuration to file old.final " << endl;
	WriteConf.open("old.final");
	
	for (int i=0; i<npart; ++i){
		WriteConf << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
	}
	WriteConf.close();
	return;
}



void MolDynamics::ConfFinal(void){ //Write final configuration
	ofstream WriteConf;
	
	cout << "Print final configuration to file config.final " << endl << endl;
	WriteConf.open("config.final");
	
	for (int i=0; i<npart; ++i){
		WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
	}
	WriteConf.close();
	return;
}

void MolDynamics::ConfXYZ(int nconf){ //Write configuration in .xyz format
	ofstream WriteXYZ;
	
	WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
	WriteXYZ << npart << endl;
	WriteXYZ << "This is only a comment!" << endl;
	for (int i=0; i<npart; ++i){
		WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
	}
	WriteXYZ.close();
}

double MolDynamics::Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
	return r - box * rint(r/box);
}

double MolDynamics::Error(double sum, double sum2){
	if (iblk==1)
		return 0.0;
	else
		return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}
