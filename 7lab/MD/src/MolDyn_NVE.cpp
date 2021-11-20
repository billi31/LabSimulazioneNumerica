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


void MolDynamics::Input(void){ //Prepare all stuff for the simulation
	ifstream ReadInput,ReadConf;
	
	
	cout << "Classic Lennard-Jones fluid        " << endl;
	cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
	cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
	cout << "The program uses Lennard-Jones units " << endl;
	
	
	ReadInput.open("input.dat"); //Read input
	if(ReadInput.fail()){
		cerr << "Unable to open input.dat" << endl ;
		exit(1);
	}
	cout << "Reading setup from file input.dat" << endl << endl;
	
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
	ReadInput >> nblk;
	ReadInput >> iprint;
	ReadInput >> nequil;
	ReadInput >> nequil_step;
	
	cout << "The program integrates Newton equations with the Verlet method " << endl;
	cout << "Time step = " << delta << endl;
	cout << "Number of blocks = " << nblk << endl;
	cout << "Number of steps per block = " << nstep << endl;
	cout << "Number of equilibration step = " << nequil << endl;
	cout << "Number of step per equilibration = " << nequil_step << endl << endl;
	ReadInput.close();
	
	//Prepare array for measurements
	iv = 0; //Potential energy
	ik = 1; //Kinetic energy
	ie = 2; //Total energy
	it = 3; //Temperature
	ip = 4; //Pressure
	
	n_props = 5; //Number of observables
	
	igofr = 5;
	nbins = 100;
	n_props = n_props + nbins;
	bin_size = (box/2.0)/(double)nbins;
	
	
	
	//Read initial configuration
	ReadConf.open("config.0");
	if(ReadConf.fail()){
		cerr << "Unable to open config.0" << endl ;
		exit(1);
	}
	cout << "Read initial configuration from file config.0" << endl << endl;

	for (int i=0; i<npart; ++i){
		ReadConf >> x[i] >> y[i] >> z[i];
		x[i] = x[i] * box;
		y[i] = y[i] * box;
		z[i] = z[i] * box;
	}
	ReadConf.close();
	
	
	if (mode == "start"){
		Random rnd;
		rnd.Initialize(rnd);	
		
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
	}
	else if (mode == "repeat"){
		
		//Read initial configuration
		ReadConf.open("old.0");
		if(ReadConf.fail()){
			cerr << "Unable to open old.0" << endl ;
			exit(1);
		}
		
		cout << "Read initial old configuration from file old.0" << endl << endl;
		
		for (int i=0; i<npart; ++i){
			ReadConf >> xold[i] >> yold[i] >> zold[i];
			xold[i] = xold[i] * box;
			yold[i] = yold[i] * box;
			zold[i] = zold[i] * box;
		}
		ReadConf.close();
		
	}
	
	Equilibrate();
	
	return;
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
		
		vx[i] = Pbc(xnew[i] - x[i])/delta;						//velocities at time t+dt/2
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
		
		//use $\vec{v}_s(t)$ to estimate a novel old spatial configuration: 
		//$\vec{r}_{new}(t) = \vec{r}(t+dt) - dt \vec{v}_s$
		//use $\vec{r}_{new}(t)$ and $\vec{r}(t+dt)$ to start the simulation
		
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

double MolDynamics::Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
	double f=0.0;
	double dvec[3], dr;
	
	for (int i=0; i<npart; ++i){
		if(i != ip){
			dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
			dvec[1] = Pbc( y[ip] - y[i] );
			dvec[2] = Pbc( z[ip] - z[i] );
			
			dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
			dr = sqrt(dr);
			
			if(dr < rcut){
				f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
			}
		}
	}
	
	return f;
}

void MolDynamics::Equilibrate(void){ //Equilibrate the system
	
	//Possibility to rescale the velocity
	string ans;
	do {
		cout << "Would you like to equilibrate the system in order to match the temperature in file input.dat? ( y / n )\n";
		cin >> ans;
	} while (ans != "y" && ans != "n");
	
	if (ans == "y"){
		for (int i=0; i<nequil; i++){
			Rescale();
			for (int j=0; j<nequil_step; j++){
				Move();
				Measure();
				PrintMeasure(1);
			}
			
		}
	}
	
	return;
}

void MolDynamics::Measure(void){ //Properties measurement
	int bin;
	double v, t, vij, p;
	double dx, dy, dz, dr;
	
	v = 0.0; //reset observables
	t = 0.0;
	p = 0.0;
	
	//reset the hystogram of g(r)
	for (int k=igofr; k<igofr+nbins; ++k)
		walker[k]=0.0;
	
	//cycle over pairs of particles
	for (int i=0; i<npart-1; ++i){
		for (int j=i+1; j<npart; ++j){
			
			dx = Pbc( xold[i] - xold[j] ); // here I use old configurations [old = r(t)]
			dy = Pbc( yold[i] - yold[j] ); // to be compatible with EKin which uses v(t)
			dz = Pbc( zold[i] - zold[j] ); // => EPot should be computed with r(t)
			
			dr = dx*dx + dy*dy + dz*dz;
			dr = sqrt(dr);
			
			//update histogram of g(r)
			bin = (int) (dr/bin_size);
			
			if (bin < 100)
				walker[igofr+bin] = walker[igofr+bin] + 2.0;
			
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
	
	walker[iv] = v;
	walker[ik] = 0.5 * t;
	walker[ie] = v + t;
	walker[it] = t / 3.;
	walker[ip] = rho * (t + 48. * p) / 3.;
	
	return;
}

void MolDynamics::Accumulate(void){ //Update block averages
	
	for(int i=0; i<n_props; ++i)
	blk_av[i] = blk_av[i] + walker[i];
	
	blk_norm = blk_norm + 1.0;
}

void MolDynamics::Averages(void){ //Print results for current block
	
	double r, gdir;
	ofstream Gofr, Gave, Epot, Ekin, Temp, Etot, Pres;
	const int wd=16;
	
	cout << "Block number " << iblk << endl;
	
	Epot.open("output/"+state+"/output.epot.0",ios::app);
	Ekin.open("output/"+state+"/output.ekin.0",ios::app);
	Etot.open("output/"+state+"/output.etot.0",ios::app);
	Temp.open("output/"+state+"/output.temp.0",ios::app);
	Pres.open("output/"+state+"/output.pres.0",ios::app);
	Gofr.open("output/"+state+"/output.gofr.0",ios::app);
	Gave.open("output/"+state+"/output.gave.0",ios::app);
	
	stima_pot = blk_av[iv]/blk_norm/(double)npart;  //Potential energy
	glob_av[iv] += stima_pot;
	glob_av2[iv] += stima_pot*stima_pot;
	err_pot = Error(glob_av[iv],glob_av2[iv]);
	
	stima_kin = blk_av[ik]/blk_norm/(double)npart;  //Kinetic energy
	glob_av[ik] += stima_kin;
	glob_av2[ik] += stima_kin*stima_kin;
	err_kin = Error(glob_av[ik],glob_av2[ik]);
	
	stima_etot = blk_av[ie]/blk_norm/(double)npart; //Total energy
	glob_av[ie] += stima_etot;
	glob_av2[ie] += stima_etot*stima_etot;
	err_etot = Error(glob_av[ie],glob_av2[ie]);
	
	stima_temp = blk_av[it]/blk_norm/(double)npart; //Temperature
	glob_av[it] += stima_temp;
	glob_av2[it] += stima_temp*stima_temp;
	err_temp = Error(glob_av[it],glob_av2[it]);
	
	stima_pres = blk_av[ip]/blk_norm/(double)npart; //Pressure
	glob_av[ip] += stima_pres;
	glob_av2[ip] += stima_pres*stima_pres;
	err_pres = Error(glob_av[ip],glob_av2[ip]);
	
	
	
	//Potential energy per particle
	Epot << setw(wd) << iblk <<  setw(wd) << stima_pot << setw(wd) << glob_av[iv]/(double)iblk << setw(wd) << err_pot << endl;
	//Kinetic energy per particle
	Ekin << setw(wd) << iblk <<  setw(wd) << stima_kin << setw(wd) << glob_av[ik]/(double)iblk << setw(wd) << err_kin << endl;
	//Total energy per particle
	Etot << setw(wd) << iblk <<  setw(wd) << stima_etot << setw(wd) << glob_av[ie]/(double)iblk << setw(wd) << err_etot << endl;
	//Temperature
	Temp << setw(wd) << iblk <<  setw(wd) << stima_temp << setw(wd) << glob_av[it]/(double)iblk << setw(wd) << err_temp << endl;
	//Pressure
	Pres << setw(wd) << iblk <<  setw(wd) << stima_pres << setw(wd) << glob_av[ip]/(double)iblk << setw(wd) << err_pres << endl;
	
	
	
	
	//g(r)
	for (int k=igofr; k<igofr+nbins; ++k){
		r = (k-igofr)*bin_size;
		double dv = 4.*pi*(pow(r+bin_size,3)-pow(r,3))/3.;
		
		gdir = blk_av[k]/blk_norm/(double)npart/rho/dv;
		glob_av[k] += gdir;
		glob_av2[k] += gdir*gdir;
		err_gdir = Error(glob_av[k],glob_av2[k]);
		
		Gofr << setw(wd) << iblk << "   " << setw(wd) << k - igofr << "   " << setw(wd) << glob_av[k]/(double)iblk << "   " << setw(wd) << err_gdir << endl;
		
		if (iblk == nblk){
			Gave << setw(wd) << r << "   " << setw(wd) << (double) glob_av[k]/nblk << "   " << setw(wd) << err_gdir << endl;
		}
	}
	
	
	cout << "----------------------------" << endl << endl;
	
	Epot.close();
	Ekin.close();
	Etot.close();
	Temp.close();
	Pres.close();
	Gofr.close();
	Gave.close();
}

void MolDynamics::PrintMeasure(int equil){
	ofstream Epot, Ekin, Etot, Temp, Pres;
	
	if (equil==1){
		Epot.open("output/" + state + "/equilibration.epot.0",ios::app);
		Ekin.open("output/" + state + "/equilibration.ekin.0",ios::app);
		Temp.open("output/" + state + "/equilibration.temp.0",ios::app);
		Etot.open("output/" + state + "/equilibration.etot.0",ios::app);
		Pres.open("output/" + state + "/equilibration.pres.0",ios::app);
	}
	else{
		Epot.open("output/" + state + "/instantaneus.epot.0",ios::app);
		Ekin.open("output/" + state + "/instantaneus.ekin.0",ios::app);
		Temp.open("output/" + state + "/instantaneus.temp.0",ios::app);
		Etot.open("output/" + state + "/instantaneus.etot.0",ios::app);
		Pres.open("output/" + state + "/instantaneus.pres.0",ios::app);
	}

	stima_pot = walker[iv]/(double)npart;  //Potential energy per particle
	stima_kin = walker[ik]/(double)npart;  //Kinetic energy per particle
	stima_temp = walker[it]/(double)npart; //Temperature
	stima_etot = walker[ie]/(double)npart; //Total energy per particle
	stima_pres = walker[ip]/(double)npart; //Pressure
	
	Epot << stima_pot  << endl;
	Ekin << stima_kin  << endl;
	Temp << stima_temp << endl;
	Etot << stima_etot << endl;
	Pres << stima_pres << endl;
	
	
	Epot.close();
	Ekin.close();
	Temp.close();
	Etot.close();
	Pres.close();
	
	return;
}

void MolDynamics::Reset(void){ //Reset block averages
	
	if(iblk == 1){
		for(int i=0; i<n_props; ++i){
			glob_av[i] = 0;
			glob_av2[i] = 0;
		}
	}
	
	for (int i=0; i<n_props; ++i)
		blk_av[i] = 0;
	
	blk_norm = 0;
}



void MolDynamics::ConfOld(void){ //Write old configuration
	ofstream WriteConf;
	
	cout << "Print old configuration to file output/"+state+"/old.final " << endl;
	WriteConf.open("old.final");
	
	for (int i=0; i<npart; ++i){
		WriteConf << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
	}
	WriteConf.close();
	return;
}

void MolDynamics::ConfFinal(void){ //Write final configuration
	ofstream WriteConf;
	
	cout << "Print final configuration to file output/"+state+"/config.final " << endl << endl;
	WriteConf.open("config.final");
	
	for (int i=0; i<npart; ++i){
		WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
	}
	WriteConf.close();
	return;
}

void MolDynamics::ConfXYZ(int nconf){ //Write configuration in .xyz format
	ofstream WriteXYZ;
	
	WriteXYZ.open("output/"+state+"/frames/config_" + to_string(nconf) + ".xyz");
	WriteXYZ << npart << endl;
	WriteXYZ << "This is only a comment!" << endl;
	for (int i=0; i<npart; ++i){
		WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
	}
	WriteXYZ.close();
}

double MolDynamics::Error(double sum, double sum2){
	if (iblk==1)
		return 0.0;
	else
		return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

double MolDynamics::Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
	return r - box * rint(r/box);
}
