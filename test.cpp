#include<iostream>
#include<fstream>
#include"laser.h"
ofstream pofferr;

int main() {

	pofferr.open("error_messages.txt");

	LaserPulse pump_laser, probe_laser;  // 2 pulses

	double pump_angle = 30;  // degrees
	double pump_wavelength = 0.4; //in mkm
	int    pump_polarization = 1; // p-polarizarion
	char   pump_direction ='r';  // right direction 
	string pump_profile = "1e5+1e13*exp(-log(16)*time*time / 0.01)"; // pulse profile in terms of vll

	pump_laser.allocate(pump_wavelength,
			            pump_angle,
			            pump_polarization,
			            pump_direction,
			            pump_profile);

	double probe_angle = 45;  // degrees
	double probe_wavelength = 0.8;
	int    probe_polarization = 0;  // s-polarizarion
	char   probe_direction ='r';
	string probe_profile = "10.0";

	probe_laser.allocate(probe_wavelength,
			             probe_angle,
			             probe_polarization,
			             probe_direction,
			             probe_profile);



	int ncells = 1000;
	vector<double> xl, xr, eps1, eps2, labs;
	xl.resize(ncells);
	xr.resize(ncells);
	eps1.resize(ncells);
	eps2.resize(ncells);
	labs.resize(ncells);
	

    double time = -0.5;
	double dt = 0.001;
	double dt_las = 0.001;



	while(time < 0.5) {

		for(int i = 0; i < labs.size(); ++i)
			labs[i] = 0.0;


		for(int i = 0; i < ncells; ++i) {
			xl[i] = i; // coordinate of left bound of i-th cell in nm
			xr[i] = i+1; // coordinate of right bound of i-th cell in nm
			eps1[i] = 1; // real part of dielectric func
		    eps2[i] = 0; // imag part of dielectric func
			if(i > 100 && i < 950) {
				eps1[i] = -30; // real part of dielectric func
				eps2[i] = 10; // imag part of dielectric func
			}
		}

		pump_laser.LaserAbsorption(ncells, xl, xr, eps1, eps2, time, dt, labs, dt_las);

		for(int i = 0; i < ncells; ++i) {
			xl[i] = i;
			xr[i] = i+1;
			eps1[i] = -50;
			eps2[i] = 40;
		}

		probe_laser.LaserAbsorption(ncells, xl, xr, eps1, eps2, time, dt, labs, dt_las);

		time += dt;
        
		if(fabs(time) < 1e-7) {
			ofstream qlas("qlas.dat");
			for(int i = 0; i < labs.size(); ++i)
				qlas<<i<<" " << labs[i] << endl;
			qlas.close();
		}
	}

	return 0;


}