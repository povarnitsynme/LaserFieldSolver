#ifndef _LASER_H_
#define _LASER_H_


#include <fstream>
#include <string>
#include <complex>
#include <vector>

#include "Expression.h"
#include "mayday.h"
#include "matrix2x2.h"

#define PI (3.141592653)
#define LIGHTSPEED_CONST (2.99e+05)
#define FIELD_DAMPING	 20 // exp(-FIELD_DAMPING)
#define EPSILON  1e-15

using std::endl;
using std::ofstream;
using std::string;

static int pcount = 0;
typedef complex<double> cmp_dbl;
extern ofstream pofferr;

struct Ecell {
	double xl, xr;
	cmp_dbl epsilon;
};

double modul(const cmp_dbl a)  {

	return sqrt(a.real() * a.real() + a.imag() * a.imag());

}

cmp_dbl cmplroot(const cmp_dbl a) {

	double mda = modul(a);

	cmp_dbl b(sqrt(0.5 * (mda + a.real())), sqrt(0.5 * (mda - a.real())));

	if(a.imag() < 0)
		return cmp_dbl(b.real(), -b.imag());
		

	return b;

}




double modul_2(const cmp_dbl a)  {

	return a.real() * a.real() + a.imag() * a.imag();

}



class LaserPulse {
private:

	char m_pulse_direction;
	int m_pulse;
	double m_wavelength;
	double m_theta;
	string m_profile;
	char   m_polarization;
	ofstream m_diag;

	cmp_dbl m_rinit;
	double m_phiinit;
    int m_haha;
	double m_savethis;
	double m_integralAbs;
	double m_integralPulse;

	double EE(const double z, const double delta, const cmp_dbl kl, const cmp_dbl epsl, const cmp_dbl Bplus, const cmp_dbl Bminus) {
		if(z < 0 || z > 1) {
			pofferr << "wrong z < 0 || z > 1 in function EE" << endl;
			mayDayAbort("error");
		}

		const cmp_dbl im(0, 1);

		const double sin_theta = sin(m_theta * PI / 180.0);
		const double phi1   = (delta * z * kl * im).real(); 
		const double phi2   = (delta * z * kl * im).imag();
		const cmp_dbl eiphi = cmp_dbl(exp(phi1) * cos(phi2), exp(phi1) * sin(phi2));

              cmp_dbl Ez, Ey;

		if (m_polarization == 'p') {
			Ez = sin_theta / epsl * (Bplus * eiphi + Bminus / eiphi );
			Ey = - cmplroot(epsl - sin_theta * sin_theta) / epsl * (Bplus * eiphi - Bminus / eiphi);
		} else {
			Ez = (Bplus * eiphi + Bminus / eiphi);
			Ey = 0;
		}

		return modul_2(Ez) + modul_2(Ey);

	}


	double Runge5(const double delta, const cmp_dbl kl, const cmp_dbl epsl, const cmp_dbl Bplus, const cmp_dbl Bminus) {

		double dz = 1;
		double cur_zpos = 0;
		double result = 0;
		const double errval = 1e-5;
		double k1, k3, k4, k5;
		while(cur_zpos < 1 && dz > 1e-5) {

			k1 = dz / 3.0 * EE(cur_zpos, delta, kl, epsl, Bplus, Bminus);
			k3 = dz / 3.0 * EE(cur_zpos + dz / 3.0, delta, kl, epsl, Bplus, Bminus);
			k4 = dz / 3.0 * EE(cur_zpos + 0.5 * dz, delta, kl, epsl, Bplus, Bminus);
			k5 = dz / 3.0 * EE(cur_zpos + dz, delta, kl, epsl, Bplus, Bminus);
			double ERR = k1 - 4.5 * k3 + 4 * k4 - 0.5 * k5;

			if(ERR < 5.0 / 32 * errval) {
				result += 0.5 * (k1 + 4 * k4 + k5);
				cur_zpos += dz;
				dz *= 1.1;
			} else if (ERR  > 5 * errval) {
				dz *= 0.5;
				continue;
			} else {
				result += 0.5 * (k1 + 4 * k4 + k5);
				cur_zpos += dz;
			}

			if(cur_zpos + dz > 1) dz = 1 - cur_zpos;

		}

              return result;
	}


public:

	LaserPulse() {
		m_phiinit = 98765.4321;
		m_integralAbs = 0;
		m_integralPulse = 0;
	}


	void allocate(double wavelength,
			double theta,
			char   polarization,
			char   pulse_direction,
			string profile) {

		m_phiinit = 98765.4321;

		m_pulse = pcount++;
		m_pulse_direction = pulse_direction;
		m_wavelength = wavelength;
		m_theta = theta;
		m_polarization = polarization;
		m_profile = profile;

        double omega = (1e9 * 2.0 * PI * LIGHTSPEED_CONST / wavelength);

		char buffer[5]; // 3 digit buffer

		sprintf(buffer, "%d", m_pulse);

		
		string name = string("pulse_") + buffer + string(".dat");

		m_diag.open(name.c_str());
            if(!m_diag.good()) {
			pofferr << "can't create laser output file " << name << endl;
			mayDayAbort("error");
		}
		m_diag << "wavelength = " << m_wavelength <<
			 "; theta = " << m_theta << 
			 "; polarization = " << m_polarization << 
			 "; profile = " << m_profile << endl;


		m_diag << "time   I_las   Re{R0}   Im{R0}   |R0|^2   Phase   |T0|^2   ERR0  IntAbs |R/R0|" << endl;

	
	}


/*
	inputs:
    ncells = number of cells for EM field solution [ ]
	xl = vector of coords of left bounds of cells [nm]
	xl = vector of coords right bounds of cells  [nm]
	eps1 = vector of real parts of permittivity [ ]
	eps2 = vector of imag parts of permittivity [ ] 
	time = current time [ps] 
	dt = timestep [ps]
	outputs:
	labs = vector of laser heat in this time step [kW/cm^3]
	dt_las = recommendation on next timestep [ps] is based on laser time profile
*/
	void LaserAbsorption(const int ncells,
		                 const vector<double>& xl,
		                 const vector<double>& xr,
						 const vector<double>& eps1,
                         const vector<double>& eps2,
						 const double time,
						 const double dt,
						 vector<double>& labs,
						 double& dt_las) {

		const double theta = m_theta * PI / 180.0;	// in radians
		const double k0 = 2 * PI / (m_wavelength * 1e3);		// (nm)^-1

		const cmp_dbl sin2theta(sin(theta) * sin(theta), 0);
		const cmp_dbl un(1.0, 0);
		const cmp_dbl im(0, 1.0);


		int err_code;
	   	Expression expr(m_profile);
	   	double I_las = expr.calculate(time + 0.5 * dt, &err_code);
        if(I_las < 1){
			dt_las = 1e30;
			return;
		} else {
			dt_las = 1.1 * dt;
		}


		double I1 = expr.calculate(time + dt, &err_code);
		double I4 = expr.calculate(time + dt + 0.5 * dt_las, &err_code);
		double I5 = expr.calculate(time + dt + dt_las, &err_code);

		double ERR = fabs((I1 + 4.0 * I4 + I5) / I4 / 6.0 - 1);
		
		while (ERR > 0.05 && I1 + I4 + I5 > 1e3) {
			dt_las *= 0.5;
			I4 = expr.calculate(time + dt + 0.5 * dt_las, &err_code);
			I5 = expr.calculate(time + dt + dt_las, &err_code);
			ERR = fabs((I1 + 4.0 * I4 + I5) / I4 / 6.0 - 1.0);
			if(dt_las < 1e-3) break;
			 
		}
		 
        

	// analysis of domain
		int count = ncells;
		for(int i = 0; i < ncells - 1; ++i) {
			if(fabs(xr[i] - xl[i + 1]) > 1e-10)
				++count;
		}


	// initial electric mesh
		vector<Ecell> electric_mesh(count);



	// electromagnetic cells are filled here
		count = 0;
		cmp_dbl teps;


		for(int i = 0; i < ncells; ++i) {
			electric_mesh[count].xl = xl[i];
			electric_mesh[count].xr = xr[i];
			teps = complex<double>(eps1[i], eps2[i]);

			electric_mesh[count].epsilon = teps;
			++count;


	 // add a void cells
			if(i < ncells - 1)
				if(xr[i] != xl[i + 1]) {
					electric_mesh[count].xl = xr[i];
					electric_mesh[count].xr = xl[i + 1];
					electric_mesh[count].epsilon = un;
					++count;

				}
		}


	// add one void cell ahead and one behind
		Ecell first;
		electric_mesh.insert(electric_mesh.begin(), 1, first);
		electric_mesh[0].xr = electric_mesh[1].xl; 
		electric_mesh[0].xl = -1e9; // the point on the right where the reflection and phase are measured (BR)
		electric_mesh[0].epsilon = un;

		electric_mesh.insert(electric_mesh.end(), 1, first);
		electric_mesh[electric_mesh.size() - 1].xr = 1e9;
		electric_mesh[electric_mesh.size() - 1].xl = electric_mesh[electric_mesh.size() - 2].xr; // (TR) position for thin foils
		electric_mesh[electric_mesh.size() - 1].epsilon = un;

		vector<Ecell> electric_mesh_calc(electric_mesh.size());

		for(int i = 0; i < electric_mesh.size(); ++i) {
			if(m_pulse_direction == 'r') {
				electric_mesh_calc[i].xr = electric_mesh[i].xr;
				electric_mesh_calc[i].xl = electric_mesh[i].xl;
				electric_mesh_calc[i].epsilon = electric_mesh[i].epsilon;
			} else if(m_pulse_direction == 'l'){
				electric_mesh_calc[electric_mesh.size() - 1 - i].xl = -electric_mesh[i].xr;
				electric_mesh_calc[electric_mesh.size() - 1 - i].xr = -electric_mesh[i].xl;
				electric_mesh_calc[electric_mesh.size() - 1 - i].epsilon = electric_mesh[i].epsilon;
			} else {
				pofferr << "wrong pulse direction" << endl;
				mayDayAbort("error");
			}
		}



		Matrix P(2,2), C(2,2), E(2,2), F(1, 2);


		cmp_dbl temp[4], epsl, epsr, kl, kr, eiphi;
		double phi1, phi2, delta;


		temp[0] = un;
		temp[1] = 0;
		temp[2] = 0;
		temp[3] = un;
		E.copy(temp);


		cmp_dbl BR, BT, B0 = un;

		int Ecutoff = electric_mesh_calc.size() - 1;

			
		for(int i = 0; i < electric_mesh_calc.size() - 1; ++i) {

			epsl = electric_mesh_calc[i].epsilon;
			epsr = electric_mesh_calc[i + 1].epsilon;

			kl = cmplroot(epsl - sin2theta) * k0;
			kr = cmplroot(epsr - sin2theta) * k0;

			delta = electric_mesh_calc[i].xr - electric_mesh_calc[i].xl;

			phi1 = (delta * kl * im).real(); 
			phi2 = (delta * kl * im).imag();
			eiphi = cmp_dbl(exp(phi1) * cos(phi2), exp(phi1) * sin(phi2));


			temp[0] = eiphi;
			temp[1] = un / eiphi;
			temp[2] = eiphi;
			temp[3] = - temp[1];
			P.copy(temp);

			temp[0] = 0.5 * un;
			temp[1] = 0.5 * kl / kr;
			if(m_polarization == 'p') temp[1] *= epsr / epsl;
			temp[2] = 0.5 * un;
			temp[3] = - temp[1];
			C.copy(temp);

			E = (C * P) * E;

			BR = -E.getComp(0, 1) * B0 / E.getComp(1, 1);
			BT = E.getComp(0, 0) * B0 +  E.getComp(1, 0) * BR;

			if(modul_2(BT) < exp(-double(FIELD_DAMPING))) {
				Ecutoff = i + 1;
				break;
			}
		}


		cmp_dbl* Bplus   =  new cmp_dbl[Ecutoff + 1];
		cmp_dbl* Bminus  = new cmp_dbl[Ecutoff + 1];
		

		temp[0] = B0;
		temp[1] = BR;
		F.copy(temp);

		Bplus[0] = B0;
		Bminus[0] = BR;
		

		for(int i = 0; i < Ecutoff; ++i) {

			kl = cmplroot(electric_mesh_calc[i].epsilon - sin2theta) * k0;
			kr = cmplroot(electric_mesh_calc[i + 1].epsilon - sin2theta) * k0;

			epsl = electric_mesh_calc[i].epsilon;
			epsr = electric_mesh_calc[i + 1].epsilon;

			delta = electric_mesh_calc[i].xr - electric_mesh_calc[i].xl;

			phi1 = (delta * kl * im).real(); 
			phi2 = (delta * kl * im).imag();
			eiphi = cmp_dbl(exp(phi1) * cos(phi2), exp(phi1) * sin(phi2));


			temp[0] = eiphi;
			temp[1] = un / eiphi;
			temp[2] = eiphi;
			temp[3] = - temp[1];
			P.copy(temp);


			temp[0] = 0.5 * un;
			temp[1] = 0.5 * kl / kr;
			if(m_polarization == 'p') temp[1] *= epsr / epsl;
			temp[2] = 0.5 * un;
			temp[3] = - temp[1];
			C.copy(temp);

			F = (C * P) * F;

			Bplus[i + 1] = F.getComp(0, 0);
			Bminus[i + 1] = F.getComp(0, 1);
		
		}


		double Q0 = 0;

		if(m_pulse_direction == 'r')
			count = 0;
		else
			count = int(labs.size() - 1);
		
		for(int i = 1; i < Ecutoff; ++i) {

			delta = electric_mesh_calc[i].xr - electric_mesh_calc[i].xl;
			epsl = electric_mesh_calc[i].epsilon;
			kl = cmplroot(epsl - sin2theta) * k0;

			double Q = k0 * delta * epsl.imag() * Runge5(delta, kl, epsl, Bplus[i], Bminus[i]);


			if(epsl.imag() < 0) {
				pofferr << "negative imaginary part of epsilon, i = " <<  i << endl;
				mayDayAbort("error");

			}


			if(m_pulse_direction == 'r') {
				labs[count] += 1e4 * Q * I_las / delta; // kW/cm^3
				++count;
			} else {
				labs[count] += 1e4 * Q * I_las / delta; // kW/cm^3
				--count;
			}
				
			Q0 += Q / cos(theta);

		}

		double R0 = modul_2(BR / B0);
		double T0 = modul_2(BT);
		double ERR0 = 1.0 - R0 - Q0 - T0;

		m_integralAbs += (1 - R0) * I_las * dt;
		m_integralPulse += I_las * dt;

		if(fabs(BR.real()) < EPSILON)
			BR = BR + EPSILON;

		if(m_phiinit == 98765.4321) {
			m_haha = 0;
			m_rinit = BR;
			m_phiinit = atan(BR.imag()/BR.real());
			m_savethis = atan(BR.imag()/BR.real());
		}

		

		if(m_savethis < -1 && atan(BR.imag()/BR.real()) > 1)
			--m_haha;
		if(m_savethis > 1 && atan(BR.imag()/BR.real()) < -1)
			++m_haha;

		m_savethis = atan(BR.imag()/BR.real());

	
		m_diag << time << "   " << I_las << "   " << BR.real() << "   " << BR.imag() <<  "   " << 
			modul_2(BR) << "   " << (m_pulse_direction == 'r' ? 1 : -1) * 
			(atan(BR.imag()/BR.real()) + PI * m_haha - m_phiinit) << "   " << 
			T0 << "  " << ERR0 <<  "  " << m_integralAbs / m_integralPulse << " "  << modul(BR/(m_rinit)) << endl;

		delete [] Bplus;
		delete [] Bminus;


	}


};

#endif




