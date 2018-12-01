/*
	 Earth density model from Nuclear Physics B 521 (1998) 3-36.
	 This can be modified in future.

	 Weiran, Feb 25, 2018
*/

#ifndef EARTH_H
#define EARTH_H

#include "TMath.h"
#include "TFile.h"
#include "TH1D.h"
#include <string>

using namespace std;

class Earth
{
	public:
		Earth(int Land = 1 /* 1 for land and others for ocean. */);
		~Earth(){};

		/* Get electron and neutron density in unit: NA/cm^3. */
		double GetNe(double radius); 
		double GetNn(double radius);

		double GetDensity(double radius); //Unit: g/cm^3

		TH1D* GetLivetime(string experiment); // Get livetime distribution for Jinping Laboratory.

		/* Get a group of density and the according L in a track. */
		void Intersect(double cos_zenith, int &nseg, double* L, double* edens, double* ndens);


	private:
		const static double NA = 6.022141e23;
		double r[7];
		double Ne[7];
		double Nn[7];
		double rho[7];
		TH1D* livetime[3];

};

#endif

