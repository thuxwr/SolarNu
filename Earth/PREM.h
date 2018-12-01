/*
	 Preliminary Reference Earth Model for matter density,  Physics of the Earth and Planetary Interiors, 25 (1981) 297-356.
	 The chemical composition is taken the same as SNO, i.e. Z/A = 0.468 inside R=3480km and Z/A = 0.497 outside.

	 No modification is considered for surface rock density at CJPL.

	 Weiran, Apr. 13, 2018.

*/

#ifndef PREM_H
#define PREM_H

#include <string>
#include "TH1D.h"
#include "TFile.h"

using namespace std;

class PREM
{
	public:
		PREM(int Land = 1 /* 1 for land and others for ocean.*/);
		~PREM(){};

		double GetNe(double radius); //NA/cm^3
		double GetNn(double radius);

		double GetDensity(double radius); //Units: g/cm^3

		void Intersect(double cos_zenith, int &nseg, double* L, double* edens, double* ndens);

		TH1D* GetLivetime(string experiment); // Get livetime distribution for Jinping Laboratory.

	private:
		const static double NA = 6.022141e23;
		int layers; //layers for core and mantle.
		double Boundary[10];
		int land;
		double Rearth; //km
		TH1D* livetime[3];

};

#endif
