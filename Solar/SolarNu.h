/*
	Get all solar neutrino spectrum, flux and error of different models.

	The original version supports the following model:
	  1. BP 2000
		2. BP 2004
		3. BS 2005
		4. GS 98
		5. AGS 09

	To shorten time for reading root files, here only GS98 and AGS09 is included. One can simply modify data/DrawSpectrum.C to add models.

	Components:
		1. pp  2. pep  3. hep  4. Be7-384  5. Be7-862  6. B8  7. N13  8. O15  9. F17

	Fluxes are in units 10^{10} /(cm2 s MeV).

	Weiran, Feb 25, 2018
*/

#ifndef SOLAR_NU
#define SOLAR_NU

#include "TFile.h"
#include "TH1D.h"
#include <string>
#include "TGraph.h"

using namespace std;

class SolarNu
{
	public:
		SolarNu();
		~SolarNu();

		/* Units for flux and error: 10^10 /(cm^2 s). Model: 1. GS98  2. AGS09. */
		double GetFlux(int model, int comp);
		double GetError(int model, int comp);
		double GetModelSeparation(int comp);
		double ModelConstrain(int model, double* nuflux);

		string GetModelName(int model);
		string GetCompName(int comp);
		TH1D* GetSpec(int model, int comp);
		TH1D* GetSpec(int comp, double flux, int option); //option is only used to distinguish this from former one. Flux in unit : *10^10 /cm2 /s.

		/* Electron and neutron number density are regarded unchanged due to different models. Radius is normalized. */
		double GetEDensity(double radius); //Unit: cm^{-3}
		double GetNDensity(double radius);
		double GetPDensity(double radius); //proton number density->quark number density, for NSI study.

		double GetResonantDensity(double Energy, double sin_2_theta, double ms);

		TH1D* GetFluxDist(int model, int comp);

		/* Only 2 models for binned distribution. 1:GS98, 2:AGS09. */
		TH1D* GetBinnedFluxDist(int model, int comp); //Binned flux distribution for numeric calculation.

	private:
		TFile* NuSpec;
		TFile* Density;
		TFile* Dist[2];
		TH1D* Distribution[2][9];  
		TFile* BinnedDist[2];
		TH1D* BinnedDistribution[2][9];
		TGraph* EDensity;
		TGraph* NDensity;
		TGraph* PDensity;
		string ModelName[2];
		string Model[2];
		string CompName[9];
		double Flux[5][9];
		double Error[5][9];
		double Be7_384Ratio;
		double NA; //Avogadro constant

};

#endif

