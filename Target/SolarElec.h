/*
	 Detect solar neutrino via elastic scattering.
	 Add CC and NC together.

	 Weiran, Mar.7, 2018.

*/

#ifndef SOLARELEC_H
#define SOLARELEC_H

#include "../Solar/SolarNu.h"
#include "../SurvProb/SurvProb.h"
#include "../Detector/FastSim.h"
#include "../CrossSection/NuElasticCS.h"
#include "TGraph.h"
#include "TF1.h"
#include "TH1D.h"
#include <string>

class SolarElec
{
	public:
		SolarElec(int generation=3, int Option=0 /* 0:oscillation, 1:non-oscillation */);
		~SolarElec(){};

		/* Electron spectrum from solar neutrino. */
		TH1D* GetElecSpec(int daynight, int model, int comp); // For model input. Not supported for sterile neutrinos.
		TH1D* GetElecSpec(int model, int comp); //Day-night averaged.
		TH1D* GetElecSpec(double sin_2_theta, double ms, int daynight, int comp, double flux, string experiment = "JP");
		TH1D* GetElecSpec(double sin_2_theta, double ms, int comp, double flux, string experiment = "JP");
		TH1D* GetElecSpec(TF1* SurvProb, int comp, double flux);
		TH1D* GetElecSpec(TGraph* SurvProb, int comp, double flux);

		TH1D* GetElecTotal(int daynight, int model, double sin_2_theta = 0.327, double ms = 4.8e-5);
		TH1D* GetElecTotal(double sin_2_theta, double ms, int daynight, double* flux);

		double GetTotalCount(int model, int comp);
		double GetDetCount(int model, int comp, double thres);

		/* For NSI. One should first set up NSI and then use GetElecSpec. */
		void SetupNSI(double el, double er, double tl, double tr);

	private:
		bool IsNSI;
		int option;
		int gene;
		SolarNu solar;
		SurvProb* surv;
		NuElasticCS scat;
		TH1D* spec[9];
		double binwidth; //Output binwidth 5keV.
		




};

#endif
