/*
	 Get Gallium capture rate as function of fluxes and mixing parameters.
	 Gallium threshold is 0.2327MeV, according to SAGE, GALLEX and GNO.

	 Weiran, Mar. 17, 2018.

*/

#ifndef GACAPTURE_H
#define GACAPTURE_H

#include <iostream>
#include "../SurvProb/SurvProb.h"
#include "TGraph.h"
#include "../Solar/SolarNu.h"
#include "TH1D.h"
#include <string>
#include "TFile.h"

using namespace std;

class GaCapture
{
	public:
		GaCapture(int generation=3) 
		{
			string nu = getenv("neutrino");
			if(nu=="")
			{
				cout << "Environment variable 'neutrino' undefined." << endl;
				exit(0);
			}

			string filename = nu + "/CrossSection/data/GaAbsorpCS.root";
			TFile* file = new TFile(filename.c_str(),"READ");
			GaCS = (TGraph*)file->Get("CrossSection");
			GaCSerr = (TGraph*)file->Get("Error");

			threshold = 0.2327; 
			gene = generation;
			surv = new SurvProb(gene);
		}

		~GaCapture()
		{
			delete GaCS;
			delete GaCSerr;
			delete surv;
		}

		/* Get capture rate. Unit: SNU. */
		/* Since Gallex did not observe daynight asymmetry, here use average of daynight capture. */
		void SetCapture(double sin_2_theta, double ms)
		{
			TGraph* prob[9];
			TH1D* spec[9];

			for(int comp=1; comp<=9; comp++)
			{
				GaDet[comp-1]=0;
				GaErr[comp-1]=0;

				prob[comp-1] = surv->GetProb(sin_2_theta, ms, comp, "e");
				spec[comp-1] = solar.GetSpec(comp, 1, 0);
				int nbins = spec[comp-1]->GetNbinsX();
				double binwidth = spec[comp-1]->GetBinWidth(1);

				for(int bin=1; bin<=nbins; bin++)
				{
					double energy = spec[comp-1]->GetBinCenter(bin);
					if(energy<=threshold) continue;

					double dflux = spec[comp-1]->GetBinContent(bin);
					double Prob = prob[comp-1]->Eval(energy);
					double AbsorpCS = GaCS->Eval(energy);
					double AbsorpCSerr = GaCSerr->Eval(energy);
					GaDet[comp-1] += dflux * Prob * binwidth * AbsorpCS;
					GaErr[comp-1] += dflux * Prob * binwidth * AbsorpCSerr;
				}
				delete prob[comp-1]; delete spec[comp-1];
			}
		}

		double GetFlux(double* flux)
		{
			double DetTot = 0;
			for(int comp=1; comp<=9; comp++)
			{
				DetTot += GaDet[comp-1] * flux[comp-1];
			}
			return DetTot;
		}

		double GetError(double* flux)
		{
			double ErrTot = 0;
			for(int comp=1; comp<=9; comp++)
			{
				ErrTot += GaErr[comp-1] * flux[comp-1];
			}
			return ErrTot;
		}


	private:
		int gene;
		double GaDet[9], GaErr[9];
		TGraph* GaCS;
		TGraph* GaCSerr;
		SurvProb* surv;
		SolarNu solar;
		double threshold;

};

#endif
