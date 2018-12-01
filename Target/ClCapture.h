/*
	 Get Chlorine capture rate as function of fluxes and mixing parameters.
	 Chlorine threshold is 0.814 MeV, according to Homestake.

	 Weiran, Mar. 18, 2018.

*/

#ifndef CLCAPTURE_H
#define CLCAPTURE_H

#include <iostream>
#include "../SurvProb/SurvProb.h"
#include "TGraph.h"
#include "../Solar/SolarNu.h"
#include "TH1D.h"
#include <string>
#include "TFile.h"

using namespace std;

class ClCapture
{
	public:
		ClCapture(int generation=3)
		{
			string nu = getenv("neutrino");
			if(nu=="")
			{
				cout << "Environment variable 'neutrino' undefined." << endl;
				exit(0);
			}

			string filename = nu + "/CrossSection/data/ClAbsorpCS.root";
			TFile* file = new TFile(filename.c_str(),"READ");
			ClCS = (TGraph*)file->Get("CrossSection");
			ClCSerr = (TGraph*)file->Get("Error");

			threshold = 0.814; 
			gene = generation;
			surv = new SurvProb(gene);
		}

		~ClCapture(){};

		/* Get capture rate. Unit: SNU. */
		void SetCapture(double sin_2_theta, double ms)
		{
			TGraph* prob[9];
			TH1D* spec[9];

			for(int comp=1; comp<=9; comp++)
			{
				ClDet[comp-1]=0;
				ClErr[comp-1]=0;

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
					double AbsorpCS = ClCS->Eval(energy);
					double AbsorpCSerr = ClCSerr->Eval(energy);
					ClDet[comp-1] += dflux * Prob * binwidth * AbsorpCS;
					ClErr[comp-1] += dflux * Prob * binwidth * AbsorpCSerr;
				}
				delete prob[comp-1]; delete spec[comp-1];
			}
		}

		double GetFlux(double* flux)
		{
			double DetTot = 0;
			for(int comp=1; comp<=9; comp++)
			{
				DetTot += ClDet[comp-1] * flux[comp-1];
			}
			return DetTot;
		}

		double GetError(double* flux)
		{
			double ErrTot = 0;
			for(int comp=1; comp<=9; comp++)
			{
				ErrTot += ClErr[comp-1] * flux[comp-1];
			}
			return ErrTot;
		}


	private:
		int gene;
		double ClDet[9], ClErr[9];
		TGraph* ClCS;
		TGraph* ClCSerr;
		SurvProb* surv;
		SolarNu solar;
		double threshold;

};

#endif


