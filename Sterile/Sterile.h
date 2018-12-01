/*
	 1. Calculate ratio of NC flux and total flux.

	 Weiran, June 15, 2018.
*/
#ifndef STERILE_H
#define STERILE_H

#include "TMath.h"
#include <iostream>
#include "../SurvProb/SurvProb.h"
#include "../Solar/SolarNu.h"
#include "../Experiment/SNO/SNO.h"
#include "TGraph.h"
#include "TH1D.h"

class Sterile
{
	public: 
		Sterile(){
			surv = new SurvProb(4);
			snos = new SNO(4);
		};
		~Sterile(){};

		/* Only B8 from SNO NC measurement. */
		double ncratio(double sin_2_theta, double ms)
		{
			//TGraph* probe = surv.GetProb(sin_2_theta, ms, 6, "e");
			TGraph* probs = surv->GetProb(sin_2_theta, ms, 6, "s");
			TH1D* specday = snos->MCspec[0];
			TH1D* specnight = snos->MCspec[1];
			TH1D* spec = (TH1D*)specday->Clone();
			spec->Add(specday, specnight);

			int nbins = spec->GetNbinsX();

			double nctotal = 0;
			for(int bin=1; bin<=nbins; bin++)
			{
				double energy = spec->GetBinCenter(bin);
				double dflux = spec->GetBinContent(bin);
				double ncprob = 1. - probs->Eval(energy);
				nctotal += ncprob * dflux;
			}
			double fluxtotal = spec->Integral();

			double NCratio = nctotal / fluxtotal;

			delete probs; delete specday; delete specnight; delete spec;
			return NCratio;
		}


	private:
		SurvProb* surv;
		SolarNu solar;
		SNO* snos;

};

#endif
