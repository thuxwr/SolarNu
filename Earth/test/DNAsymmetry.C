/*
	 See the effect of earth MSW.
	 Only B8 is calculated.
	 Draw 2*(Flux_night - Flux_day)/(Flux_night + Flux_day) as function of theta and mass.
	 Use GS98 for flux info.

	 Weiran, Mar. 15, 2018.
	
*/

#include "../../SurvProb/SurvProb.h"
#include "../../Solar/SolarNu.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include <iostream>
#include "TMath.h"

SurvProb surv(3);
SolarNu solar;

using namespace std;

double nueflux(double sin_2_t12, double ms12, int daynight);

/* Drawing 100 points costs about 1 sec. So only 1 core for this work. */
int main()
{
	int nbinstheta = 300;
	int nbinsmass = 300;

	TH2D* earth = new TH2D("earth", "", nbinstheta, 0.1, 0.9, nbinsmass, 1e-5, 15e-5);	
	for(int binx=1; binx<=nbinstheta; binx++) for(int biny=1; biny<=nbinsmass; biny++)
	{
		double tan_2_t12 = earth->GetXaxis()->GetBinCenter(binx);
		double sin_2_t12 = tan_2_t12 / (1 + tan_2_t12);
		double ms12 = earth->GetYaxis()->GetBinCenter(biny);
		double fluxday = nueflux(sin_2_t12, ms12, 0);
		double fluxnight = nueflux(sin_2_t12, ms12, 1);
		double dif = 2 * (fluxnight-fluxday)/(fluxday+fluxnight);
		if(dif<=1e-4) dif = 1e-4; //Since I want to set logz, dif should not be zero.

		earth->SetBinContent(binx, biny, dif);
	}
	earth->GetXaxis()->SetTitle("tan^{2}#theta_{12}");
	earth->GetYaxis()->SetTitle("#Delta m^{2}_{12} [eV^{2}]");
	earth->SetStats(kFALSE);
	earth->SaveAs("./result/SKEarthAffectRegion.root");
	return 0;
}

double nueflux(double sin_2_t12, double ms12, int daynight)
{
	double flux = 0;
	int model = 1;
	for(int comp=6; comp<=6; comp++)
	{
		TGraph* survprob = surv.GetProb(sin_2_t12, ms12, comp, daynight, "e", "SK");
		TH1D* spec = solar.GetSpec(model, comp);
		int nbins = spec->GetNbinsX();
		double binwidth = spec->GetBinWidth(1);
		for(int bin=1; bin<=nbins; bin++)
		{
			double energy = spec->GetBinCenter(bin);
			double dflux = spec->GetBinContent(bin);
			double prob = survprob->Eval(energy);
			flux += dflux * binwidth * prob;
		}
		delete survprob;
		delete spec;
	}
	return flux;
}



	


