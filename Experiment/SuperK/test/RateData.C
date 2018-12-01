/*
	 Draw SK's data with statistical error plus energy-uncorrelated systematic uncertainty.
	 Error has been symmetrized.

	 Weiran, Mar. 24, 2018.

*/

#include "../SK.h"
#include "TH1D.h"
#include <string>
#include "TMath.h"

using namespace std;

string PhaseName[4]={"SK1", "SK2", "SK3", "SK4"};

int main()
{
	SK SuperK;
	string path = "../data/Rate/figure";
	for(int phase=1; phase<=4; phase++)
	{
		int nbins = SuperK.GetNbinsX(phase);
		TH1D* th = new TH1D("spectrum","",nbins,SuperK.GetXbins(phase));
		for(int bin=1; bin<=nbins; bin++)
		{
			double kinenergy = th->GetBinCenter(bin);
			double rate = SuperK.GetRate(bin, 2, phase);
			double mcrate = SuperK.GetRateMC(bin, phase);
			double ratio = rate / mcrate;
			double error1 = SuperK.GetStatError(bin, 2, phase)/rate;
			double error2 = SuperK.GetUncorreSysError(kinenergy, phase)*0.01;
			double error = sqrt(error1*error1 + error2*error2);
			th->SetBinContent(bin, ratio);
			th->SetBinError(bin, error*ratio);
		}
		th->SetStats(kFALSE);
		th->SetTitle((PhaseName[phase-1]+" Spectrum").c_str());
		th->GetXaxis()->SetTitle("E_{kin} in MeV");
		th->GetYaxis()->SetTitle("Data/MC(unoscillated)");
		th->SaveAs((path+"/"+PhaseName[phase-1]+".root").c_str());
		delete th;
	}
	return 0;
}
		


