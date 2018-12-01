/*
	 Test correlated systematic error in SK.

	 Weiran, Mar. 22, 2018.

*/

#include "../SK.h"
#include <iostream>
#include "TH1D.h"
#include "TMath.h"
#include <string>
#include "TCanvas.h"

string PhaseName[4] = {"SK1","SK2","SK3","SK4"};
using namespace std;

int main()
{
	SK SuperK;
	string path = "../data/Systematic/figure";
	for(int phase=1; phase<=4; phase++)
	{
		TH1D* s1p = new TH1D("","",SuperK.GetNbinsX(phase),SuperK.GetXbins(phase));
		TH1D* s1m = new TH1D("","",SuperK.GetNbinsX(phase),SuperK.GetXbins(phase));
		TH1D* s2p = new TH1D("","",SuperK.GetNbinsX(phase),SuperK.GetXbins(phase));
		TH1D* s2m = new TH1D("","",SuperK.GetNbinsX(phase),SuperK.GetXbins(phase));
		TH1D* s3p = new TH1D("","",SuperK.GetNbinsX(phase),SuperK.GetXbins(phase));
		TH1D* s3m = new TH1D("","",SuperK.GetNbinsX(phase),SuperK.GetXbins(phase));

		for(int bin=1; bin<=SuperK.GetNbinsX(phase); bin++)
		{
			double eng = s1p->GetBinCenter(bin);
			s1p->SetBinContent(bin, SuperK.GetCorreSysError(eng, 1, phase, 0));
			s1m->SetBinContent(bin, SuperK.GetCorreSysError(eng, 1, phase, 1));
			s2p->SetBinContent(bin, SuperK.GetCorreSysError(eng, 2, phase, 0));
			s2m->SetBinContent(bin, SuperK.GetCorreSysError(eng, 2, phase, 1));
			s3p->SetBinContent(bin, SuperK.GetCorreSysError(eng, 3, phase, 0));
			s3m->SetBinContent(bin, SuperK.GetCorreSysError(eng, 3, phase, 1));
		}

		s1p->GetYaxis()->SetRangeUser(-0.15, 0.15);
		TCanvas* c1 = new TCanvas();
		s1p->SetStats(kFALSE);
		s1p->SetTitle(("Correlated Systematic Uncertainty for " + PhaseName[phase-1]).c_str());
		s1p->GetXaxis()->SetTitle("E_{kin} [MeV]");
		s1p->GetYaxis()->SetTitle("Uncertainty");
		s1p->SetLineColor(kRed);
		s1m->SetLineColor(kRed);
		s1p->Draw();
		s1m->Draw("same");
		s2p->SetLineColor(kBlue);
		s2p->Draw("same");
		s2m->SetLineColor(kBlue);
		s2m->Draw("same");
		s3p->SetLineColor(kBlack);
		s3p->Draw("same");
		s3m->SetLineColor(kBlack);
		s3m->Draw("same");
		c1->SaveAs((path+"/"+PhaseName[phase-1]+".root").c_str());
	}
	return 0;
}

