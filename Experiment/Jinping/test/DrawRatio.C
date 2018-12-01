/*
	 Draw detect over no oscillation spectrum.

	 Weiran, May 28, 2018.

*/

#include "../JP.h"
#include "../../../Target/SolarElec.h"
#include "../../../Detector/FastSim.h"
#include "../../../Solar/SolarNu.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"

using namespace std;

int main()
{
	SolarElec noosci(3, 1);
	SolarElec osci(3, 0);
	cout << "OK" << endl;
	SolarNu solar;
	FastSim sim;

	TFile* file = new TFile("../Simulation/result/solar/GS98-500PE.root");
	TH1D* Day = (TH1D*)file->Get("Jinping2kPredDetday");
	TH1D* Night = (TH1D*)file->Get("Jinping2kPredDetnight");

	double scale = 0.005*24*60*60*3.307e31*20*1500;

	TH1D* noosciT = noosci.GetElecTotal(0, 1);
	TH1D* osciTday = osci.GetElecTotal(0, 1);
	TH1D* osciTnight = osci.GetElecTotal(1, 1);

	TH1D* noosciD;
	TH1D* osciDday;
	TH1D* osciDnight;
	sim.Convert(noosciD, noosciT, 2);
	sim.Convert(osciDday, osciTday, 2);
	sim.Convert(osciDnight, osciTnight, 2);

	int nbins = noosciD->GetNbinsX();
	double xmin = noosciD->GetXaxis()->GetXmin();
	double xmax = noosciD->GetXaxis()->GetXmax();

	TH1D* ratio0 = new TH1D("","",nbins, xmin, xmax);
	ratio0->SetName("ratio");

	for(int i=1; i<=nbins; i++)
	{
		if(noosciD->GetBinContent(i)!=0)
			ratio0->SetBinContent(i, (osciDday->GetBinContent(i)+osciDnight->GetBinContent(i))/(2.*noosciD->GetBinContent(i)));
	}

	TGraph* ratio = new TGraph(ratio0);
	//ratio->SetStats(kFALSE);
	ratio->GetXaxis()->SetRangeUser(0.1, 10);
	ratio->GetYaxis()->SetRangeUser(0.4,0.73);

	ratio->GetXaxis()->SetTitle("Electron Kinetic Energy [MeV]");
	ratio->GetYaxis()->SetTitle("Detectable Ratio");
	ratio->GetXaxis()->SetTitleSize(0.045);
	ratio->GetYaxis()->SetTitleSize(0.045);
	ratio->SetLineColor(kBlue);
	ratio->SetLineWidth(2);
	const int npt = 19;
	double binlimit[npt+1] = {0.101, 0.201, 0.301, 0.401, 0.501, 0.601, 0.701, 0.801, 0.901, 1.001,
		1.101, 1.201, 1.301, 1.501, 2.701, 3.901, 5.101, 6.301, 7.501, 8.701};
	double eCenter[npt], eError[npt];
	double pCenter[npt], pError[npt];
	
	for(int bin=0; bin<npt; bin++)
	{
		eCenter[bin] = (binlimit[bin]+binlimit[bin+1])/2;
		eError[bin] = (binlimit[bin+1]-binlimit[bin])/2;

		int startbin = ratio0->FindBin(binlimit[bin]);
		int endbin = ratio0->FindBin(binlimit[bin+1]);
		int centerbin = ratio0->FindBin(eCenter[bin]);

		double Ntot = Day->Integral(startbin, endbin) + Night->Integral(startbin, endbin);
		double Nnoos = noosciD->Integral(startbin, endbin) * scale;
		pCenter[bin] = ratio0->GetBinContent(centerbin);
		pError[bin] = sqrt(Ntot) / Nnoos; 
	}

	TGraphAsymmErrors *err = new TGraphAsymmErrors(npt, eCenter, pCenter, eError, eError, pError, pError);
	err->SetMarkerColor(kBlack);
	err->SetMarkerStyle(8);
	err->SetMarkerSize(1);
	err->SetLineWidth(2);

	/* For all bins, calculate systematic error from theoretical uncertainty. */
	TH1D* ElecSpec[9];
	TH1D* DetSpec[9];
	for(int comp=1; comp<=9; comp++)
	{
		ElecSpec[comp-1] = osci.GetElecSpec(1, comp);
		sim.Convert(DetSpec[comp-1], ElecSpec[comp-1], 2);
	}

	double sysError[npt];
	for(int bin=0; bin<npt; bin++)
	{
		int lowbin = noosciD->FindBin(binlimit[bin]);
		int upbin = noosciD->FindBin(binlimit[bin+1])-1;

		sysError[bin] = 0;
		double TotalEvents = osciDday->Integral(lowbin, upbin) + osciDnight->Integral(lowbin, upbin);
		for(int comp=1; comp<=9; comp++)
		{
			double events = DetSpec[comp-1]->Integral(lowbin, upbin);
			double percent = 2. * events/TotalEvents;
			double fluxerr = solar.GetModelSeparation(comp)/solar.GetFlux(1, comp);
			if(comp==6) fluxerr = 0.2e-4/5.25e-4;
			sysError[bin] += pow(percent * fluxerr, 2);
		}
		sysError[bin] = sqrt(sysError[bin]);
		cout << "bin: " << bin << "  sysError:  " << sysError[bin] << endl;
	}
			
	double syserr[npt];
	for(int i=0; i<npt; i++) syserr[i] = pCenter[i] * sysError[i];
	TGraphAsymmErrors* sysErr = new TGraphAsymmErrors(npt, eCenter, pCenter, eError, eError, syserr, syserr);
	sysErr->SetFillStyle(1001);
	sysErr->SetFillColorAlpha(16, 0.3);
	sysErr->SetLineWidth(0);
	sysErr->SetLineColor(kWhite);

	TCanvas* c1 = new TCanvas;
	c1->SetLogx();
	TGraph* graph = (TGraph*)ratio->Clone();
	graph->GetXaxis()->SetRangeUser(0.1,10);
	graph->GetYaxis()->SetRangeUser(0.4, 0.8);
	graph->SetLineWidth(0);
	graph->Draw("al");
	ratio->GetXaxis()->SetRangeUser(0.1,10);
	ratio->GetYaxis()->SetRangeUser(0.4, 0.8);
	ratio->Draw("lsame");
	sysErr->Draw("2same");
	err->Draw("Psame");

	TLegend* leg = new TLegend(.5,.6,.85,.85);
	leg->SetBorderSize(0);
	leg->AddEntry(ratio, "Predicted ratio", "l");
	leg->AddEntry(sysErr, "Systematic error from model separation", "f");
	leg->AddEntry(err, "Statistical error", "lep");
	leg->Draw("same");
	

	c1->SaveAs("ratiotest.root");
	return 0;
}



