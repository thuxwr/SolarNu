/*
	 Draw JP's spectrum for CNO energy range.

	 Weiran, May 28, 2018.

*/

#include "../../../../Target/SolarElec.h"
#include "../../../../Detector/FastSim.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"

using namespace std;

int main()
{
	SolarElec elec;
	FastSim sim;
	double scale = 0.005*24*60*60*3.307e31*20*1500;

	TFile* file = new TFile("../../Simulation/result/solar/GS98-500PE.root");

	TH1D* day = (TH1D*)file->Get("Jinping2kPredDetday");
	TH1D* night = (TH1D*)file->Get("Jinping2kPredDetnight");

	TH1D* total = (TH1D*)day->Clone();
	total->Add(day, night);

	TH1D* Det;
	sim.Convert(Det, total, 2);

	/* We draw from 0.7MeV to 1.5MeV. */
	TH1D* Bi210 = (TH1D*)file->Get("Jinping2kBi210");
	TH1D* C11 = (TH1D*)file->Get("Jinping2kC11");
	TH1D* ExtTl208 = (TH1D*)file->Get("Jinping2kExtTl208");
	TH1D* pep = elec.GetElecSpec(1, 2);
	TH1D* N13 = elec.GetElecSpec(1, 7);
	TH1D* O15 = elec.GetElecSpec(1, 8);
	TH1D* B8 = elec.GetElecSpec(1,6);
	TH1D* Be7 = elec.GetElecSpec(1,5);

	TH1D* pepDet0; TH1D* Be7Det0; TH1D* B8Det0; TH1D* N13Det0; TH1D* O15Det0; 
	sim.Convert(pepDet0,pep, 2);
	sim.Convert(Be7Det0,Be7, 2);
	sim.Convert(B8Det0,B8, 2);
	sim.Convert(N13Det0,N13, 2);
	sim.Convert(O15Det0,O15, 2);

	pepDet0->Scale(scale);
	Be7Det0->Scale(scale);
	B8Det0->Scale(scale);
	N13Det0->Scale(scale);
	O15Det0->Scale(scale);
	Bi210->Scale(2);
	C11->Scale(2);
	ExtTl208->Scale(2);

	TGraph* pepDet = new TGraph(pepDet0);
	TGraph* Be7Det = new TGraph(Be7Det0);
	TGraph* B8Det = new TGraph(B8Det0);
	TGraph* N13Det = new TGraph(N13Det0);
	TGraph* O15Det = new TGraph(O15Det0);
	TGraph* Bi210Det = new TGraph(Bi210);
	TGraph* C11Det = new TGraph(C11);
	TGraph* ExtTl208Det = new TGraph(ExtTl208);


	TCanvas* c1 = new TCanvas;
	c1->SetLogy();
	total->GetXaxis()->SetRangeUser(0.6,1.8);
	total->GetYaxis()->SetRangeUser(1,1e5);
	total->SetLineWidth(2);
	total->SetStats(kFALSE);
	total->SetTitle("");
	total->GetXaxis()->SetTitle("Energy [MeV]");
	total->GetXaxis()->SetTitleSize(0.045);
	total->GetYaxis()->SetTitle("Events / 5keV");
	total->GetYaxis()->SetTitleSize(0.045);
	total->Draw();

	pepDet->SetLineWidth(2);
	pepDet->SetLineColor(kBlue);
	Be7Det->SetLineWidth(2);
	Be7Det->SetLineColor(kGreen);
	O15Det->SetLineWidth(2);
	O15Det->SetLineColor(kRed);
	N13Det->SetLineWidth(2);
	N13Det->SetLineColor(6);
	B8Det->SetLineWidth(2);
	B8Det->SetLineColor(42);


	Bi210Det->SetLineWidth(2);
	Bi210Det->SetLineColor(kOrange);
	Bi210Det->SetLineStyle(7);
	C11Det->SetLineWidth(2);
	C11Det->SetLineColor(7);
	C11Det->SetLineStyle(7);
	ExtTl208Det->SetLineWidth(2);
	ExtTl208Det->SetLineColor(38);
	ExtTl208Det->SetLineStyle(7);
	
	TLegend* leg = new TLegend(.4,.6,.85,.85);
	leg->SetNColumns(2);
	leg->SetFillColor(0);
	leg->SetColumnSeparation(.1);
	leg->SetBorderSize(0);
	leg->AddEntry(total,"Total", "l");
	leg->AddEntry(pepDet,"pep","l");
	leg->AddEntry(Be7Det,"^{7}Be","l");
	leg->AddEntry(B8Det,"^{8}B","l");
	leg->AddEntry(N13Det,"^{13}N","l");
	leg->AddEntry(O15Det,"^{15}O","l");
	leg->AddEntry(Bi210Det,"^{210}Bi","l");
	leg->AddEntry(C11Det,"^{11}C","l");
	leg->AddEntry(ExtTl208Det,"External ^{208}Tl","l");
	leg->Draw("same");

	pepDet->Draw("same");
	Be7Det->Draw("same");
	B8Det->Draw("same");
	N13Det->Draw("same");
	O15Det->Draw("same");
	Bi210Det->Draw("same");
	C11Det->Draw("same");
	ExtTl208Det->Draw("same");

	c1->SaveAs("test.root");
	return 0;
}
