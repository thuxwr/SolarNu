/*
	 Draw Fig.16 in SK IV, in order to test ES cross section.

	 Weiran, Mar. 25, 2018.

*/

#include "../NuElasticCS.h"
//#include "../Backup/NuElasticCS.h"
#include <iostream>
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"

using namespace std;

int main()
{
	NuElasticCS nu;
	NuElasticCS nuer;
	NuElasticCS nuel;
	NuElasticCS nuer2;
	NuElasticCS nuel2;
	nuer.SetupNSI(0,0.2,0,0);
	nuel.SetupNSI(0.2,0,0,0);
	nuer2.SetupNSI(0,-0.2,0,0);
	nuel2.SetupNSI(-0.2,0,0,0);
	TCanvas c1;
	TGraph* graph = new TGraph(1000);
	TGraph* graph2 = new TGraph(1000);
	TGraph* graph3 = new TGraph(1000);
	TGraph* graph4 = new TGraph(1000);
	TGraph* graph5 = new TGraph(1000);
	for(int i=0; i<1000; i++)
	{
		double Te = 0.001 + 0.862/1000. * i;
		graph->SetPoint(i, Te, nu.dsdTe(0.862,Te,0,0));
		graph2->SetPoint(i, Te, nuer.dsdTe(0.862,Te,0,0));
		graph3->SetPoint(i, Te, nuel.dsdTe(0.862,Te,0,0));
		graph4->SetPoint(i, Te, nuer2.dsdTe(0.862,Te,0,0));
		graph5->SetPoint(i, Te, nuel2.dsdTe(0.862,Te,0,0));
	}

	graph->SetTitle("");
	graph->SetLineWidth(2);
	graph->GetXaxis()->SetTitle("Electron Kinetic Energy [MeV]");
	graph->GetYaxis()->SetTitle("Cross Section [cm^{2}/MeV]");
	graph->GetXaxis()->SetTitleSize(0.045);
	graph->GetYaxis()->SetTitleSize(0.045);
	graph->GetYaxis()->SetRangeUser(0,20e-45);
	graph->Draw("al");
	graph2->SetLineColor(kBlue);
	graph2->SetLineWidth(2);
	graph2->SetLineStyle(7);
	graph2->Draw("lsame");
	graph3->SetLineColor(kRed);
	graph3->SetLineWidth(2);
	graph3->SetLineStyle(9);
	graph3->Draw("lsame");
	graph4->SetLineColor(kOrange);
	graph4->SetLineWidth(2);
	graph4->SetLineStyle(7);
	graph4->Draw("lsame");
	graph5->SetLineColor(kViolet);
	graph5->SetLineWidth(2);
	graph5->SetLineStyle(9);
	graph5->Draw("lsame");

	TLegend* leg = new TLegend(0.5,0.7,0.85,0.88);
	leg->SetBorderSize(0);
	leg->SetFillColor(kWhite);
	leg->AddEntry(graph, "Standard", "l");
	leg->AddEntry(graph2, "#varepsilon_{eR}=0.2", "l");
	leg->AddEntry(graph4, "#varepsilon_{eR}=-0.2", "l");
	leg->AddEntry(graph3, "#varepsilon_{eL}=0.2", "l");
	leg->AddEntry(graph5, "#varepsilon_{eL}=-0.2", "l");
	leg->Draw("same");
	//c1.SaveAs("test.pdf");
	c1.SaveAs("test.root");
	return 0;
}
