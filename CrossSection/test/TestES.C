/*
	 Draw Fig.16 in SK IV, in order to test ES cross section.

	 Weiran, Mar. 25, 2018.

*/

//#include "../NuElasticCS.h"
#include "../Backup/NuElasticCS.h"
#include <iostream>
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"

using namespace std;

int main()
{
	NuElasticCS nu;
	TCanvas c1;
	TGraph* graph = new TGraph(1000);
	TGraph* graph2 = new TGraph(1000);
	for(int i=0; i<1000; i++)
	{
		double Te = 0.001 + 10./1000. * i;
		graph->SetPoint(i, Te, nu.dsdTe(10,Te,0,0));
		graph2->SetPoint(i, Te, nu.dsdTe(10,Te,1,0));
	}

	graph->SetTitle("");
	graph->GetXaxis()->SetTitle("Electron Kinetic Energy [MeV]");
	graph->GetYaxis()->SetTitle("Cross Section [cm^{2}/MeV]");
	graph->GetXaxis()->SetTitleSize(0.045);
	graph->GetYaxis()->SetTitleSize(0.045);
	graph->GetYaxis()->SetRangeUser(0,13e-45);
	graph->Draw("al");
	graph2->SetLineColor(kBlue);
	graph2->SetLineStyle(7);
	graph2->Draw("lsame");

	TLegend* leg = new TLegend(0.5,0.7,0.85,0.88);
	leg->SetBorderSize(0);
	leg->SetFillColor(kWhite);
	leg->AddEntry(graph, "#nu_{#font[0]{e}}", "l");
	leg->AddEntry(graph2, "#nu_{#mu,#tau}", "l");
	leg->Draw("same");
	c1.SaveAs("test.pdf");
	//graph->SaveAs("gr.root");
	return 0;
}
