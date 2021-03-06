#include "../SurvProb.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLegend.h"
#include <iostream>

using namespace std;

int main()
{
	SurvProb surv(4);
	TGraph* graph0 = surv.GetProb(1e-4, 4e-6, 6, 0);
	graph0->SetName("g1");
	graph0->SetLineWidth(2);
	TGraph* graph1 = surv.GetProb(1e-4, 9e-6, 6, 0);
	graph1->SetName("g2");
	graph1->SetLineColor(kBlue);
	graph1->SetLineWidth(2);
	graph1->SetLineStyle(7);
	TGraph* graph2 = surv.GetProb(0.0001, 1.2e-5, 6, 0);
	graph2->SetLineColor(kRed);
	graph2->SetName("g3");
	graph2->SetLineWidth(2);
	graph2->SetLineStyle(8);
	//TGraph* graph1 = surv.GetProb(0.307, 1.5e-5, 6, 1, "e", "JP");
	//graph1->SetName("JPnight");
	//TGraph* graph2 = surv.GetProb(0.307, 1.5e-5, 6, 1, "e", "SNO");
	//graph2->SetName("SNOnight");
	//TGraph* graph3 = surv.GetProb(0.307, 1.5e-5, 6, 1, "e", "SK");
	//graph3->SetName("SKnight");
	TCanvas* c1 = new TCanvas;
	graph0->GetYaxis()->SetRangeUser(0.2,0.6);
	graph0->SetTitle("");
	graph0->GetYaxis()->SetTitle("Survival Probability");
	graph0->GetXaxis()->SetTitleSize(0.045);
	graph0->GetXaxis()->SetRangeUser(0,20);
	graph0->GetYaxis()->SetTitleSize(0.045);
	graph0->GetYaxis()->SetTitleOffset(1.0);
	graph0->Draw("al");
	graph1->Draw("lsame");
	graph2->Draw("lsame");

	TLegend* leg = new TLegend(.4,.65,.85,.85);
	leg->AddEntry(graph0, "#Delta#font[12]{m}^{2}_{1s} = 4#kern[0.3]{#times}#kern[0.05]{10}^{-6}#kern[0.2]{eV}^{2}, sin^{2}(#it{#theta}_{1s}) = 1#kern[0.3]{#times}#kern[0.05]{10}^{-4}", "l");
	leg->AddEntry(graph1, "#Delta#font[12]{m}^{2}_{1s} = 9#kern[0.3]{#times}#kern[0.05]{10}^{-6}#kern[0.2]{eV}^{2}, sin^{2}(#it{#theta}_{1s}) = 1#kern[0.3]{#times}#kern[0.05]{10}^{-4}", "l");
	leg->AddEntry(graph2, "#Delta#font[12]{m}^{2}_{1s} = 1.2#kern[0.3]{#times}#kern[0.05]{10}^{-5}#kern[0.2]{eV}^{2}, sin^{2}(#it{#theta}_{1s}) = 1#kern[0.3]{#times}#kern[0.05]{10}^{-4}", "l");
	leg->SetBorderSize(0);
	leg->SetFillColor(0);
	leg->Draw("same");

	c1->SaveAs("test.pdf");

	//TFile* file = new TFile("test.root","RECREATE");
	//graph0->Write();
	//graph1->Write();
	//graph2->Write();
	//graph3->Write();
	//file->Close();

	return 0;
}
