#include <iostream>
#include "TApplication.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "../FastSim.C"

using namespace std;

int main()
{
  cout<<"Test fast detector response simulation"<<endl;

  FastSim fs;

  cout<<"Create testing true spectrum"<<endl;

  // test A
  double energy = 1.5;

	TApplication *app = new TApplication("app",0,0);

  TH1D* DetAS = new TH1D("DetAS","DetAS", 100, 0, 3);
  fs.Convert( DetAS, energy, 3);

  // test B
  TH1D* TrueS = new TH1D("TrueS","TrueS", 100, 0, 3);
  TrueS->SetBinContent(1, 1);
  TrueS->SetBinContent(100, 1);
  for( int b=20; b<=30; b++ )  {
    TrueS->SetBinContent(b, 1);
  }
  for( int b=60; b<=80; b++ )  {
    TrueS->SetBinContent(b, b/70.);
  }
  
  TH1D* DetBS;  
  fs.Convert(DetBS, TrueS, 1);

  // test C
  TH1D* TrueSC = new TH1D("TrueSC","TrueSC", 100, 0, 3.0);
  double center = 2.614;
  double slope = 0.4;
  for( int b=0; b<=100; b++ )  {
    double e = TrueSC->GetBinCenter(b);
    double con = exp(-(center-e)/slope);
    if(e>center) con=0;
    TrueSC->SetBinContent(b, con);
  }

  TH1D* DetCS;
  fs.Convert(DetCS, TrueSC, 1);

  // plot them
  TCanvas c("Detector", "Detector Response", 1000, 750);;
  c.Divide(3,2);

  c.cd(1);
  DetAS->Draw("");
  DetAS->Fit("gaus");

  c.cd(2);
  TrueS->Draw();

  c.cd(5);
  DetBS->Draw();


  TCanvas *c3 = new TCanvas("ExtGam","ExtGam",900,400);

  c3->Divide(2,1);
  TVirtualPad* p1 = c3->cd(1);
  p1->SetLeftMargin(0.15);
  p1->SetBottomMargin(0.15);
  p1->SetRightMargin(0.05);

  TrueSC->Scale(1/TrueSC->Integral());
  TrueSC->Draw();
  TrueSC->SetLineColor(kBlue);
  TrueSC->SetLineWidth(2);
  TrueSC->GetXaxis()->SetTitle("Visible energy [MeV]");
  TrueSC->GetYaxis()->SetTitle("Strength [Arbitary unit]");
  TrueSC->GetYaxis()->SetTitleOffset(1.6);

  TVirtualPad* p2 = c3->cd(2);
  p2->SetLeftMargin(0.15);
  p2->SetBottomMargin(0.15);
  p2->SetRightMargin(0.05);
  
  DetCS->Draw();
  DetCS->SetLineColor(kBlue);
  DetCS->SetLineWidth(2);
  DetCS->GetXaxis()->SetTitle("Visible energy [MeV]");
  DetCS->GetYaxis()->SetTitle("Strength [Arbitary unit]");
  DetCS->GetYaxis()->SetTitleOffset(2.1);

  c3->SaveAs("ExtGamSpec.png");
	//c3->Print("ExtGamSpec.ps");
}
