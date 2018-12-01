#include <iostream>
#include "TH1D.h"
#include "TGraph.h"
#include "../FastSim.h"
#include "TCanvas.h"

using namespace std;

int main()
{
  cout<<"Test detector resolution function"<<endl;

  FastSim fs;

  double TeMax = 20;
  TH1D * h = new TH1D( "Resolution", "Resolution", 500, 0, TeMax);

  for( int b=1; b<=500; b++ ) {
    double Te = (b-0.5)*TeMax/500;
    double res = fs.Resolution( Te, 2 );
    h->SetBinContent( b, res );
  }

	TGraph* g = new TGraph(h);
	g->SetTitle("");
	g->GetXaxis()->SetTitle("Energy [MeV]");
	g->GetYaxis()->SetTitle("Energy Resolution [MeV]");
	g->GetXaxis()->SetTitleSize(0.045);
	g->GetYaxis()->SetTitleSize(0.045);
	g->SetLineWidth(2);

	g->GetXaxis()->SetRangeUser(0,20);

  TCanvas * c = new TCanvas;
  g->Draw("al");
	c->SaveAs("test.pdf");
}
