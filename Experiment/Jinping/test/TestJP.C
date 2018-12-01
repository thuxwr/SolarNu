#include "../JP.h"
#include "TF1.h"
#include "TCanvas.h"
#include <iostream>
#include <string>
#include "TH1D.h"

using namespace std;

int main()
{
	JP jinping;
	TF1* f = new TF1("","0.3+0.001*(x-10)+0.001*(x-10)*(x-10)",0,20);
	jinping.SetupParameter(f);
	jinping.SetupParameter(0.307, 7.5e-5);
	TH1D* SigDay = jinping.hsig[0][5];
	TH1D* SigNight = jinping.hsig[1][5];
	TH1D* SigTot = jinping.hsigtot[5];

	TCanvas* c1 = new TCanvas;
	SigTot->Draw();
	SigDay->SetLineColor(kRed);
	SigDay->Draw("same");
	SigNight->SetLineColor(kBlue);
	SigNight->Draw("same");

	c1->SaveAs("test.root");

	return 0;
}
