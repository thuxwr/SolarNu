#include "TGraph.h"
#include "TCanvas.h"
#include <iostream>
#include "TMath.h"

using namespace std;

int Draw()
{
	double x[20];
	for(int i=0; i<20; i++)
		x[i] = i+1;
	/* Solar Best. */
	//double y[20] = {0.275, 0.389, 0.477, 0.551, 0.616, 0.674, 0.728, 0.779, 0.826, 0.871, 0.913, 0.954, 0.993, 1.03, 1.066, 1.101, 1.135, 1.168, 1.2, 1.231};
	double y[20] = {0.0746, 0.1452, 0.2121, 0.2735, 0.334, 0.392, 0.446, 0.5, 0.549, 0.599, 0.648, 0.692, 0.738, 0.779, 0.822, 0.865, 0.902, 0.942, 0.978, 1.017}; //JP B8
	//double y[20] = {1.555, 1.9118, 2.1518, 2.3335, 2.4828, 2.6024, 2.71, 2.801, 2.879, 2.957, 3.028, 3.08, 3.139, 3.19, 3.237, 3.276, 3.322, 3.361, 3.4, 3.432};
	/* Global Best. */
	//double z[20] = {0.708, 1.001, 1.226, 1.416, 1.583, 1.734, 1.874, 2.003, 2.124, 2.239, 2.348, 2.453, 2.553, 2.649, 2.742, 2.832, 2.919, 3.004, 3.086, 3.166};
	double z[20] = {0.396, 0.686, 0.971, 1.193, 1.392, 1.574, 1.741, 1.902, 2.047, 2.190, 2.312, 2.435, 2.557, 2.669, 2.776, 2.878, 2.977, 3.079, 3.179, 3.259}; //JP B8
	//double z[20] = {2.3136, 2.975, 3.4531, 3.7921, 4.0674, 4.2721, 4.4625, 4.629, 4.7718, 4.8951, 5.0212, 5.1224, 5.2263, 5.3237, 5.434, 5.5238, 5.6288, 5.7198, 5.8107, 5.9069};

	TGraph* g1 = new TGraph(20, x, y);
	TGraph* g2 = new TGraph(20, x, z);

	TCanvas* c1 = new TCanvas;
	g2->GetYaxis()->SetRangeUser(0,3.8);
	g2->GetXaxis()->SetRangeUser(0,20);
	g2->GetXaxis()->SetTitle("Year");
	g2->GetYaxis()->SetTitle("Sensitivity(sigma)");
	g2->GetXaxis()->SetTitleSize(0.045);
	g2->GetYaxis()->SetTitleSize(0.045);
	g2->SetLineWidth(2);
	g2->SetLineStyle(7);
	g2->SetTitle("");
	g2->SetMarkerStyle(8);
	g2->Draw("alp");
	g1->SetLineWidth(2);
	g1->SetMarkerStyle(8);
	g1->Draw("lpsame");

	TLegend* leg = new TLegend(.11,.65,.4,.85);
	leg->SetBorderSize(0);
	leg->AddEntry(g1, "Solar Best", "l");
	leg->AddEntry(g2, "Global Best", "l");
	leg->Draw("same");

	c1->SaveAs("test.root");
	return 0;
}

