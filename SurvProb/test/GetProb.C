#include "../SurvProb.h"
#include "TGraph.h"
#include "TFile.h"

using namespace std;

int main()
{
	SurvProb surv(4);
	TGraph* g = surv.GetProb(0.334,4.8e-5,6,0);
	g->SetName("SK");
	TGraph* g2 = surv.GetProb(0.307,7.5e-5, 6, 0);
	g2->SetName("global");

	TGraph* g3 = surv.GetProb(1.2e-2,1e-4, 6, 0, "e");
	g3->SetName("day");

	TGraph* g4 = surv.GetProb(1.2e-2,1e-4, 6, 0, "s");
	g4->SetName("night");

	TGraph* g5 = new TGraph(100);
	for(int i=0; i<100; i++)
	{
		double x = 0.05 + i * 20. / 100.;
		double y = g3->Eval(x)/(1.-g4->Eval(x));
		g5->SetPoint(i,x,y);
	}
	g5->SetName("over");

	TFile* file = new TFile("surv.root", "RECREATE");
	//g->Write();
	//g2->Write();
	g3->Write();
	g4->Write();
	g5->Write();
	file->Close();
	return 0;
}
