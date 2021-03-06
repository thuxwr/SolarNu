#include "../SurvNSI.h"
#include <iostream>
#include "TGraph.h"
#include "TFile.h"

using namespace std;

int main()
{
	SurvNSI surv;
	//double prob = surv.SurvProb(0.307, 7.5e-5, 6, 1,-0.22,-0.3, 10);

	//cout << prob << endl;
	double sin_2_t12 = 0.307;
	double ms12 = 7.5e-5;
	double ed = 0;
	double en = 1;
	TGraph* graph1 = new TGraph(100);
	TGraph* graph2 = new TGraph(100);
	TGraph* graph3 = new TGraph(100);
	for(int point=0; point<100; point++)
	{
		double x = 0.05 + point * 20. / 100;
		double y = surv.SurvProb(0.307, 7.5e-5, 6, 1, ed, en, x);
		graph1->SetPoint(point, x, y);
		y = surv.SurvProb(0.32, 7.35e-5, 6, 2, -0.12, -0.16, x);
		graph2->SetPoint(point, x, y);
		y = surv.SurvProb(0.31, 7.5e-5, 6, 0, 0, 0, x);
		graph3->SetPoint(point, x, y);

	}
	TFile* file = new TFile("hehe.root", "RECREATE");
	graph1->SetName("upNSI");
	graph2->SetName("downNSI");
	graph3->SetName("MSW");
	graph1->Write();
	//graph2->Write();
	//graph3->Write();
	return 0;
}
