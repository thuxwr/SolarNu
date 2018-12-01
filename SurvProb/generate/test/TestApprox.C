/*
	 Test approximation using jump probability.

	 Weiran, Mar.7, 2018.

*/

#include "../SurvProb.h"
#include "TGraph.h"

using namespace std;

int main()
{
	SurvProb surv(3);
	double x[200], y[200];
	for(int i=0; i<200; i++)
	{
		x[i] = 0.05 + (double)(i*20)/200.;
		y[i] = surv.GetProb(x[i]);
	}
	TGraph* gr = new TGraph(200, x, y);
	gr->SaveAs("approx.root");
	return 0;
}
