#include "../SurvProb.h"
#include <iostream>
#include "TGraph.h"

using namespace std;

int main()
{
	SurvProb surv(3);
	double ms12 = 1e-5;
	surv.SetupProb(10,0.001,0.0241,ms12);
	cout << surv.GetProbFromCalculation(1,1,1) << endl;
	return 0;
	int n=100;
	double x[100],y[100];
	for(int i=0; i<n; i++)
	{
		x[i] = 0.05+i*20./100.;
		surv.SetupProb(x[i], 0.307, 0.0241, ms12);
		y[i] = surv.GetProbFromCalculation(0,1,1);
	}
	TGraph* graph = new TGraph(n,x,y);
	graph->SaveAs("test.root");

	return 0;
}
