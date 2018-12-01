/*
	 Survival Probability using adiabatic condition, assuming no interference.
	 First check 3 generation.

*/

#include "../Osci.h"
#include "TGraph.h"

using namespace std;

int main()
{
	Oscillation os(3);
	os.SetupParameters();
	int n = 100;
	double x[n],y[n];
	for(int i=0; i<n; i++)
	{
		x[i] = 0.05 + 20./n * i;
		y[i] = os.SurvProb(x[i],6e25);
	}
	TGraph* sp = new TGraph(n,x,y);
	sp->SetTitle("Survival Probability");
	sp->GetXaxis()->SetTitle("Energy[MeV]");
	sp->GetYaxis()->SetTitle("SurvProb");
	sp->SaveAs("SurvivalProbability.root");
	return 0;
}

