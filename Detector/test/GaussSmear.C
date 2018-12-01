#include "../Response.h"
#include "TH1D.h"
#include "TRandom.h"

using namespace std;

int main()
{
	Response sim;
	TRandom rand;
	TH1D* GaussRandom = new TH1D("GaussRandom","GaussRandom",1000,0,6);
	for(int n=0; n<1e5; n++)
	{
		GaussRandom->Fill(rand.Gaus()+3);
	}
	TH1D* GaussDet;
	sim.Convert(GaussDet, GaussRandom, "SKI");
	GaussDet->SaveAs("GaussDet.root");
	return 0;
}
	
