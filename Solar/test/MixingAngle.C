/*
	 Test theta12 in matter.
	 By changing osci.SetupParameters, one can achieve theta12 for any set of parameter.

	 Weiran, Feb 26, 2018
*/

#include "../SolarNu.h"
#include "../../Oscillation/Osci.h"
#include "TGraph.h"
#include "TMath.h"
#include "TCanvas.h"
#include <iostream>

using namespace std;

int main()
{
	double t12 = 0.4;
	double m12 = 7.5e-5;
	SolarNu solar;
	Oscillation osci(3);
	osci.SetupParameters(t12,0.0241,0.5,1e-3,m12,2.5e-3,1e-4);
	double Energy = 10;
	int pointnum = 1000;
	double rad[pointnum], theta[pointnum];
	for(int i=0; i<pointnum; i++)
	{	
		rad[i] = ((double)i)/pointnum;
		double ne = solar.GetEDensity(rad[i]);
		double nn = solar.GetNDensity(rad[i]);
		MatrixXcd Um = osci.GetMatterMixingMatrix(Energy, ne, nn);

		/*
		double s13 = abs(Um(0,2).real());
		double c13 = sqrt(1 - s13*s13);
		double s12 = abs(Um(0,1).real()/c13);
		*/
		double t12 = abs(abs(Um(0,1))/abs(Um(0,0)));
		theta[i] = atan(t12);
		//double t14 = abs(Um(0,3).real()/Um(0,0).real());
		//theta[i] = atan(t14);
	}

	TCanvas* c1 = new TCanvas;
	TGraph* graph = new TGraph(pointnum, rad, theta);
	graph->Draw("apl");

	c1->SaveAs("theta.root");

	return 0;
}

