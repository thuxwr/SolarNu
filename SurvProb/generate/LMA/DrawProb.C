/*
	 Draw out prob for e neu to be detected as e, mu, tau, s neu.
	 I separate codes for 2nu, 3nu and 4nu so one can modify and run one of these easily.

	 Weiran, Mar.3, 2018
*/

#include "./SurvProb.h"
#include "../../../Solar/SolarNu.h"
#include <sstream>
#include "TFile.h"
#include <string>
#include "TMath.h"
#include <iostream>

using namespace std;

SolarNu solar;

int main(int argc, char** argv)
{
	stringstream ss1;
	ss1 << argv[1];
	stringstream ss2;
	ss2 << argv[2];

	int theta;
	ss1 >> theta;
	int mass;
	ss2 >> mass;

	char thetanum[5], massnum[5];
	sprintf(thetanum, "%d", theta);
	sprintf(massnum, "%d", mass);
	
	string Theta = thetanum;
	string Mass = massnum;

	string flavorname[4];
	flavorname[0] = "e"; flavorname[1] = "mu"; flavorname[2] = "tau"; flavorname[3] = "s";
	string flavorlatex[4];
	flavorlatex[0] = "e"; flavorlatex[1] = "#mu"; flavorlatex[2] = "#tau"; flavorlatex[3] = "s";
	string daynightname[4];
	daynightname[0] = "day"; daynightname[1] = "nightJP"; daynightname[2] = "nightSNO"; daynightname[3] = "nightSK";

	/* Scan for (theta12, ms12). */
//	if(0)
	if(1)
	{
		/* 3 flavors */
		SurvProb surv(3);
		string path = "../../LMA_inv/theta" + Theta + "/mass" + Mass;
		double tan_2_t12 = 0.1 + (theta+0.5) * 0.8/200.;
		double ms12 = 1e-5 + (mass+0.5) * 14e-5/200.;

		double sin_2_t12 = tan_2_t12 / (1+tan_2_t12);
		double sin_2_t13 = 0.0241;

		TFile* file = new TFile((path+"/SurvProb.root").c_str(),"RECREATE");

		int Npoints;

		/* Choose whether to add Nbins. */
		double L = 4 * TMath::Pi() * 10 / ms12;
		if(L*50 > 7.58122e11) Npoints = 2000; //In VAC the oscillation is strenuous.
		else Npoints = 200;

		/* Draw graphs. */
		int nflavor = 2;
		TGraph* graph[9][nflavor][4]; // 4: day, JP night, SNO night, SK night.
		for(int comp=1; comp<=9; comp++)
		{
			for(int finflavor=1; finflavor<=nflavor; finflavor++)
			{
				for(int daynight=0; daynight<4; daynight++)
				{
					graph[comp-1][finflavor-1][daynight] = new TGraph(Npoints);
				}
			}
		}

		double eng = 0.05;
		double binwidth = 20./Npoints;
		int point = 0;
		while(true)
		{
			double minusms12 = 0.-ms12;
			surv.SetupProb(eng, sin_2_t12, sin_2_t13, minusms12);
			for(int comp=1; comp<=9; comp++)
			{
				for(int finflavor=1; finflavor<=nflavor; finflavor++)
				{
					double prob = surv.GetProbFromCalculation(0, finflavor, comp, 0);
					graph[comp-1][finflavor-1][0]->SetPoint(point, eng, prob);
					prob = surv.GetProbFromCalculation(1, finflavor, comp, 0);
					graph[comp-1][finflavor-1][1]->SetPoint(point, eng, prob);
					prob = surv.GetProbFromCalculation(1, finflavor, comp, 1);
					graph[comp-1][finflavor-1][2]->SetPoint(point, eng, prob);
					prob = surv.GetProbFromCalculation(1, finflavor, comp, 2);
					graph[comp-1][finflavor-1][3]->SetPoint(point, eng, prob);
				}
			}
			eng += binwidth;
			point += 1;
			if(point >= Npoints) break;
		}

		/* Update name and axis title. */
		for(int comp=1; comp<=9; comp++)
		{
			for(int finflavor=1; finflavor<=nflavor; finflavor++)
			{
				for(int daynight=0; daynight<4; daynight++)
				{
					string name = solar.GetCompName(comp) + daynightname[daynight] + flavorname[finflavor-1];
					graph[comp-1][finflavor-1][daynight]->SetTitle(("Probability for #nu_{e} to be detected as #nu_{" + flavorlatex[finflavor-1] + "}").c_str());
					graph[comp-1][finflavor-1][daynight]->SetName(name.c_str());
					graph[comp-1][finflavor-1][daynight]->GetXaxis()->SetTitle("Energy/[MeV]");
					graph[comp-1][finflavor-1][daynight]->GetYaxis()->SetTitle("Probability");
					graph[comp-1][finflavor-1][daynight]->Write();
				}
			}
		}
		file->Close();
	}

	return 0;
}







