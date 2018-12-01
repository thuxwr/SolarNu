/*
	 Draw out prob for e neu to be detected as e, mu, tau, s neu.
	 I separate codes for 2nu, 3nu and 4nu so one can modify and run one of these easily.

	 Weiran, Mar.3, 2018
*/

#include "./SurvProb.h"
#include "../../Solar/SolarNu.h"
#include <sstream>
#include "TFile.h"
#include <string>
#include "TMath.h"

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
	string daynightname[2];
	daynightname[0] = "day"; daynightname[1] = "night";

	/* Scan for (theta12, ms12). */
	if(0)
//	if(1)
	{
		/* 2 flavors */
		SurvProb surv(2);
		string path = "../2nu/theta" + Theta + "/mass" + Mass;
		double log_theta12 = -4 + theta * 5./199.;
		double log_mass12 = -12 + mass * 9./199.;

		double tan_2_t12 = pow(10, log_theta12);
		double sin_2_t12 = tan_2_t12 / (1+tan_2_t12);
		double ms12 = pow(10, log_mass12);
		double sin_2_t13 = 0.0241;

		TFile* file = new TFile((path+"/SurvProb.root").c_str(),"RECREATE");

		int Npoints;

		/* Choose whether to add Nbins. */
		double L = 4 * TMath::Pi() * 10 / ms12;
		if(L*50 > 7.58122e11) Npoints = 2000; //In VAC the oscillation is strenuous.
		else Npoints = 200;

		/* Draw graphs. */
		int nflavor = 2;
		TGraph* graph[9][nflavor][2];
		for(int comp=1; comp<=9; comp++)
		{
			for(int finflavor=1; finflavor<=nflavor; finflavor++)
			{
				for(int daynight=0; daynight<2; daynight++)
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
			surv.SetupProb(eng, sin_2_t12, sin_2_t13, ms12);
			for(int comp=1; comp<=9; comp++)
			{
				for(int finflavor=1; finflavor<=nflavor; finflavor++)
				{
					for(int daynight=0; daynight<2; daynight++)
					{
						double prob = surv.GetProbFromCalculation(daynight, finflavor, comp);
						graph[comp-1][finflavor-1][daynight]->SetPoint(point, eng, prob);
					}
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
				for(int daynight=0; daynight<2; daynight++)
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

	/* Scan for (theta12, ms12). */
	if(0)
//	if(1)
	{
		/* 3 flavors */
		SurvProb surv(3);
		string path = "../3nu/12/theta" + Theta + "/mass" + Mass;
		double log_theta12 = -4 + theta * 5./199.;
		double log_mass12 = -12 + mass * 9./199.;

		double tan_2_t12 = pow(10, log_theta12);
		double sin_2_t12 = tan_2_t12 / (1+tan_2_t12);
		double ms12 = pow(10, log_mass12);
		double sin_2_t13 = 0.0241;

		TFile* file = new TFile((path+"/SurvProb.root").c_str(),"RECREATE");

		int Npoints;

		/* Choose whether to add Nbins. */
		double L = 4 * TMath::Pi() * 10 / ms12;
		if(L*50 > 7.58122e11) Npoints = 2000; //In VAC the oscillation is strenuous.
		else Npoints = 200;

		/* Draw graphs. */
		int nflavor = 3;
		TGraph* graph[9][nflavor][2];
		for(int comp=1; comp<=9; comp++)
		{
			for(int finflavor=1; finflavor<=nflavor; finflavor++)
			{
				for(int daynight=0; daynight<2; daynight++)
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
			surv.SetupProb(eng, sin_2_t12, sin_2_t13, ms12);
			for(int comp=1; comp<=9; comp++)
			{
				for(int finflavor=1; finflavor<=nflavor; finflavor++)
				{
					for(int daynight=0; daynight<2; daynight++)
					{
						double prob = surv.GetProbFromCalculation(daynight, finflavor, comp);
						graph[comp-1][finflavor-1][daynight]->SetPoint(point, eng, prob);
					}
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
				for(int daynight=0; daynight<2; daynight++)
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



	/* Scan for (theta13, ms13). */
	/* To be continued. */


	/* Scan for (theta14, ms14). */
//	if(0)
	if(1)
	{
		/* 4 flavors */
		SurvProb surv(4);
		string path = "../4nu/solar/theta" + Theta + "/mass" + Mass;
		double log_theta14 = -6 + theta * 7./199.;
		double log_mass14 = -12 + mass * 9./199.;

		double tan_2_t14 = pow(10, log_theta14);
		double sin_2_t14 = tan_2_t14 / (1+tan_2_t14);
		double ms14 = pow(10, log_mass14);
		//double sin_2_t12 = 0.307;
		double sin_2_t12 = 0.327;
		//double ms12 = 7.54e-5;
		double ms12 = 4.8e-5;
		
		double sin_2_t13 = 0.0241;
		double ms13 = 2.5e-3;

		TFile* file = new TFile((path+"/SurvProb.root").c_str(),"RECREATE");

		int Npoints;

		/* Choose whether to add Nbins. */
		double L = 4 * TMath::Pi() * 10 / ms12;
		if(L*50 > 7.58122e11) Npoints = 2000; //In VAC the oscillation is strenuous.
		else Npoints = 200;

		/* Draw graphs. */
		int nflavor = 4;
		TGraph* graph[9][nflavor][2];
		for(int comp=1; comp<=9; comp++)
		{
			for(int finflavor=1; finflavor<=nflavor; finflavor++)
			{
				for(int daynight=0; daynight<2; daynight++)
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
			surv.SetupProb(eng, sin_2_t12, sin_2_t13, ms12, ms13, sin_2_t14, ms14);
			for(int comp=1; comp<=9; comp++)
			{
				for(int finflavor=1; finflavor<=nflavor; finflavor++)
				{
					for(int daynight=0; daynight<2; daynight++)
					{
						double prob = surv.GetProbFromCalculation(daynight, finflavor, comp);
						graph[comp-1][finflavor-1][daynight]->SetPoint(point, eng, prob);
					}
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
				for(int daynight=0; daynight<2; daynight++)
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







