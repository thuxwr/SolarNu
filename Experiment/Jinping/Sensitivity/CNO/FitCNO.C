/*
	 The same as in LOI.

	 Weiran, May 7, 2018.
*/

#include "../../JP.h"
#include "../../../GALLEX/GaConstrain.h"
#include "../../../../Solar/SolarNu.h"
#include "TMinuit.h"
#include <string>
#include "TMath.h"
#include <iostream>
#include <iomanip>
#include "TFile.h"
#include "TH2D.h"

using namespace std;

TFile* file;
TH2D* all;
void fcn(int &npar, double* gin, double &f, double* par, int iflag);

double Be7Ratio = 10.52/89.48;
int model=1;
JP jinping(3, model, 310); //GS
GaConstrain Ga;
SolarNu solar;
double sin_2_theta=0;
double ms=0;

int main()
{
	/* Setup parameter constrain. */
	file = new TFile("../../../AllSolar/AllSolar.root");
	all = (TH2D*)file->Get("AllSolar");

	TMinuit minuit;
	//minuit.SetPrintLevel(-1);

	int ierflg;
	for(int bkg=0; bkg<8; bkg++)
	{
		if(bkg==3) continue;
		minuit.mnparm(bkg, (jinping.GetBkgName(bkg)+"day").c_str(), jinping.GetBkgFlux(bkg), sqrt(jinping.GetBkgFlux(bkg)),0,0,ierflg);
		minuit.mnparm(bkg+8, (jinping.GetBkgName(bkg)+"night").c_str(), jinping.GetBkgFlux(bkg), sqrt(jinping.GetBkgFlux(bkg)), 0,0,ierflg);
	}

	/*Test*/
	minuit.mnparm(3, "C14day", jinping.GetBkgFlux(3)*1e-10,sqrt(jinping.GetBkgFlux(3))*1e-10, 0, 0, ierflg);
	minuit.mnparm(11, "C14night", jinping.GetBkgFlux(3)*1e-10,sqrt(jinping.GetBkgFlux(3))*1e-10, 0, 0, ierflg);

	for(int sig=0; sig<9; sig++)
		minuit.mnparm(sig+16, solar.GetCompName(sig+1), solar.GetFlux(model, sig+1), solar.GetError(model, sig+1),0,10,ierflg);

	//for(int i=0; i<25; i++) minuit.FixParameter(i);
	minuit.FixParameter(18); //Fix hep.
	minuit.mnparm(25, "sin_2_theta", 0.327, 0.01, 0.01, 0.8, ierflg);
	minuit.mnparm(26, "ms", 4.8e-5, 1e-5, 5e-6, 3e-4, ierflg);

	minuit.FixParameter(24); //Fix F17.
	minuit.FixParameter(19); //Fix Be7_384.
	minuit.FixParameter(6); //Fix Be11.
	minuit.FixParameter(14);
	
	minuit.mnparm(27, "penalty", 0, 1, 0, 0, ierflg);

	minuit.SetFCN(fcn);
	minuit.SetErrorDef(1);
	double arglist[10];
	arglist[0] = 1e8;
	minuit.mnexcm("MIGRAD", arglist, 0, ierflg);

	double BkgFlux[2][8], BkgErr[2][8], NuFlux[9], NuErr[9];
	for(int bkg=0; bkg<8; bkg++)
	{
		minuit.GetParameter(bkg, BkgFlux[0][bkg], BkgErr[0][bkg]);
		minuit.GetParameter(bkg+8, BkgFlux[1][bkg], BkgErr[1][bkg]);
	}

	/*Test*/
	BkgFlux[0][3] *= 1e10; BkgErr[0][3] *= 1e10;
	BkgFlux[1][3] *= 1e10; BkgErr[1][3] *= 1e10;

	for(int sig=0; sig<9; sig++)
		minuit.GetParameter(sig+16, NuFlux[sig], NuErr[sig]);
	NuFlux[3] = Be7Ratio * NuFlux[4];
	NuErr[3] = Be7Ratio * NuErr[4];

	cout << "-----------------------------------------------------------------" << endl;
	cout << "Fit results:" << endl;
	for(int sig=0; sig<9; sig++) cout << setw(10) << solar.GetCompName(sig+1) << setw(20) << NuFlux[sig]
		<< "  \\pm" << setw(15) << NuErr[sig] << setw(20) << NuErr[sig]/NuFlux[sig] << endl;
	for(int bkg=0; bkg<8; bkg++) cout << setw(10) << jinping.GetBkgName(bkg) << setw(20) << BkgFlux[0][bkg]+BkgFlux[1][bkg] 
		<< "  \\pm" << setw(15) << sqrt(pow(BkgErr[0][bkg],2)+pow(BkgErr[1][bkg],2)) << setw(20) 
			<< sqrt(pow(BkgErr[0][bkg],2)+pow(BkgErr[1][bkg],2)) / (BkgFlux[0][bkg]+BkgFlux[1][bkg])<< endl;

	return 0;
}



void fcn(int &npar, double* gin, double &f, double* par, int iflag)
{
	if(!(sin_2_theta==par[25] && ms==par[26]))
	{
		jinping.SetupParameter(par[25], par[26]);
		Ga.SetupParameter(par[25], par[26]);
		sin_2_theta = par[25]; ms = par[26];
	}

	double* BkgFlux[2];
	double BkgFluxday[8];
	double BkgFluxnight[8];
	for(int bkg=0; bkg<8; bkg++)
	{
		BkgFluxday[bkg] = par[bkg];
		BkgFluxnight[bkg] = par[bkg+8];
	}

	/*Test*/
	BkgFluxday[3] *= 1e10;
	BkgFluxnight[3] *= 1e10;


	BkgFlux[0] = BkgFluxday; BkgFlux[1] = BkgFluxnight;

	double NuFlux[9];
	for(int sig=0; sig<9; sig++) NuFlux[sig] = par[sig+16];
	NuFlux[3] = Be7Ratio * NuFlux[4];
	NuFlux[8] = 0; //F17 degenerate with O15.
	
	f = jinping.chi2(BkgFlux, NuFlux, par[27]);
	f += jinping.B8Constrain(NuFlux[5]);
	f += Ga.chi2(NuFlux);
	f += all->Interpolate(par[25], par[26]);
}

