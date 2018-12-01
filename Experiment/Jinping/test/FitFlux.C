/*
	 The same as in LOI.

	 Weiran, May 7, 2018.
*/

#include "../JP.h"
#include "../../../Solar/SolarNu.h"
#include "TMinuit.h"
#include <string>
#include "TMath.h"
#include <iostream>

using namespace std;

void fcn(int &npar, double* gin, double &f, double* par, int iflag);

double Be7Ratio = 10.52/89.48;
int model=2;
JP jinping(3, model); //GS
SolarNu solar;

int main()
{
	jinping.SetupParameter(0.307, 7.5e-5);

	TMinuit minuit;

	int ierflg;
	for(int bkg=0; bkg<8; bkg++)
	{
		minuit.mnparm(bkg, (jinping.GetBkgName(bkg)+"day").c_str(), jinping.GetBkgFlux(bkg), sqrt(jinping.GetBkgFlux(bkg)),0,0,ierflg);
		minuit.mnparm(bkg+8, (jinping.GetBkgName(bkg)+"night").c_str(), jinping.GetBkgFlux(bkg), sqrt(jinping.GetBkgFlux(bkg)), 0,0,ierflg);
	}

	for(int sig=0; sig<9; sig++)
		minuit.mnparm(sig+16, solar.GetCompName(sig+1), solar.GetFlux(model, sig+1), solar.GetError(model, sig+1),0,0,ierflg);

	minuit.FixParameter(24); //Fix F17.
	minuit.FixParameter(19); //Fix Be7_384.
	minuit.FixParameter(18); //Fix hep.
	minuit.FixParameter(6); //Fix Be11.
	minuit.FixParameter(14);

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

	for(int sig=0; sig<9; sig++)
		minuit.GetParameter(sig+16, NuFlux[sig], NuErr[sig]);
	NuFlux[3] = Be7Ratio * NuFlux[4];
	NuErr[3] = Be7Ratio * NuErr[4];

	cout << "-----------------------------------------------------------------" << endl;
	cout << "Fit results:" << endl;
	for(int sig=0; sig<9; sig++) cout << solar.GetCompName(sig+1) << "  " << NuFlux[sig]
		<< "  \\pm  " << NuErr[sig] << "  " << NuErr[sig]/NuFlux[sig] << endl;
	for(int bkg=0; bkg<8; bkg++) cout << jinping.GetBkgName(bkg) << "  " << BkgFlux[0][bkg]+BkgFlux[1][bkg] 
		<< "  \\pm  " << sqrt(pow(BkgErr[0][bkg],2)+pow(BkgErr[1][bkg],2)) << "  " 
			<< sqrt(pow(BkgErr[0][bkg],2)+pow(BkgErr[1][bkg],2)) / (BkgFlux[0][bkg]+BkgFlux[1][bkg])<< endl;

	return 0;
}



void fcn(int &npar, double* gin, double &f, double* par, int iflag)
{
	double* BkgFlux[2];
	double BkgFluxday[8];
	double BkgFluxnight[8];
	for(int bkg=0; bkg<8; bkg++)
	{
		BkgFluxday[bkg] = par[bkg];
		BkgFluxnight[bkg] = par[bkg+8];
	}
	BkgFlux[0] = BkgFluxday; BkgFlux[1] = BkgFluxnight;

	double NuFlux[9];
	for(int sig=0; sig<9; sig++) NuFlux[sig] = par[sig+16];
	NuFlux[3] = Be7Ratio * NuFlux[4];
	NuFlux[8] = 0; //F17 degenerate with O15.
	
	f = jinping.chi2(BkgFlux, NuFlux);
}

