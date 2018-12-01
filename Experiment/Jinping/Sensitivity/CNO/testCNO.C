/*
	 The same as in LOI.

	 Weiran, May 7, 2018.
*/

#include "../../JP.h"
#include "../../../../Solar/SolarNu.h"
#include "../../../GALLEX/GaConstrain.h"
#include "TMinuit.h"
#include <string>
#include "TMath.h"
#include <iostream>
#include <iomanip>

using namespace std;

void fcn(int &npar, double* gin, double &f, double* par, int iflag);
void fcn2(int &npar, double* gin, double &f, double* par, int iflag);

double Be7Ratio = 10.52/89.48;
int model=1;
JP jinping(3, model, 310); //GS
SolarNu solar;
GaConstrain Ga;
double sin_2_theta=0;
double ms=0;

int main()
{
	TMinuit min;
	int ierflg;
	min.mnparm(0, "sin_2_theta", 0.324, 0.01, 0.1, 0.9, ierflg);
	min.mnparm(1, "ms", 4.8e-5, 1e-5, 1e-5, 2e-4, ierflg);

	min.SetFCN(fcn2);
	min.SetErrorDef(1);
	min.Migrad();

	return 0;



}

void fcn2(int &npar, double* gin, double &f, double* par, int iflag)
{
	jinping.SetupParameter(par[0], par[1]);
	Ga.SetupParameter(par[0], par[1]);

	TMinuit* minuit = new TMinuit;
	minuit->SetPrintLevel(-1);

	int ierflg;
	for(int bkg=0; bkg<8; bkg++)
	{
		minuit->mnparm(bkg, (jinping.GetBkgName(bkg)+"day").c_str(), jinping.GetBkgFlux(bkg), sqrt(jinping.GetBkgFlux(bkg)),0,0,ierflg);
		minuit->mnparm(bkg+8, (jinping.GetBkgName(bkg)+"night").c_str(), jinping.GetBkgFlux(bkg), sqrt(jinping.GetBkgFlux(bkg)), 0,0,ierflg);
	}

	for(int sig=0; sig<9; sig++)
		minuit->mnparm(sig+16, solar.GetCompName(sig+1), solar.GetFlux(model, sig+1), solar.GetError(model, sig+1),0,0,ierflg);

	//for(int i=0; i<25; i++) minuit.FixParameter(i);
	minuit->FixParameter(18); //Fix hep.

	minuit->FixParameter(24); //Fix F17.
	minuit->FixParameter(19); //Fix Be7_384.
	minuit->FixParameter(6); //Fix Be11.
	minuit->FixParameter(14);

	minuit->SetFCN(fcn);
	minuit->SetErrorDef(1);
	double arglist[10];
	arglist[0] = 1e8;
	minuit->mnexcm("MIGRAD", arglist, 0, ierflg);

	double amin, b, c;
	int d, e, g;
	minuit->mnstat(amin, b, c, d, e, g);
	f = amin;
	cout << "sin_2_theta  " << par[0] << endl;
	cout << "ms  " << par[1] << endl;
	cout << "f  " << f << endl;
	delete minuit;
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
	f += jinping.B8Constrain(NuFlux[5]);
	f += Ga.chi2(NuFlux);
}

