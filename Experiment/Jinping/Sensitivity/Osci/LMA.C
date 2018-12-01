/*
	 JP's result in LMA region.

	 Weiran, May 13, 2018.

*/

#include "../../JP.h"
#include "../../../GALLEX/GaConstrain.h"
#include "../../../../Solar/SolarNu.h"
#include <iostream>
#include "TMinuit.h"
#include <sstream>
#include <fstream>
#include <string>
#include <iomanip>
#include "TMath.h"

using namespace std;

int thetatotal = 300;
int masstotal = 300;
int model = 1;

SolarNu solar;
GaConstrain Ga;

double Be7Ratio = 10.52/89.48;

double Getsin_2_theta(int thetabin);
double Getms(int massbin);
JP jinping(3, model, 310);

void fcn(int &npar, double* gin, double &f, double* par, int iflag);

int main(int argc, char** argv)
{
	stringstream ss1;
	ss1 << argv[1];
	stringstream ss2;
	ss2 << argv[2];

	int thetabin, massbin;
	ss1 >> thetabin;
	ss2 >> massbin;

	double sin_2_theta = Getsin_2_theta(thetabin);
	double ms = Getms(massbin);

	jinping.SetupParameter(sin_2_theta, ms);
	Ga.SetupParameter(sin_2_theta, ms);
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

	minuit.mnparm(25, "penalty", 0, 1, 0, 0, ierflg);

	minuit.SetFCN(fcn);
	minuit.SetErrorDef(1);
	double arglist[10];
	arglist[0] = 1e8;
	minuit.mnexcm("MIGRAD", arglist, 0, ierflg);

	double chi2 = 0;
	double edm, errdef;
	int a, b, c;
	minuit.mnstat(chi2, edm, errdef, a, b, c);

	char corenum[5];
	sprintf(corenum, "%d", thetabin);
	string Core = corenum;
	string filepath = "./result/tmp2/core" + Core + ".dat";

	ofstream fout(filepath.c_str(), ios::app);
	fout << left << setw(5) << thetabin << "\t" << setw(5) << massbin << setw(20) << chi2 << endl;
	fout.close();

	return 0;
}



double Getsin_2_theta(int thetabin)
{
	/* tan_2_theta from 0.1 to 0.9 */
	double tan_2_theta = 0.8/thetatotal * (thetabin + 0.5) + 0.1;
	double sin_2_theta = tan_2_theta / (1 + tan_2_theta);
	return sin_2_theta;
}

double Getms(int massbin)
{
	/* mass from 1e-5 to 2e-4 */
	double ms = 1e-5 + 14e-5/masstotal * (massbin+0.5);
	return ms;
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
	
	f = jinping.chi2(BkgFlux, NuFlux, par[25]);
	f += jinping.B8Constrain(NuFlux[5]);

	/* Ga constrain. */
	f += Ga.chi2(NuFlux);

	/* Flux constrain from model. */
	//f += jinping.ModelConstrain(NuFlux);
}

