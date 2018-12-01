/*
	 JP's result in LMA region.

	 Weiran, May 13, 2018.

*/

#include "../../../JP.h"
#include "../../../../../Solar/SolarNu.h"
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

double Getsin_2_theta(int thetabin);
double Getms(int massbin);
JP jinping(3, model);

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
	TMinuit minuit;
	int ierflg;

	int Bkg[5] = {2,4,5,6,7};
	for(int bkg=0; bkg<5; bkg++)
	{
		minuit.mnparm(bkg, (jinping.GetBkgName(Bkg[bkg])+"day").c_str(), jinping.GetBkgFlux(Bkg[bkg]), sqrt(jinping.GetBkgFlux(Bkg[bkg])),0,0,ierflg);
		minuit.mnparm(bkg+5, (jinping.GetBkgName(Bkg[bkg])+"night").c_str(), jinping.GetBkgFlux(Bkg[bkg]), sqrt(jinping.GetBkgFlux(Bkg[bkg])), 0,0,ierflg);
	}

	int Sig[2] = {5,2};
	for(int sig=0; sig<2; sig++)
		minuit.mnparm(sig+10, solar.GetCompName(Sig[sig]+1), solar.GetFlux(model, Sig[sig]+1), solar.GetError(model, Sig[sig]+1),0,0,ierflg);

	minuit.FixParameter(11); //Fix hep.
	minuit.FixParameter(3); //Fix Be11.
	minuit.FixParameter(1);
	minuit.FixParameter(8);
	minuit.FixParameter(6);

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
	string filepath = "./result/tmp/core" + Core + ".dat";

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
	double BkgFluxday[5];
	double BkgFluxnight[5];
	for(int bkg=0; bkg<5; bkg++)
	{
		BkgFluxday[bkg] = par[bkg];
		BkgFluxnight[bkg] = par[bkg+5];
	}
	BkgFlux[0] = BkgFluxday; BkgFlux[1] = BkgFluxnight;

	double NuFlux[2];
	for(int sig=0; sig<2; sig++) NuFlux[sig] = par[sig+10];
	
	f = jinping.chi2he(BkgFlux, NuFlux);
	f += jinping.B8Constrain(NuFlux[0]);
}

