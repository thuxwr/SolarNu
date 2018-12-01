/*
	 JP's result in LMA region.

	 Weiran, May 13, 2018.

*/

#include "../Experiment/Jinping/JP.h"
#include "../Experiment/GALLEX/GaConstrain.h"
#include "../Solar/SolarNu.h"
#include "../SurvProb/SurvNSI.h"
#include <iostream>
#include "TMinuit.h"
#include <sstream>
#include <fstream>
#include <string>
#include <iomanip>
#include "TMath.h"
#include "TGraph.h"

using namespace std;

int edtotal = 300;
int entotal = 300;
int model = 1;

SolarNu solar;
//GaConstrain Ga;
SurvNSI nsi;

double Be7Ratio = 10.52/89.48;

double Geted(int edbin);
double Geten(int enbin);
JP jinping(2, model, 310);

void fcn(int &npar, double* gin, double &f, double* par, int iflag);

int main(int argc, char** argv)
{
	stringstream ss1;
	ss1 << argv[1];
	stringstream ss2;
	ss2 << argv[2];

	int edbin, enbin;
	ss1 >> edbin;
	ss2 >> enbin;

	double ed = Geted(edbin);
	double en = Geten(enbin);

	TGraph* prob = nsi.GetProb(0.307, 7.5e-5, 1, 1, ed, en);
	jinping.SetupParameter(prob);
	//Ga.SetupParameter(sin_2_theta, ms);
	TMinuit minuit;
	int ierflg;

	int Bkg[5] = {2,4,5,6,7};
	for(int bkg=0; bkg<5; bkg++)
	{
		minuit.mnparm(bkg, jinping.GetBkgName(Bkg[bkg]), jinping.GetBkgFlux(Bkg[bkg])*2, sqrt(jinping.GetBkgFlux(Bkg[bkg])*2),0,0,ierflg);
	}
	minuit.FixParameter(3);
	//for(int i=0; i<5; i++) minuit.FixParameter(i);

	for(int sig=0; sig<1; sig++)
		minuit.mnparm(sig+5, solar.GetCompName(6), solar.GetFlux(model, 6), solar.GetError(model, 6),0,0,ierflg);

	minuit.mnparm(6, "hep", solar.GetFlux(model, 3), 0.1,0,0,ierflg);
	minuit.FixParameter(6);

	minuit.mnparm(7, "penalty", 0, 1, 0, 0, ierflg);
	minuit.FixParameter(7);

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
	sprintf(corenum, "%d", edbin);
	string Core = corenum;
	string filepath = "./result/tmp/core" + Core + ".dat";

	ofstream fout(filepath.c_str(), ios::app);
	fout << left << setw(5) << edbin << "\t" << setw(5) << enbin << setw(20) << chi2 << endl;
	fout.close();

	return 0;
}



double Geted(int edbin)
{
	return -1.5 + 3.0*edbin/edtotal;
}

double Geten(int enbin)
{
	return -1.5 + 3.0*enbin/entotal;
}

void fcn(int &npar, double* gin, double &f, double* par, int iflag)
{
	double BkgFlux[5];
	for(int bkg=0; bkg<5; bkg++)
		BkgFlux[bkg] = par[bkg];

	double NuFlux[2];
	NuFlux[0] = par[5];
	NuFlux[1] = par[6];
	
	f = jinping.chi2he(BkgFlux, NuFlux, par[7]);
	f += jinping.B8Constrain(NuFlux[0]);

	/* Ga constrain. */
	//f += Ga.chi2(NuFlux);

	/* Flux constrain from model. */
	//f += jinping.ModelConstrain(NuFlux);
}

