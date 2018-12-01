/*
	 JP's result in LMA region.

	 Weiran, May 13, 2018.

*/

#include "../Experiment/Jinping/JP.h"
#include "../Experiment/GALLEX/GaConstrain.h"
#include "../Solar/SolarNu.h"
#include "../SurvProb/SurvNSI.h"
#include "../Experiment/KamLAND/chi2Kam.h"
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
chi2Kam kam;

double Be7Ratio = 10.52/89.48;

double Geted(int edbin);
double Geten(int enbin);
JP jinping(2, model, 310);

double ed, en;

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

	ed = Geted(edbin);
	en = Geten(enbin);

	TGraph* prob = nsi.GetProb(0.327, -4.8e-5, 1, 1, ed, en);
	jinping.SetupParameter(prob);
	//Ga.SetupParameter(sin_2_theta, ms);
	TMinuit* minuit = new TMinuit;
	minuit->SetPrintLevel(-1);
	int ierflg;

	for(int bkg=0; bkg<8; bkg++)
	{
		minuit->mnparm(bkg, jinping.GetBkgName(bkg), jinping.GetBkgFlux(bkg)*2, sqrt(jinping.GetBkgFlux(bkg)*2),0,0,ierflg);
	}

	for(int sig=0; sig<9; sig++)
		minuit->mnparm(sig+8, solar.GetCompName(sig+1), solar.GetFlux(model, sig+1), solar.GetError(model, sig+1),0,0,ierflg);

	minuit->FixParameter(16); //Fix F17.
	minuit->FixParameter(11); //Fix Be7_384.
	minuit->FixParameter(10); //Fix hep.
	minuit->FixParameter(6); //Fix Be11.
	minuit->FixParameter(13); //Fix B8
	minuit->FixParameter(4);
	minuit->FixParameter(5);

	minuit->mnparm(17, "penalty", 0, 1, 0, 0, ierflg);
	minuit->FixParameter(17);

	minuit->SetFCN(fcn);
	minuit->SetErrorDef(1);
	double arglist[10];
	arglist[0] = 1e8;
	minuit->mnexcm("MIGRAD", arglist, 0, ierflg);

	double chi2 = 0;
	double edm, errdef;
	int a, b, c;
	minuit->mnstat(chi2, edm, errdef, a, b, c);

	char corenum[5];
	sprintf(corenum, "%d", edbin);
	string Core = corenum;
	string filepath = "./result/tmp2/core" + Core + ".dat";

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
	double BkgFlux[8];
	for(int bkg=0; bkg<8; bkg++)
		BkgFlux[bkg] = par[bkg];

	double NuFlux[9];
	for(int sig=0; sig<9; sig++) NuFlux[sig] = par[sig+8];
	NuFlux[3] = Be7Ratio * NuFlux[4];
	NuFlux[8] = 0; //F17 degenerate with O15.
	
	f = jinping.chi2le(BkgFlux, NuFlux, par[17]);
	//f += jinping.B8Constrain(NuFlux[5]);

	/* Ga constrain. */
	//f += Ga.chi2(NuFlux);

	/* Flux constrain from model. */
	f += jinping.ModelConstrain(NuFlux);
}

