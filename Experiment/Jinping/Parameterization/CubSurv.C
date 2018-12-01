/*
	 Parameterized survival probability of B8 neutrinos for JP experiment.
	 o(2), o(3) polynomials and exponential fcns are used here.

	 Energy range is selected as [1.75, 15]MeV, in which only B8 and hep neutrinos exist.

	 Weiran, May 8, 2018.

*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include "TMath.h"
#include "../JP.h"
#include "TMinuit.h"
#include <string>
#include "../../../Solar/SolarNu.h"

using namespace std;

void probfcn(int &npar, double* gin, double &f, double* par, int iflag);
void fluxfcn(int &npar, double* gin, double &f, double* par, int iflag);
double Be7Ratio = 10.52/89.48;
TF1* SurvProb;

int model=1;
JP jinping(3, model, 310);
SolarNu solar;

int main()
{
	TMinuit minuit;
	//minuit.SetPrintLevel(-1);
	int ierflg;

	/* SurvProb = c0 + c1 * (E-10) + c2 * (E-10)^2 + c3 * (E-10)^3. */
	minuit.mnparm(0, "c0", 0.321, 0.019, 0, 0, ierflg);
	minuit.mnparm(1, "c1", 0, 0.0078, 0, 0, ierflg);
	minuit.mnparm(2, "c2", 0, 0.0034, 0, 0, ierflg);
	minuit.mnparm(3, "c3", 0, 0.001, 0, 0, ierflg);

	/* Null hypothesis of c0=0 and c1=0. */
	//minuit.FixParameter(1);
	//minuit.FixParameter(2);

	minuit.SetFCN(probfcn);
	minuit.SetErrorDef(1);
	double arglist[10];
	arglist[0] = 1e6;
	minuit.mnexcm("MIGRAD", arglist, 0, ierflg);

	/* Output results. */
	double c0, c0err, c1, c1err, c2, c2err;
	minuit.GetParameter(0, c0, c0err);
	minuit.GetParameter(1, c1, c1err);
	minuit.GetParameter(2, c2, c2err);
	//string filepath = "./result/solar/Quadratic.log";
	//ofstream fout(filepath.c_str(), ios::app);
	//cout << left << setw(10) << c0 << "\t" << setw(10) << c0err << "\t" << setw(10) << c1 << "\t"
	//	<< setw(10) << c1err << "\t" << setw(10) << c2 << "\t" << setw(10) << c2err << endl;
	//minuit.mnmatu(1);

	return 0;
}

void probfcn(int &npar, double* gin, double &f, double* par, int iflag)
{
	TF1* func = new TF1("SurvProb", "[0]+[1]*(x-10)+[2]*(x-10)*(x-10)+[3]*(x-10)*(x-10)*(x-10)", 0, 20);
	for(int i=0; i<4; i++) func->SetParameter(i, par[i]);
	SurvProb = func;
	jinping.SetupParameter(SurvProb);

	TMinuit* Min = new TMinuit;
	Min->SetPrintLevel(-1);
	int ierflg;

	int Bkg[5] = {2,4,5,6,7};
	for(int bkg=0; bkg<5; bkg++)
		Min->mnparm(bkg, jinping.GetBkgName(Bkg[bkg]), 2*jinping.GetBkgFlux(Bkg[bkg]), sqrt(2*jinping.GetBkgFlux(Bkg[bkg])), 0, 0, ierflg);

	int Sig[2] = {5,2};
	for(int sig=0; sig<2; sig++)
		Min->mnparm(sig+5, solar.GetCompName(Sig[sig]+1), solar.GetFlux(model, Sig[sig]+1), solar.GetError(model, Sig[sig]+1), 0, 0, ierflg);

	/* Add penalty term. */
	Min->mnparm(7, "penalty", 0, 1, 0, 0, ierflg);

	Min->FixParameter(6); //Fix hep.
	Min->FixParameter(3); //Fix Be11.
	Min->FixParameter(1); //Fix C10.

	/*Debug*/
//	for(int bkg=0; bkg<5; bkg++)
//		Min->FixParameter(bkg);

	Min->SetFCN(fluxfcn);
	Min->SetErrorDef(1);
	double arglist[10];
	arglist[0] = 1e6;
	Min->mnexcm("MIGRAD", arglist, 0, ierflg);

	double chi2;
	double edm, errdef;
	int nvpar, nparx, icstat;
	Min->mnstat(chi2, edm, errdef, nvpar, nparx, icstat);
	//Min->mnmatu(1);

	f = chi2;

	delete Min; delete func;
}


void fluxfcn(int &npar, double* gin, double &f, double* par, int iflag)
{
	double chi2 = 0;
	double BkgFlux[5];
	for(int bkg=0; bkg<5; bkg++) BkgFlux[bkg] = par[bkg];

	double NuFlux[2];
	for(int sig=0; sig<2; sig++) NuFlux[sig] = par[sig+5];

	chi2 = jinping.chi2he(BkgFlux, NuFlux, par[7]);
	chi2 += pow((5.25e-4-NuFlux[0])/0.2e-4, 2); //B8 flux constrain.

	f = chi2;
}



