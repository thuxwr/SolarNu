/*
	 Test JP's sensitivity for B8 upturn.
	 Use Asimov data set, use horizontal line to fit detected spectrum.

	 Weiran, June 3, 2018.

*/

#include <iostream>
#include "TMath.h"
#include "TMinuit.h"
#include <string>
#include "../../JP.h"
#include "../../../../Solar/SolarNu.h"
#include "TF1.h"

using namespace std;

void probfcn(int &npar, double* gin, double &f, double* par, int iflag);
TF1* SurvProb;

int model=1;
JP jinping(3, model, 310);
SolarNu solar;

int main()
{
	TMinuit minuit;
	//minuit.SetPrintLevel(-1);
	int ierflg;

	/* SurvProb = c0 + c1 * (E-10) + c2 * (E-10)^2. */
	minuit.mnparm(0, "c0", 0.321, 0.019, 0, 0, ierflg);

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
	//string filepath = "./result/solar/Quadratic.log";
	//ofstream fout(filepath.c_str(), ios::app);
	//cout << left << setw(10) << c0 << "\t" << setw(10) << c0err << "\t" << setw(10) << c1 << "\t"
	//	<< setw(10) << c1err << "\t" << setw(10) << c2 << "\t" << setw(10) << c2err << endl;
	//minuit.mnmatu(1);

	return 0;
}

void probfcn(int &npar, double* gin, double &f, double* par, int iflag)
{
	TF1* func = new TF1("SurvProb", "[0]", 0, 20);
	for(int i=0; i<1; i++) func->SetParameter(i, par[i]);
	SurvProb = func;
	jinping.SetupParameter(SurvProb);

	int Bkg[5] = {2,4,5,6,7};
	double BkgFlux[5];
	for(int bkg=0; bkg<5; bkg++)
		BkgFlux[bkg] = 2*jinping.GetBkgFlux(Bkg[bkg]);

	double NuFlux[2];
	NuFlux[0] = 5.25e-4;
	NuFlux[1] = solar.GetFlux(model, 3);


	double chi2;
	chi2 = jinping.chi2he(BkgFlux, NuFlux, 0);

	f = chi2;
	delete func;
}

