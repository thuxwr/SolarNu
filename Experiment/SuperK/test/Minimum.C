/*
	 Minimum for spectrum chi2. Each phase needs a minimum to include a total systematic uncertainty for flux measurement.
	 The best fit point is sin_2_theta = 0.310, ms12 = 4.8e-5.

	 Weiran, Mar. 27, 2018.

*/

#include "../chi2SK.h"
#include "TH2D.h"
#include <iostream>
#include "TMinuit.h"
#include <sstream>
#include <fstream>
#include <string>
#include <iomanip>
#include "TMath.h"
#include "TFile.h"

using namespace std;

chi2SK chi;

double sin_2_theta = 0.310, ms = 4.8e-5;

void fcn(int &npar, double* gin, double &f, double* par, int iflag);

int main()
{
	chi.SetupParameter(sin_2_theta, ms);
	TMinuit minuit(2);
	int ierflg;
	minuit.mnparm(0,"B8flux",5.25e6,0.2e6,0,0,ierflg);
	minuit.mnparm(1,"hepflux",8e3, 16e3, 0,8e4, ierflg);
	minuit.SetFCN(fcn);
	minuit.SetErrorDef(1);

	double arglist[10];
	arglist[0] = 1e5;
	minuit.mnexcm("MIGRAD",arglist,0,ierflg);

	double chi2 = 0;
	double edm, errdef;
	int nvpar, nparx, icstat;
	minuit.mnstat(chi2, edm, errdef, nvpar, nparx, icstat);

	cout << "Chi2:   " << chi2 << endl;
	double B8flux, B8err, hepflux, heperr;

	minuit.GetParameter(0,B8flux, B8err);
	minuit.GetParameter(1,hepflux, heperr);

	double chi2spec = chi2 - pow((B8flux-5.25e6)/0.2e6,2) - pow((hepflux-8e3)/16e3,2);
	cout << "Chi2 without constraint:   " << chi2spec << endl;
	return 0;
}

void fcn(int &npar, double* gin, double &f, double* par, int iflag)
{
	double chi2 = 0;
	chi2 += chi.chi2spec(par[0], par[1]);
	chi2 += pow((par[0]-5.25e6)/0.2e6,2);
	chi2 += pow((par[1]-8e3)/16e3,2);

	f = chi2;
}


