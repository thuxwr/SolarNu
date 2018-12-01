/*
	 Analyze sterile constrain from SK.

	 Weiran, May 15, 2018.

*/

#include "../../Experiment/SuperK/chi2SK.h"
#include "../Sterile.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <iomanip>
#include "TMath.h"
#include "TMinuit.h"

using namespace std;

int thetatotal = 300;
int masstotal = 300;

double Getsin_2_theta(int thetabin);
double Getms(int massbin);

Sterile sterile;
chi2SK SuperK(4);
double ncratio;

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

	SuperK.SetupParameter(sin_2_theta, ms);
	ncratio = sterile.ncratio(sin_2_theta, ms);

	TMinuit minuit;
	int ierflg;
	minuit.mnparm(0, "B8flux", 5.25e6, 0.2e6, 0, 0, ierflg);
	minuit.mnparm(1, "hepflux", 8e3, 16e3, 0, 0, ierflg);
	minuit.SetFCN(fcn);
	minuit.SetErrorDef(1);
	minuit.SetMaxIterations(500);

	minuit.Migrad();

	double chi2 = 0;
	double edm, errdef;
	int nvpar, nparx, icstat;
	minuit.mnstat(chi2, edm, errdef, nvpar, nparx, icstat);

	double B8flux, B8err, hepflux, heperr;
	minuit.GetParameter(0, B8flux, B8err);
	minuit.GetParameter(1, hepflux, heperr);

	SuperK.chi2spec(B8flux, hepflux);
	chi2 += SuperK.chi2tv();

	char corenum[5];
	sprintf(corenum, "%d", thetabin);
	string Core = corenum;
	string filepath = "./tmp/core" + Core + ".dat";

	ofstream fout(filepath.c_str(), ios::app);
	fout << left << setw(5) << thetabin << "\t" << setw(5) << massbin << setw(20) << chi2 << endl;
	fout.close();


	return 0;

}

void fcn(int &npar, double* gin, double &f, double* par, int iflag)
{
	double chi2 = 0;
	chi2 += SuperK.chi2spec(par[0], par[1]);
	chi2 += pow((par[0]*ncratio-5.25e6)/0.2e6, 2);
	chi2 += pow((par[1]-8e3)/16e3, 2);
	f = chi2;
//	cout << "B8flux " << par[0] << "  hepflux " << par[1] << " chi2 " << chi2 << endl;
}

double Getsin_2_theta(int thetabin)
{
	double log_tan_2_theta = 7./thetatotal * (thetabin+0.5) - 6.;
	double tan_2_theta = pow(10, log_tan_2_theta);
	double sin_2_theta = tan_2_theta / (1+tan_2_theta);
	return sin_2_theta;
}

double Getms(int massbin)
{
	double log_ms = -12. + 9./masstotal * (massbin+0.5);
	double ms = pow(10, log_ms);
	return ms;
}
