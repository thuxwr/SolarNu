/*
	 Try reproducing SK's result.

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

int thetatotal = 300;
int masstotal = 300;

double chi2daynight = 0;

double Getsin_2_theta(int thetabin);
double Getms(int massbin);
chi2SK chi;

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

	/* Get daynight asymmetry from JP's latitude. This is only an approximation. */
	TFile* file = new TFile("../../../Earth/test/result/EarthAffectRegion.root", "READ");
	TH2D* earth = (TH2D*)file->Get("earth");
	double xearth = 5./thetatotal * (thetabin + 0.5) - 4.;
	double yearth = -12. + 9./masstotal * (massbin+0.5);
	int nx = earth->GetXaxis()->FindBin(xearth);
	int ny = earth->GetYaxis()->FindBin(yearth);
	double ADN = earth->GetBinContent(nx, ny);
	chi2daynight = pow((3.3*0.01-ADN)/(1.118*0.01), 2);

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

	char corenum[5];
	sprintf(corenum, "%d", thetabin);
	string Core = corenum;
	string filepath = "../result/tmp/core" + Core + ".dat";

	ofstream fout(filepath.c_str(), ios::app);
	fout << left << setw(5) << thetabin << "\t" << setw(5) << massbin << setw(20) << chi2 << endl;
	fout.close();

	return 0;
}

void fcn(int &npar, double* gin, double &f, double* par, int iflag)
{
	double chi2 = 0;
	chi2 += chi.chi2spec(par[0], par[1]);
	chi2 += pow((par[0]-5.25e6)/0.2e6,2);
	chi2 += pow((par[1]-8e3)/16e3,2);

	/* Add daynight asymmetry. */
	chi2 += chi2daynight;
	f = chi2;
}

double Getsin_2_theta(int thetabin)
{
	/* theta from -4 to 1 */
	double log_tan_2_theta = 5./thetatotal * (thetabin + 0.5) - 4.;
	double tan_2_theta = pow(10, log_tan_2_theta);
	double sin_2_theta = tan_2_theta / (1 + tan_2_theta);
	return sin_2_theta;
}

double Getms(int massbin)
{
	/* mass from -12 to -3. */
	double log_ms = -12. + 9./masstotal * (massbin+0.5);
	double ms = pow(10, log_ms);
	return ms;
}
	

