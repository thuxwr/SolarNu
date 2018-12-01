/*
	 Draw sterile preclude region from Ga and Cl experiments.

	 Weiran, June 24, 2018.

*/

#include "../../Experiment/GALLEX/GaConstrain.h"
#include "../../Experiment/Homestake/ClConstrain.h"
#include "../../Solar/SolarNu.h"
#include "TH2D.h"
#include "TMath.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>
#include "TMinuit.h"

using namespace std;

const int model=1;
double Getsin_2_theta(int thetabin);
double Getms(int massbin);
double Be7Ratio = 10.52/89.48;

void fcn(int &npar, double* gin, double &f, double* par, int iflag);
GaConstrain ga(4);
ClConstrain cl(4);
SolarNu solar;
const int thetatotal = 300;
const int masstotal = 300;

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

	ga.SetupParameter(sin_2_theta, ms);
	cl.SetupParameter(sin_2_theta, ms);

	TMinuit minuit;
	int ierflg;
	for(int comp=1; comp<=9; comp++)
		minuit.mnparm(comp-1, solar.GetCompName(comp), solar.GetFlux(model, comp), solar.GetError(model, comp), 0, 0, ierflg);

	minuit.FixParameter(3);
	minuit.SetFCN(fcn);
	minuit.Migrad();

	double chi2 = 0;
	double edm, errdef;
	int a,b,c;
	minuit.mnstat(chi2, edm, errdef, a, b, c);

	char corenum[5];
	sprintf(corenum, "%d", thetabin);
	string Core = corenum;
	string filepath = "./tmp/core" + Core + ".dat";

	ofstream fout(filepath.c_str(), ios::app);
	fout << left << setw(5) << thetabin << "\t" << setw(5) << massbin << setw(20) << chi2 << endl;
	fout.close();
	//cout << "chi2: " << chi2 << endl;

	return 0;

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

void fcn(int &npar, double* gin, double &f, double* par, int iflag)
{
	double* nuflux = par;
	nuflux[3] = Be7Ratio * nuflux[4];
	double chi2 = 0;
	chi2 += solar.ModelConstrain(model, nuflux);
	chi2 += ga.chi2(nuflux);
	chi2 += cl.chi2(nuflux);
	f = chi2;
}
