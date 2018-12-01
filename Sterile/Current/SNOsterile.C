/*
	 Try reproducing SK's result in LMA region.
	 Linear axis.

	 Weiran, Mar. 28, 2018.

*/

#include "../../Experiment/SNO/SNO.h"
#include "../Sterile.h"
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

double Getsin_2_theta(int thetabin);
double Getms(int massbin);
Sterile sterile;
SNO sno(4);
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

	sno.SetupParameter(sin_2_theta, ms);
	ncratio = sterile.ncratio(sin_2_theta, ms);
	TMinuit minuit(1);
	int ierflg;
	minuit.mnparm(0,"B8flux",5.25e6,0.2e6,0,0,ierflg);
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
	string filepath = "./tmp2/core" + Core + ".dat";

	ofstream fout(filepath.c_str(), ios::app);
	fout << left << setw(5) << thetabin << "\t" << setw(5) << massbin << setw(20) << chi2 << endl;
	fout.close();

	return 0;
}

void fcn(int &npar, double* gin, double &f, double* par, int iflag)
{
	double chi2 = 0;
	chi2 += sno.chi2(par[0]*ncratio);

	f = chi2;
}

double Getsin_2_theta(int thetabin)
{
	/* tan_2_theta from 0.1 to 0.9 */
	double log_tan_2_theta = 7./thetatotal * (thetabin+0.5) - 6.;
	double tan_2_theta = pow(10, log_tan_2_theta);
	double sin_2_theta = tan_2_theta / (1+tan_2_theta);
	return sin_2_theta;
}

double Getms(int massbin)
{
	/* mass from 1e-5 to 2e-4 */
	double log_ms = -12. + 9./masstotal * (massbin+0.5);
	double ms = pow(10, log_ms);
	return ms;
}
	

