/*
	 All solar result except Borexino.
	 Focus on LMA region.
	 Low energy fluxes are constrained by model.

	 Weiran, Apr. 4, 2018.

*/

#include "../SNO/SNO.h"
#include "../SuperK/chi2SK.h"
#include "../GALLEX/GaConstrain.h" 
#include "../Homestake/ClConstrain.h"
#include "../JP.h"
#include "../../Solar/SolarNu.h"
#include <iostream>
#include "TMinuit.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include "TMath.h"
#include <iomanip>

using namespace std;

int thetatotal = 300;
int masstotal = 300;
double Getsin_2_theta(int thetabin);
double Getms(int massbin);

double Be7Ratio = 10.52/89.48;

/* All experiments, 3nu assumption. */
SNO sno;
chi2SK SuperK;
JP Jinping(3, 1, 310);
GaConstrain Ga;
ClConstrain Cl;

SolarNu solar;

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
	SuperK.SetupParameter(sin_2_theta, ms);
	Ga.SetupParameter(sin_2_theta, ms);
	Cl.SetupParameter(sin_2_theta, ms);

	TMinuit minuit(6);
	int ierflg;
	minuit.mnparm(0,"ppflux",solar.GetFlux(1,1),solar.GetError(1,1),0,0,ierflg);
	minuit.mnparm(1,"pepflux",solar.GetFlux(1,2),solar.GetError(1,2),0,0,ierflg);
	minuit.mnparm(2,"Be7_862flux",solar.GetFlux(1,5),solar.GetError(1,5),0,0,ierflg);
	minuit.mnparm(3,"B8flux",solar.GetFlux(1,6),solar.GetError(1,6),0,0,ierflg);
	minuit.mnparm(4,"N13flux",solar.GetFlux(1,7),solar.GetError(1,7),0,0,ierflg);
	minuit.mnparm(5,"O15flux",solar.GetFlux(1,8)+solar.GetFlux(1,9),solar.GetError(1,8),0,0,ierflg);
	minuit.SetFCN(fcn);
	minuit.SetErrorDef(1);

	/*Debug.*/
	minuit.FixParameter(0);
	minuit.FixParameter(1);
	minuit.FixParameter(2);
	minuit.FixParameter(4);
	minuit.FixParameter(5);


	double arglist[10];
	arglist[0] = 1e5;
	minuit.mnexcm("MIGRAD",arglist,0,ierflg);

	double chi2 = 0;
	double edm, errdef;
	int nvpar, nparx, icstat;
	minuit.mnstat(chi2, edm, errdef, nvpar, nparx, icstat);

	double B8flux, B8err;
	minuit.GetParameter(3, B8flux, B8err);
	SuperK.chi2spec(B8flux*1e10, 7.88e3);
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
	double flux[9]; //Flux for Ga and Cl experiments.
	flux[0] = par[0]; flux[1] = par[1]; flux[2] = 7.88e-7; flux[3] = Be7Ratio * par[2];
	flux[4] = par[2]; flux[5] = par[3]; flux[6] = par[4]; flux[7] = par[5]; flux[8] = 0;
	double chi2 = 0;

	/* SNO */
	chi2 += sno.chi2(par[3]*1e10);

	/* SK */
	chi2 += SuperK.chi2spec(par[3]*1e10, 7.88e3);

	/* Ga */
	//chi2 += Ga.chi2(flux);

	/* Cl */
	//chi2 += Cl.chi2(flux);

	/* Constraint from pp chain flux. */

	f = chi2;
}

double Getsin_2_theta(int thetabin)
{
	/* theta from -4 to 1 */
	double tan_2_theta = 0.8/thetatotal * (thetabin + 0.5) + 0.1;
	double sin_2_theta = tan_2_theta / (1 + tan_2_theta);
	return sin_2_theta;
}

double Getms(int massbin)
{
	/* mass from -12 to -3. */
	double ms = 1e-5 + 14e-5/masstotal * (massbin+0.5);
	return ms;
}
	


