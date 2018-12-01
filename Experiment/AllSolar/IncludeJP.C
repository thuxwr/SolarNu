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
#include "../Jinping/JP.h"
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
JP jinping(3, 1, 310);
GaConstrain Ga;
ClConstrain Cl;

SolarNu solar;

double skchi2=0;
double skpar = 0;

void fcn(int &npar, double* gin, double &f, double* par, int iflag);
void B8fcn(int &npar, double* gin, double &f, double* par, int iflag);

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
	jinping.SetupParameter(sin_2_theta, ms);

	TMinuit min;
	int ierflg;
	min.mnparm(0, "B8flux", 5.25e-4, 0.2e-4, 3.25e-4, 7.25e-4, ierflg);
	min.SetErrorDef(1);
	min.SetFCN(B8fcn);
	min.Migrad();

	double chi2 = 0;
	double edm, errdef;
	int nvpar, nparx, icstat;
	min.mnstat(chi2, edm, errdef, nvpar, nparx, icstat);

	double B8flux, B8err, hepflux, heperr;
	min.GetParameter(0, B8flux, B8err);
	//min.GetParameter(2, hepflux, heperr);
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
	
void B8fcn(int &npar, double* gin, double &f, double* par, int iflag)
{
	TMinuit* minuit = new TMinuit;
	minuit->SetPrintLevel(-1);
	int ierflg;

	for(int bkg=0; bkg<8; bkg++)
	{
		minuit->mnparm(bkg, (jinping.GetBkgName(bkg)+"day").c_str(), jinping.GetBkgFlux(bkg), sqrt(jinping.GetBkgFlux(bkg)), 0,0,ierflg);
		minuit->mnparm(bkg+8, (jinping.GetBkgName(bkg)+"night").c_str(), jinping.GetBkgFlux(bkg), sqrt(jinping.GetBkgFlux(bkg)), 0,0,ierflg);
	}

	for(int sig=0; sig<9; sig++)
	{
		if(sig==5) //B8
			minuit->mnparm(sig+16, solar.GetCompName(sig+1), par[0], 0, 0, 0, ierflg);
		else
			minuit->mnparm(sig+16, solar.GetCompName(sig+1), solar.GetFlux(1, sig+1), solar.GetError(1, sig+1), 0, 0, ierflg);
	}

	minuit->FixParameter(19); //Be7_384
	minuit->FixParameter(24); //F17
	minuit->FixParameter(18); //hep
	minuit->FixParameter(6);
	minuit->FixParameter(14);
	minuit->FixParameter(21); //B8, fix to par[0]. 

	minuit->mnparm(25, "penalty", 0, 1, 0, 0, ierflg);

	minuit->SetFCN(fcn);
	minuit->SetErrorDef(1);

	double arglist[10];
	arglist[0] = 1e5;
	minuit->mnexcm("MIGRAD",arglist,0,ierflg);

	double chi2 = 0;
	double edm, errdef;
	int nvpar, nparx, icstat;
	minuit->mnstat(chi2, edm, errdef, nvpar, nparx, icstat);

	f = chi2;
	cout << "B8 flux: " << par[0] << "  chi2: " << chi2 << endl;
	delete minuit;
}


void fcn(int &npar, double* gin, double &f, double* par, int iflag)
{
	double* BkgFlux[2];
	double BkgFluxday[8];
	double BkgFluxnight[8];
	for(int bkg=0; bkg<8; bkg++)
	{
		BkgFluxday[bkg] = par[bkg];
		BkgFluxnight[bkg] = par[bkg+8];
	}
	BkgFlux[0] = BkgFluxday; BkgFlux[1] = BkgFluxnight;

	double NuFlux[9]; //Flux for Ga and Cl experiments.
	for(int sig=0; sig<9; sig++) NuFlux[sig] = par[sig+16];
	NuFlux[3] = Be7Ratio * NuFlux[4];
	NuFlux[8] = 0;

	double chi2 = 0;

	/* SNO */
	chi2 += sno.chi2(par[21]*1e10);

	/* SK */
	//chi2 += SuperK.chi2spec(par[5]*1e10, par[2]*1e10);
	if(skpar!=par[21])
	{
		skpar = par[21];
		skchi2 = SuperK.chi2spec(par[21]*1e10, par[18]*1e10);
	}
	chi2 += skchi2;


	/* Ga */
	chi2 += Ga.chi2(NuFlux);

	/* Cl */
	chi2 += Cl.chi2(NuFlux);

	/* JP */
	chi2 += jinping.chi2(BkgFlux, NuFlux, par[25]);
	chi2 += jinping.ModelConstrain(NuFlux);

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
	


