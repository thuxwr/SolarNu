/*
	 Get total chi2 from SK1~SK4 combined.

	 chi2spec: return chi2 for spectrum fit.
	 chi2tv: return chi2 for time variation(unfinished).

	 The definition of chi2 and calculation can be found in published papers.

	 I assume total capture rate to be proportional to total flux, regardless of small changes in cut efficiency.

	 Weiran, Mar 24, 2018.

	 Change chi2 definition to contain daynight asymmetry without using log likelihood.
	 May 2

*/

#ifndef CHI2SK_H
#define CHI2SK_H

#include "../../Target/SolarElec.h"
#include "SK.h"
#include "TMinuit.h"
#include "TH1D.h"
#include "TH2D.h"
#include "../../Detector/Response.h"
#include <string>
#include "TFile.h"

using namespace std;

class chi2SK
{
	public:
		chi2SK(int generation = 3);
		~chi2SK();
		/* chi2 spectrum is calculated by minimization for parameters deltaB, deltaS and deltaR. */
		/* Setup should be called to create B8 and hep spec. */
		void SetupParameter(double sin_2_theta, double ms);
		double chi2spec(double B8flux, double hepflux);

		/* chi2 for time variation is now estimated to be day night asymmetry. */
		double chi2tv();

		/* Here see note. */
		TH1D* GetPredB8(double sin_2_theta, double ms, int phase, int daynight=2);
		TH1D* GetPredhep(double sin_2_theta, double ms, int phase, int daynight=2);

		/* Alpha is the coefficient of correction for additional systematic uncertainty to the rate. */
		double alpha[4];

		/*--------------------------------------------------------------------*/
		/*---------Don't use the following public member functions.-----------*/
		SK SuperK;
		double B8fluxMC;
		double hepfluxMC;
		double B8fluxmin[4];
		double hepfluxmin[4];
		TH1D* SKspec[4];

		/* B8 spec: B8_{osci} / (B8 + hep)_{MC} */
		/* These two change with theta and mass. First dimension: phase, second dimension: day, night, avg. */
		TH1D* B8spec[4][3];
		TH1D* hepspec[4][3];

		/* B8flux and hepflux. Change when calling chi2spec(). */
		double pB8flux;
		double phepflux;

		/* Get chi2 for given bin and given phase, without penalty terms. */
		double GetBinChi2(int phase, int bin, double* par, double B8flux, double hepflux);
		/* Get tv chi2 for each phase. */
		double GetPhasetvChi2(int phase, double B8flux, double hepflux);

		int mphase;
		double* mpar;

	private:
		int gene;
		static void skfcn(int &npar, double* gin, double &f, double* par, int iflag); //Final minimization function.
		static void phasefcn(int &npar, double* gin, double &f, double* par, int iflag); //First minimize this.
		double binwidth; //All spec share the same binwidth and nbins.
		double msin_2_theta, mms; 
		double mchi2;
		double mdeltaB, mdeltaS[4], mdeltaR[4];
		double sigma0_2[4]; //Uncertainty for total flux measurement.
		double sigmar_2[4];
		int nbins;
		SolarElec* ElecOsci;
		SolarElec* ElecNoOsci;
		TMinuit* minuit; //deltaB, deltaS and deltaR
		Response res;
		TH1D* GetDetSpec(double sin_2_theta, double ms, int comp, double flux, int phase, int IsOsci /* 0 for osci and 1 for non-osci */, int daynight=2);
		TH1D* UnOscispec[4];

		//TFile* file;
		//TH2D* earth;

		string ExpName(int phase);

		static chi2SK* mSK;

};

#endif
