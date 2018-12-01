/*
	 Get chi2 from JP simulation.

	 Fitting method: fit day spectra and night spectra respectively.

	 Fiducial mass: 2000 ton, run time: 1500 days, resolution: 500PE/MeV.
	 Livetime: assuming full livetime. Zenith distribution is approximated by averaging in one year.

	 Return chi square from JP only. Thus additional constraints are needed, i.e. Ga, Cl capture rate, B8 NC flux.

	 Weiran, Apr. 9, 2018.

*/

#ifndef JP_H
#define JP_H

#include "../../Solar/SolarNu.h"
#include "../../Target/SolarElec.h"
#include "../../Detector/FastSim.h"
#include "TMinuit.h"
#include "TH1D.h"
#include "TF1.h"
#include "TGraph.h"
#include "TFile.h"
#include <string>
#include <iostream>

using namespace std;

class JP
{
	public:
		/* Global(1) includes antineutrino experiments. Solar(others) does not. */
		/* Default setup is 2kt, 1500days. Use option!=0 for other setup. */
		JP(int generation=3, int model=1, int GlobalOrSolar=1, int option=0, int path=-1); //Model: 1 for GS and 2 for AGS.
		~JP(){};

		void SetupParameter(double sin_2_theta, double ms);

		/* If a parameterized survival probability is used, we calculate only day-night added chi2.*/
		/* This and "chi2he()" are used in band calculation. */
		void SetupParameter(TF1* SurvProb);
		void SetupParameter(TGraph* SurvProb);

		void SetupThreshold(double thres);
		/* Background flux should be in the following order: Kr85, Bi210, C11, C14, C10, Tl208, Be11, ExtTl208. */
		/* Bkgflux: 1st dimension for daynight, second for components. */
		double chi2(double** bkgflux, double* nuflux, double penalty);
		double chi2(double* nuflux); //Minimized with other pars.

		/* For high energy, background should be in the following order: C11, C10, Tl208, Be11, ExtTl208. */
		/* Neutrino flux should be in order of B8, hep. Day-night averaged. */
		double chi2he(double* bkgflux, double* nuflux, double penalty);
		/* Day-night separated. */
		double chi2he(double** bkgflux, double* nuflux, double penalty);

		/* For low energy parts. Day-night averaged. */
		double chi2le(double* bkgflux, double* nuflux, double penalty);

		double GetBkgFlux(int bkg); //Total flux in 750 live days.
		string GetBkgName(int bkg);

		/* B8 constrain from SNO NC. */
		double B8Constrain(double B8flux);

		/* Low energy constrain from solar model. */
		double ModelConstrain(double* nuflux);
		double B8ModelConstrain(double B8flux); //From model.

		TH1D* hsig[2][9]; //day and night detected spectra for 9 components. Normalized to flux=1.
		TH1D* hsigtot[9]; //day and night added, normalized to flux=1.

		double* mNuFlux;

		/* For NSI setup. */
		void SetupNSI(double el, double er, double tl, double tr);
		
	private:
		int gene;
		int mmodel;
		TFile* file; //high and low metallicity
		TH1D* hbkg[8];
		TH1D* hsim[2]; //day and night.
		TH1D* hsimtot; //day and night added.
		string bkgname[8];
		double msin_2_theta, mms;
		double threshold;
		double bkgflux[8];
		double ScaleFactor;
		double EngUp, EngLow; //Energy range for high energy chi2 calculation.
		SolarElec* elec;
		FastSim sim;
		SolarNu solar;

		static JP* mJP;
		static void jpfcn(int &npar, double* gin, double &f, double* par, int iflag);

};

#endif



