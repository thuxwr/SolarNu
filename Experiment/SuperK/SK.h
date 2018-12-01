/*
	 Get systematic uncertainties, detected rates and others for solar neutrino experiments in SuperK.

	 Weiran, Mar. 22, 2018.

*/

#ifndef SK_H
#define SK_H

#include <string>
#include "TH1D.h"

using namespace std;

class SK
{
	public:
		SK();
		~SK(){};

		/* Spectrum systematic error(relative). Energy here means kinetic energy. */
		double GetCorreSysError(double energy, int source, int phase, int sign); //source 1~3: B8 spectrum, energy scale and energy resolution. phase 1~4: SK-I~SK-IV. sign: 0:+, 1:-.
		double GetUncorreSysError(double energy, int phase, int sign);
		double GetUncorreSysError(double energy, int phase); //Symmetric error.

		/* Detectable rate. */
		double GetRate(int bin, int daynight, int phase); //daynight: 0-day, 1-night, 2-all.
		double GetStatError(int bin, int daynight, int phase, int sign);
		double GetStatError(int bin, int daynight, int phase); //Symmetric error.

		/* Rate from MC, assuming no oscillation. */
		double GetRateMC(int bin, int phase);
		double GetRateMCB8(int bin, int phase);
		double GetRateMChep(int bin, int phase);

		/* Flux systematic error in percentage, except for Energy scale, resolution and spectrum uncertainty. Affect flux in total. */
		double GetFluxSysError(int phase);
		/* Systematic error for daynight and seasonal variation. Unfinished. */

		/* Output for ROOT to setup histogram. */
		int GetNbinsX(int phase); //x-axis: electron kinetic energy.
		double* GetXbins(int phase);

		/* Get SK total spectrum. I rescale all spectrum and make MC flux the same. */
		/* MC flux(same as SNO): 5.25e6 for B8, 7.88e3 for hep. */
		TH1D* GetSpec(int phase, int daynight=2); //daynight: 0-day, 1-night, 2-all

	private:
		string SKpath;
		string PhaseName[4];

		int nbinsCorre[4], nbinsUncorre[4], nbinsRate[4]; 

		/* For correlated systematic error. */
		double CorreEnergyLowEdge[4][40]; //4 phases, energy bin number is smaller than 40.
		double CorreEnergyHighEdge[4][40];

		double SpecErrorHigh[4][40]; //Correlated systematic error for B8 spectrum.
		double SpecErrorLow[4][40];

		double ScaleErrorHigh[4][40]; //Correlated systematic error for energy scale.
		double ScaleErrorLow[4][40];

		double ResoErrorHigh[4][40]; //Correlated systematic error for energy resolution.
		double ResoErrorLow[4][40];

		/* For uncorrelated systematic error. */
		double UncorreEnergyLowEdge[4][40];
		double UncorreEnergyHighEdge[4][40];

		double UncorreErrorHigh[4][40];
		double UncorreErrorLow[4][40];

		/* For detected ES rate and error. */
		double RateEnergyHighEdge[4][40];
		double RateEnergyLowEdge[4][40];

		double ObservedRateAll[4][40];
		double RateAllErrorHigh[4][40];
		double RateAllErrorLow[4][40];

		double ObservedRateDay[4][40];
		double RateDayErrorHigh[4][40];
		double RateDayErrorLow[4][40];

		double ObservedRateNight[4][40];
		double RateNightErrorHigh[4][40];
		double RateNightErrorLow[4][40];

		double ExpectedB8Rate[4][40];
		double ExpectedhepRate[4][40];


};

#endif
