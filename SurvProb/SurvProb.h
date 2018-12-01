/*
	 Get survival probability from numerical calculation.
	 Provide a simple model using jump probability.

	 Weiran, Mar. 14, 2018.

	 Add refined results for LMA region.
	 Apr. 9, 2018.

*/
#ifndef SURVPROB_H
#define SURVPROB_H

#include "TGraph.h"
#include <string>
#include "../Solar/SolarNu.h"

using namespace std;

class SurvProb
{
	public: 
		SurvProb(int generation=3);
		~SurvProb(){};

		/* Flavorname is: e, mu, tau, and s. Daynight: 0 for day and 1 for night. Experiment: JP, SNO and SK. */
		TGraph* GetProb(double sin_2_theta, double ms, int comp, int daynight, string flavorname = "e", string experiment = "JP", int whichtheta = 0 /* For 3nu, 0:theta12, other:theta13 */);
		/* Average of day and night. */
		TGraph* GetProb(double sin_2_theta, double ms, int comp, string flavorname, string experiment = "JP", int whichtheta = 0);

		/* An analytical method to calculate survival probability in daytime. */
		/* See Fundamentals of Neutrino Physics and Astrophysics, Carlo Giunti. */
		/* Since flux distribution is insensitive to model, here fix it to GS98. */
		TGraph* GetProbAnaly(double sin_2_t12=0.307, double ms12=7.5e-5, double sin_2_t13=0.0241, double ms13=2.5e-3, int comp = 1);

		double GetProbAnaly(double Energy, double sin_2_t12, double ms12, double sin_2_t13, double ms13, int comp = 1);

	private:
		int gene;
		SolarNu solar;
		string path;
		string Path; //temporary save path
		string LMApath;
		string LMA2nupath;
		string IHpath;
		string Daynight[2];
		double GetThetaBin(double sin_2_theta);
		double GetThetaBinLMA(double sin_2_theta);
		double GetMassBin(double ms);
		double GetMassBinLMA(double ms);
		double Prob(double ne, double Energy, double sin_2_t12, double ms12, double sin_2_t13, double ms13);
};

#endif
