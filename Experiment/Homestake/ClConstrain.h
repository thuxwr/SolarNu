/*
	 A simple program to calculate constrain from Homestake.

	 Weiran, Apr. 4, 2018.

*/

#ifndef CLCONSTRAIN_H
#define CLCONSTRAIN_H

#include "../../Target/ClCapture.h"
#include "TMath.h"

using namespace std;

class ClConstrain
{
	public:
		ClConstrain(int generation=3)
		{
			gene = generation;
			ClCapt = new ClCapture(gene);
			ClCenter = 2.56;
			ClError = 0.23;
		}

		~ClConstrain()
		{
			delete ClCapt;
		}

		void SetupParameter(double sin_2_theta, double ms)
		{
			ClCapt->SetCapture(sin_2_theta, ms);
		}

		double chi2(double* flux) //Flux in units: *10^10 cm^-2 s^-1.
		{
			double Pred = ClCapt->GetFlux(flux);
			double PredErr = ClCapt->GetError(flux);
			double Err = sqrt(PredErr*PredErr + ClError*ClError);
			return pow((Pred-ClCenter)/Err,2);
		}


	private:
		int gene;
		ClCapture* ClCapt;
		double ClCenter, ClError;

};

#endif
