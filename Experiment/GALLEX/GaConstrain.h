/*
	 A simple program to calculate constrain from all radiative experiments using Gallium.

	 Weiran, Apr. 4, 2018.

*/

#ifndef GACONSTRAIN_H
#define GACONSTRAIN_H

#include "../../Target/GaCapture.h"
#include "TMath.h"

using namespace std;

class GaConstrain
{
	public:
		GaConstrain(int generation=3)
		{
			gene = generation;
			GaCapt = new GaCapture(gene);
			GaCenter = 66.1;
			GaError = 3.1;
		}

		~GaConstrain()
		{
			delete GaCapt;
		}

		void SetupParameter(double sin_2_theta, double ms)
		{
			GaCapt->SetCapture(sin_2_theta, ms);
		}

		double chi2(double* flux) //Flux in units: *10^10 cm^-2 s^-1.
		{
			double Pred = GaCapt->GetFlux(flux);
			double PredErr = GaCapt->GetError(flux);
			double Err = sqrt(PredErr*PredErr + GaError*GaError);
			return pow((Pred-GaCenter)/Err,2);
		}


	private:
		int gene;
		GaCapture* GaCapt;
		double GaCenter, GaError;

};

#endif
