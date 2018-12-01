/*
	 A generator of SNO.
	 Return chi square from SNO, combining phase 1~3. Parameters are B8 flux, theta and mass.

	 Weiran, Apr. 1, 2018.

*/

#ifndef SNO_H
#define SNO_H

#include "TH1D.h"
#include "../../SurvProb/SurvProb.h"
#include "TMinuit.h"

class SNO
{
	public:
		SNO(int generation=3);
		~SNO(){};

		/* hep is regarded as background by SNO. The MCs are also B8 only. */
		void SetupParameter(double sin_2_theta, double ms); //B8flux in unit cm-2 s-1.
		double chi2(double B8flux); 
		double chi2(); //Chi2 for SNO-only. Should be irrelevant to B8flux.
		TH1D* GetSpec(); //To be finished.

		/*------------------------------------------------------------*/
		/*------Don't use the following public member functions.------*/
		TH1D* MCspec[2]; //Day and night.
		int nbins;
		SurvProb* surv;

		TH1D* MCdistspec[2]; //Distorted MC spectra from SSM & numerical calculation.
	private:
		int gene;
		int ierflg;
		double* Energy;
		double CC[2][150], ESe[2][150], ESo[2][150]; //CC spec for day and night, and ES.
		double c0, c1, c2, a0, a1; 

		static void snofcn(int &npar, double* gin, double &f, double* par, int iflag); //Minimize this function.
		static SNO* mSNO;
		TMinuit* minuit; //c0, c1, c2, a0, a1

};

#endif
