/*
	 Solar neutrino survival probability from inner core to detector.
	 Solar MSW effect, daynight(Earth MSW), changes of Sun-Earth distance are considered. Neutrinos from different radius are regarded incoherent.
	 Two-step jump is not considered yet and should be tested in the future.

	 Only GS model is considered here.

	 Weiran, Feb 25, 2018
*/

#ifndef SURVPROB_H
#define SURVPROB_H

#include <iostream>
#include "TMath.h"
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <complex>
#include "TGraph.h"
#include "TFile.h"
#include "TH1D.h"
#include "../../Earth/Earth.h"
#include "../../Solar/SolarNu.h"
#include "../../Oscillation/Osci.h"

using namespace Eigen;

class SurvProb
{
	public:
		SurvProb(int generation);
		~SurvProb(){};

		/* SetupProb is used to setup survival probability for day and night in numeric calculation. */
		void SetupProb(double Energy /* MeV */, double sin_2_t12 = 0.307, double sin_2_t13 = 0.0241, double ms12 = 7.54e-5 /* eV^2 */, double ms13 = 2.5e-3, double sin_2_t14 = 1e-4, double ms14 = 1e-5);
		/* Get prob from SetupProb numeric calculation. daynight: 0-day, 1-night;  flavor: 1-e, 2-mu, 3-tau, 4-sterile;  components from 1 to 9. */
		double GetProbFromCalculation(int daynight, int flavor, int comp);

	private:
		const static int model = 1; //GS98 model.
		/* Calculate and save rotate matrix to an array with dimension nsteps. */
		void SetupRotateMatrix(double Energy, double sin_2_t12, double sin_2_t13, double ms12, double ms13, double sin_2_t14, double ms14, int &nsteps, double* radius, MatrixXcd* Rotate);
		MatrixXcd CalculateRotateMatrix(double Energy, double radius, double L /* Normalized by sun radius */);
		MatrixXcd CalculateEarthRotate(double Energy, double ne, double nn, double L /* Unit:km */);
		MatrixXcd CalculateVacuumRotate(double Energy, double L /* Normalized by Sun-Earth Distance*/);

		int gene;
		double SunEarthDistance; /* MeV/eV^2 */
		double RSun; /* MeV/eV^2 */
		double DistPeriAp; /* Distance between perihelion and aphelion. Unit: AU */ 
		Oscillation* osci;
		Earth earth;
		SolarNu solar;
		double prob[2][4][9]; // First dimension 0 for day and 1 for night. Second dimension for e, mu, tau, sterile. Third dimension for 9 components.
		TH1D* Dist[9];
		string DN[2]; //Daynight

};

#endif


