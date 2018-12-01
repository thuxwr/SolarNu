/*
	Adiabatic calculation of survival probability for NSI. No day-night asymmetry is included.

	Weiran, July 2, 2018.
*/
#ifndef SURVNSI_H
#define SURVNSI_H

#include "../Oscillation/Osci.h"
#include "../Solar/SolarNu.h"
#include <Eigen/Dense>
#include "TMath.h"
#include "TGraph.h"
#include <complex>

using namespace Eigen;

class SurvNSI
{
	public:
		SurvNSI()
		{
			osci = new Oscillation(2);
		}
		~SurvNSI(){};

		double SurvProb(double sin_2_t12, double ms12, int comp, int fermion /* 0~2 for e,u,d */, double ed, double en, double Energy)
		{
			//osci->SetupParameters(sin_2_t12, 0.0241, 0.5, 0, ms12);
			osci->SetupNSI(sin_2_t12, ms12, fermion, ed, en);
			double radius = 0.05;
			double ne = solar.GetEDensity(radius);
			double nn = solar.GetNDensity(radius);
			double np = solar.GetPDensity(radius);
			double prob = osci->AdiabaticPropagation(Energy, ne, nn, np);
			return prob;
		}

		TGraph* GetProb(double sin_2_t12, double ms12, int comp, int fermion, double ed, double en)
		{
			TGraph* prob = new TGraph(200);
			for(int point=0; point<200; point++)
			{
				double x = 0.05 + point * 20. /200;
				double y = SurvProb(sin_2_t12, ms12, comp, fermion, ed, en, x);
				prob->SetPoint(point, x, y);
			}
			return prob;
		}


	private:
		Oscillation* osci;
		SolarNu solar;




};

#endif
