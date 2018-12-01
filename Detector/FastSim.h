/*
  Fast detector response simulation.
	Neglect those events outside 5sigma.
  
  Weiran, Feb 24, 2018
*/

#ifndef _FAST_SIM_H_
#define _FAST_SIM_H_

#include "TH1D.h"
#include "Math/ProbFuncMathCore.h"

class FastSim
{
	public:
  	FastSim();
  	~FastSim(){};

  	char* GetModeName(int mode ) 
		{
    	return ModeName[mode-1];
  	};

  	// Convert a true spectrum to a detected spectrum
  	// All input spectra use MeV for x-axis
  	void Convert(TH1D*& DetS, TH1D* TrueS, int mode = 2 /* 1-4, for 1000,500,200,50PE/MeV */);

  	// Return resolution just by light yield
  	double Resolution(double Energy, int mode = 2);

	private:
		const static int nGauss = 1000000;
		double gauss_cdf(double x);
		double gauss_array[nGauss];
  	char ModeName[5][10];
  	double Yield[5];   // Light yield of 5 detector response modals
};

#endif 
