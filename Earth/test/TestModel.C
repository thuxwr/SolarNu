/*
	 Test model effect for DN asymmetry.

*/

#include "../PREM.h"
#include "../Earth.h"
#include "../../Oscillation/Osci.h"
#include <iostream>
#include <complex>
#include "Eigen/Dense"
#include "unsupported/Eigen/MatrixFunctions"
#include "TMath.h"

using namespace Eigen;

MatrixXcd CalculateEarthRotate(double Energy, double ne, double nn, double L);
Oscillation osci(3);

int main()
{
	Earth earth;
	PREM prem;

	double sin_2_t12 = 0.1;
	double ms12 = 2e-5;
	/* Simply test a muon neutrino passing through the earth. */
	Vector3cd Psi(0,1,0);
	osci.SetupParameters(sin_2_t12, 0.0241, 0.5, 0, ms12);

	double cos_zenith[5] = {-0.2,-0.4,-0.6,-0.8,-1};
	double Energy = 10;
	for(int i=0; i<5; i++)
	{
		std::cout << "cos zenith:  " << cos_zenith[i] << std::endl;
		//std::cout << "Simple earth model calculation:" << std::endl;
		double L[3000], edens[3000], ndens[3000];
		int nseg;
		{
			earth.Intersect(cos_zenith[i], nseg, L, edens, ndens);
			for(int j=0; j<nseg; j++)
			{
				//std::cout << "L:" << L[j] << "  edens:" << edens[j] << "  ndens:" << ndens[j] << std::endl;
			}
			MatrixXcd RotateForGivenZenith = MatrixXcd::Identity(3,3);
			MatrixXcd RotateForGivenShell[nseg];

			for(int EarthSeg=0; EarthSeg<nseg; EarthSeg++)
			{
				if(2*EarthSeg<nseg) RotateForGivenShell[EarthSeg] = CalculateEarthRotate(Energy, edens[EarthSeg], ndens[EarthSeg], L[EarthSeg]);
				else RotateForGivenShell[EarthSeg] = RotateForGivenShell[nseg-EarthSeg-1];
				RotateForGivenZenith = RotateForGivenZenith * RotateForGivenShell[EarthSeg];
			}

			//std::cout << "Muon neutrino regeneration as e neutrino: " << pow((RotateForGivenZenith * Psi)(0).real(),2) << std::endl;
		}

		//std::cout << std::endl;
		std::cout << "PREM model calculation: " << std::endl;
		{
			prem.Intersect(cos_zenith[i], nseg, L, edens, ndens);
			for(int j=0; j<nseg; j++)
			{
				//std::cout << "L:" << L[j] << "  edens:" << edens[j] << "  ndens:" << ndens[j] << std::endl;
			}
			MatrixXcd RotateForGivenZenith = MatrixXcd::Identity(3,3);
			MatrixXcd RotateForGivenShell[nseg];

			for(int EarthSeg=0; EarthSeg<nseg; EarthSeg++)
			{
				if(2*EarthSeg<nseg) RotateForGivenShell[EarthSeg] = CalculateEarthRotate(Energy, edens[EarthSeg], ndens[EarthSeg], L[EarthSeg]);
				else RotateForGivenShell[EarthSeg] = RotateForGivenShell[nseg-EarthSeg-1];
				RotateForGivenZenith = RotateForGivenZenith * RotateForGivenShell[EarthSeg];
			}

			std::cout << "Muon neutrino regeneration as e neutrino: " << pow((RotateForGivenZenith * Psi)(0).real(),2) << std::endl;
		}
		std::cout << std::endl;

	}

	return 0;
}

MatrixXcd CalculateEarthRotate(double Energy, double ne, double nn, double L)
{
	double coeff = 5067.73094;
	std::complex<double> img(0,-1);
	MatrixXcd Ham = osci.Hamiltonian(Energy, ne, nn);
	return (img * coeff * Ham * L).exp();
}

