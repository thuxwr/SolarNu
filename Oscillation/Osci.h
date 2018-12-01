/*
	 A general program to calculate state vector after oscillation. 
	 Hamiltonian can be achieved here.

	 2, 3 and 4 generations of neutrinos are supported.

	 Add NSI parameters.

	 Weiran, Feb 24, 2018
*/
#ifndef OSCILLATION_H
#define OSCILLATION_H

#include <Eigen/Dense>
#include <complex>

using namespace Eigen;

class Oscillation
{
	public:
		Oscillation(int generation);
		~Oscillation(){};

		void SetupParameters(double sin_2_t12 = 0.307, double sin_2_t13 = 0.0241, double sin_2_t23 = 0.5, 
												 double sin_2_t14 = 0, double ms12 = 7.54e-5, double ms13 = 2.5e-3, double ms14 = 0, double cp_phase = 0);
		/* NSI parameters: 1st and 2nd for flavor, 3rd for fermion(e,u,d), 4th for left and right. */
		void SetupNSI(double eeel, double eeer, double etel, double eter, double eeu, double eed, double euu, double eud, double etu, double etd, double uuu, double uud, double utu, double utd, double ttu, double ttd);
		void SetupNSI(double sin_2_t12, double ms12, int fermion /* 0~2 for e,u,d */, double ed, double en);

		double AdiabaticPropagation(double Energy, double ne, double nn=0, double np=0); //Neglect interference term.

		MatrixXcd GetVacuumMixingMatrix();
		MatrixXcd GetMatterMixingMatrix(double Energy, double ne, double nn=0, double np=0);

		MatrixXcd Hamiltonian(double Energy /*MeV*/, double ne /*cm^{-3}*/, double nn=0, double np=0);

		/* Adiabatic survival probability, no interference. */
		double SurvProb(double Energy, double ne, double nn=0);

		/* Get mixing angle in matter. All mixing angles are designated in range 0~pi/2. */
		/* See note for details. */
		double GetTheta12(double Energy, double ne, double nn=0, double np=0);
		double GetTheta13(double Energy, double ne, double nn=0, double np=0);
		double GetTheta14(double Energy, double ne, double nn=0);
		
		/* Get mass splitting in matter. */
		double GetMass12(double Energy, double ne, double nn=0, double np=0);
		double GetMass13(double Energy, double ne, double nn=0, double np=0);

	private:
		int gene;
		MatrixXcd Uv;
		MatrixXcd A; //A = Uv * M * Uv^{dagger}
		MatrixXcd Mv; //Mass matrix in vacuum.
		double s12, s13, s23, s14, c12, c13, c23, c14;
		double EEe, EEu, EEd, EUu, EUd, ETe, ETu, ETd, UUu, UUd, UTu, UTd, TTu, TTd;
		ComplexEigenSolver<MatrixXcd> solver;

		/* Next parameters for NSI. */
		int IsNSI;
		int Fermion;
		double Ed, En;

};

#endif
