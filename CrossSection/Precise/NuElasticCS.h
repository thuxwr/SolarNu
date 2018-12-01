/*
	 Nu_(e, mu, tau) + e scattering cross section.
	 The expression in LOI (code written by Zhe Wang) is a good approximation. Here I provide a precised one.

	 From astro-ph/9502003.

	 Weiran, Mar. 25, 2018.
*/

#ifndef NUELASTICCS_H
#define NUELASTICCS_H

#include "TH1D.h"

class NuElasticCS
{
	public:
		NuElasticCS();
		~NuElasticCS(){};

		double dsdTe(double Enu, double Te, int species /* 0,1,2 for nu_e, mu, tau */, int nu = 0 /* 0:nu, 1:anti-nu */);

		double TeMax(double Enu);

		void VisElecSpec(TH1D* NuSpec, TH1D* ElecSpec, int species, int nu);

	private:
		double m; //electron mass
		double alpha; //fine-structure constant
		double pi;
		double s0; //See LOI.
		double rhoNC;
		double s2tw; //sin^2(theta_W), theta_W is Weinberg angle.
		/* gL and gR for radiative correction. */
		double GetgL(int species, double T);
		double GetgR(int species, double T);
		double Getx(double T);
		double GetI(double T);
		double Getk(int species, double T);
		/* fp, fm and fpm for QED effects. */
		double Getfp(double Enu, double T); //f_{+}(z)
		double Getfm(double Enu, double T); //f_{-}(z)
		double Getfpm(double Enu, double T);
		double Getbeta(double T);
		double Getl(double E);
		double GetE(double T);
		double L(double x); //Spence function. Inverse of DiLog(x).

};

#endif
