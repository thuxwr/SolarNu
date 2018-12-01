/*
  Nu (e,mu,tau) e- scattering crosssection is well understood and can be 
  calculated precisely, for example, P140 of Carlo Giunti and Chung W. Kim
  Fundamentals of Neutrino Physics and Astrophysics.

  Created by Zhe Wang, May 3, 2014
  wangzhe-hep@tsinghua.edu.cn

	Add NSI effect.
	Weiran, July 2, 2018.

*/
#ifndef _NU_ELASTIC_CS_H_
#define _NU_ELASTIC_CS_H_
#include "TH1D.h"

class NuElasticCS
{
 public:
  NuElasticCS();
  ~NuElasticCS() {};

 public:
  // Differential cross-section as a function of Te in laboratory frame
  // input unit: Enu [MeV], Te [MeV]
  // output unit: ds [cm^2/MeV]
  // Eq. 5.25
  double dsdTe( double Enu, double Te, int species=0 /* 0,1,2 for nu_e, mu, tau */, int nu=0 /* 0: nu; 1: anti-nu */);
  
  // Total cross-section
  // Input unit: MeV
  // Output unit: cm^2
  // Eq. 5.32
  double s( double Enu, int species=0 /* 0,1,2 for nu_e, mu, tau */, int nu=0 /* 0: nu; 1: anti-nu */);

  // Input: MeV
  // Output: MeV
  // Eq. 5.30
  double TeMax( double Enu );

  // Convert neutrino spectrum to visible electron kinetic energy, unit: MeV.
  // Requirement: Emin of ElecSpec is zero. Emax of ElecSpec is the same as NeuSpec
  //              BinWidth of ElecSpec is 0.001 MeV.
  int VisElecSpec(TH1D* NuSpec, TH1D* ElecSpec, int species /*=0*/ /* 0,1,2 for nu_e, mu, tau */, int nu /*=0*/ /* 0: nu; 1: anti-nu */);

	void SetupNSI(double el, double er, double tl, double tr);

 private:
  double s0;
  double me;
  double g1e;
  double g2e;
  double g1mutau;
  double g2mutau;

	double eL, eR, tL, tR;
	bool IsNSI;
};

#endif // _NU_ELASTIC_CS_H_
