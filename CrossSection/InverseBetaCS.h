// Author: Chris Jillings 10/18/2004

/********************************************************************
 * Cross section for nu_e bar + p -> n + e^+
 * Formulae taken from Vogel and Beacom, PRD 60, 053003
 ********************************************************************/

// Adopted by Zhe Wang, Nov. 6, 2011
// The way this is programed is so bad to understand.
// How can gSigmaTotal be negative sometimes?

#ifndef __InverseBetaCS__H__
#define __InverseBetaCS__H__

#include "TMath.h"
#include "TF1.h"
#include "TH1.h"
#include <iostream>
#include "TH1D.h"

using namespace std;

const Double_t gkMassProton = 938.27203;  // MeV
const Double_t gkMassNeutron = 939.56536; // MeV
const Double_t gkMassElectron = 0.51099892; // MeV

class InverseBetaCS 
{
 private:
  Double_t fEpos; // positron energy (MeV)
  Double_t fEnu;  // anti neutrino energy (MeV)
  
  //utility constants set in constructor...
  Double_t fF;    // hadronic weak interaction constants
  Double_t fG;    // hadronic weak interaction constants
  Double_t fF2;   // hadronic weak interaction constants
  Double_t fCosCab; // cosine of Caibibo angle
  Double_t fDelta;  // M_n - M_p (MeV)
  Double_t fYsq;   // y defined in text below eqn 11
  Double_t fMassEsq; // electron mass squared (MeV^2)
  Double_t fR,tauN;
  Double_t fSigma0; // eqn 8,9,10
  Double_t fF2Plus3G2;
  Double_t fF2MinusG2;

  // variables that are here to save calculation
  Double_t fE0;  //set in Ee1()
  Double_t fP0;  //set in Ee1();
  Double_t fV0;  //set in Ee1();
  Double_t fE1;  //set in DSigDCosTh

  TH1D* fHTotalCrossSection;

public:
  InverseBetaCS();
  virtual ~InverseBetaCS() {;}
  
  Double_t Ee0(Double_t aEnu) { return (aEnu - fDelta); }  // eqn 6 in Vogel/Beacom 
  Double_t Ee1(Double_t aEnu, Double_t aCosTheta); // eqn 11
  Double_t GammaTerm(Double_t aCosTheta); // eqn 13 

  // Differential crosssection dSigma/dE
  // Input neutrino energy unit is MeV
  // Output crosssection unit is 1E-42 cm^2/hydrogen
  Double_t SigmaTot(Double_t aEnu); // integration of eqn 12 by Gaussian quadrature

  // Convert neutrino spectrum to visible energy, unit: MeV.
  // NeuSpec and VisSpec must have the same binning and range
  // Fine binning, for example 0.001 MeV, is expected.
  // Large energy range is expected, for example 0-20 MeV for reactor neutrino.
  int VisIbdSpec(TH1D* NeuSpec, TH1D* VisSpec);
  
  Double_t DSigDCosTh(Double_t aEnu, Double_t aCosTheta); // eqn 12
  Double_t PromptEnergyToNeutrinoEnergy(Double_t aEprompt); // incles 2*.511
  TF1* fDifferential; //! here for easy integration of diff cross section
  
private:
  void setupTotalCrossSection();

};

///global function to allow calls from a TF1.
///This is necessary to use the gaussian quadrature method
///built into TF1
/// dsig/dcos(theta)
Double_t gDSigmaByDCosTheta(Double_t* x, Double_t*a);

///global function to allow calls from a TF1.
/// a redirect to InverseBetaCS::SigmaTot()
/// Output unit 1*10^{-42} cm^2
Double_t gSigmaTotal(Double_t* x, Double_t*a);

R__EXTERN InverseBetaCS* gInverseBeta;

#endif
