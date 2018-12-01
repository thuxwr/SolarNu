/*
  Standard reactor neutrino spectrum 
  with the unit of fissions/MeV/second/GW

  Zhe Wang, Oct. 29, 2014 Prague
*/

#ifndef _RCT_SPEC_H_
#define _RCT_SPEC_H_

#include "TF1.h"
#include "TH1D.h"

class RctSpec
{
 public:
  RctSpec();
  ~RctSpec();
  
  // Function with unit of fissions/MeV/second/GW
  TF1*  Get();
  // Then [GetSpec(...)->Integral()] * [GetSpec(...)->GetBinWidth()] == Total flux
  // Unit: fissions/MeV/second/GW
  TH1D* GetSpec();
  
 private:
  double fFractionU235;
  double fFractionU238;
  double fFractionPu239;
  double fFractionPu241;

  double fHSFluxParU235[6];
  double fHSFluxParU238[6];
  double fHSFluxParPu239[6];
  double fHSFluxParPu241[6];

  double fFissionsPerSecondGW;

  TF1*  dNdE_Expected;
  TH1D* dNdE_Histo;
};

#endif // _RCT_SPEC_H_
