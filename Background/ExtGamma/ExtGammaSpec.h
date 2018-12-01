//
// External gamma spectrum
// Zhe Wang Created Sep. 22, 2014
//
// Only work for Tl208 2.6 MeV signal
//
#ifndef EXTGAMMASPEC_H
#define EXTGAMMASPEC_H

#include <string>
#include "TH1D.h"

#define MAXGAMMA  40

class TH1D;
class TGraph;

class ExtGammaSpec {
 public:
  ExtGammaSpec(std::string name="ExtTl208", double atten=0.4 /* MeV */);
  virtual ~ExtGammaSpec();
    
  void SetAttenCoeff(double atten);
  void SetIsotope(std::string name);
  void MakeSpectra();
  
  std::string isotope;
  double eGamma;  // Only the most significant gamma now, otherwise very complicated.
  
  // binwidth 1keV, Integral() == 1
  TH1D* hVisTotal;
  
  double AttenuCoeff;  /* Unit: MeV */
};

#endif  // EXTGAMMASPEC_H
