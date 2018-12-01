// External gamma spectrum
// Zhe Wang Created Sep. 22, 2014
//
#include "ExtGammaSpec.h"
#include "math.h"

ExtGammaSpec::ExtGammaSpec(std::string name, double atten /* MeV */) 
{
  SetAttenCoeff(atten);
  SetIsotope(name);
  MakeSpectra();
}

ExtGammaSpec::~ExtGammaSpec()
{}
    
void ExtGammaSpec::SetAttenCoeff(double atten)
{
  AttenuCoeff=atten;
}

void ExtGammaSpec::SetIsotope(std::string name)
{
  isotope=name;
}

void ExtGammaSpec::MakeSpectra()
{
  if( isotope=="ExtTl208" )  {
    eGamma = 2.614; // MeV
  } 

  double binWidth = 5e-3;
  double Begin = 0;
  double End = eGamma;
  int nEng = int(eGamma/binWidth);

  double Ext = eGamma*1.05;
  int nbin = int(Ext/binWidth);
  
  hVisTotal = new TH1D(isotope.c_str(), isotope.c_str(), nbin, Begin, Ext);
  
  for( int b=1; b<=nEng; b++ )  {
    double e = hVisTotal->GetBinCenter(b);
    double c = exp(-(eGamma-e)/AttenuCoeff);
    hVisTotal->SetBinContent( b, c );
  }
  
  for( int b=nEng+1; b<=nbin; b++ )  {
    hVisTotal->SetBinContent( b, 0 );
  }

  
  hVisTotal->Scale( 1/hVisTotal->Integral() );
}
