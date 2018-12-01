#include "NuElasticCS.h"
#include <math.h>
#include <iostream>
using namespace std;

NuElasticCS::NuElasticCS() 
{
  s0 = 88.06e-46;
  me = 0.511;
  g1e = 0.73;
  g2e = 0.23;
  g1mutau = -0.27;
  g2mutau = 0.23;
	IsNSI = false;
}

double NuElasticCS::dsdTe( double Enu, double Te, int species /* default =0; 0,1,2 for nu_e, mu, tau */, int nu /* 0: nu; 1: anti-nu */)
{
  double ds=0;
  double g1, g2, g;
  
  if( species==0 )  {
    g1 = g1e;
    g2 = g2e;
  } else {
    g1 = g1mutau;
    g2 = g2mutau;
  }

	if(IsNSI && species==0)  {
		g1 += eL;
		g2 += eR;
	} else if(IsNSI && species==2)  {
		g1 += tL;
		g2 += tR;
	}


  if( nu==0 )  {
  } else {
    // swap them
    g=g1; g1=g2; g2=g;
  }

  if( Te<=TeMax( Enu ) ) {
    ds = s0/me * ( g1*g1 + g2*g2*pow((1-Te/Enu),2) - g1*g2*me*Te/(Enu*Enu) );
  }
  
  return ds;
}

double NuElasticCS::s( double Enu, int species /* default =0; 0,1,2 for nu_e, mu, tau */, int nu /* 0: nu; 1: anti-nu */)
{
  double g1, g2, g;

  if( species==0 )  {
    g1 = g1e;
    g2 = g2e;
  } else {
    g1 = g1mutau;
    g2 = g2mutau;
  }

  if( nu==0 )  {
  } else {
    // swap them
    g=g1; g1=g2; g2=g;
  }

  
  return s0/me * ( (g1*g1+g2*g2)*TeMax(Enu) 
		   - (g2*g2+g1*g2*me/(2*Enu))*pow(TeMax(Enu),2)/Enu
		   + 1/3.*g2*g2*(pow(TeMax(Enu),3)/pow(Enu,2)) );
}

double NuElasticCS::TeMax( double Enu )
{
  return 2*Enu*Enu/(me+2*Enu);
}

// Convert neutrino spectrum to visible electron kinetic energy spectrum, unit: MeV.
// Requirement: Emin of ElecSpec is zero. Emax of ElecSpec is the same as NeuSpec.
//              BinWidth of ElecSpec is 0.001 MeV. 
int NuElasticCS::VisElecSpec(TH1D* NuSpec, TH1D* ElecSpec, int species/* 0,1,2 for nu_e, nu_mu, nu_tau */, int nu/* 0: nu; 1: anti-nu */)
{
  // Nue Spectrum * Cross section => true electron spectrum
  double NuEmin = NuSpec->GetXaxis()->GetXmin();
  double NuEmax = NuSpec->GetXaxis()->GetXmax();
  double NuNbin = NuSpec->GetXaxis()->GetNbins();

  double binWidth = 5e-3;   // 5 keV

	int    TeNbin = ElecSpec->GetXaxis()->GetNbins();


  for( int bin=1; bin<=NuNbin; bin++ )  {
    double c = NuSpec->GetBinContent( bin );
    double Enu = NuSpec->GetBinCenter( bin );

    for( int el=1; el<=TeNbin; el++ )  {
      double Te = ElecSpec->GetBinCenter( el );
      double curr = ElecSpec->GetBinContent( el );

      // Neutrino scatters on electron
      //            N/cm2/s/MeV * MeV                    *  [cm2/MeV]  = [1/(s MeV)]
      double Addi = c           * NuSpec->GetBinWidth(1) * dsdTe( Enu, Te, species, nu );
      ElecSpec->SetBinContent( el, curr+Addi );
    }
  }

  return 1;
}

void NuElasticCS::SetupNSI(double el, double er, double tl, double tr)
{
	IsNSI = true;
	eL = el;
	eR = er;
	tL = tl;
	tR = tr;
}
