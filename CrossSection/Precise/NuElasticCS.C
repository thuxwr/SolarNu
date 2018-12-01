#include "NuElasticCS.h"
#include "TMath.h"
#include <iostream>

using namespace std;

NuElasticCS::NuElasticCS()
{
	m = 0.511; //MeV
	s0 = 8.8059048e-45; //cm^2
	rhoNC = 1.0126; //I didn't calculate revision according to the latest top quark mass.
	s2tw = 0.2317;
	alpha = 7.2973525698e-3; //From http://pdg.lbl.gov/2014/reviews/rpp2014-rev-phys-constants.pdf
	pi = TMath::Pi();
}

double NuElasticCS::dsdTe(double Enu, double Te, int species, int nu)
{
	double ds = 0;
	double gL = GetgL(species, Te);
	double gR = GetgR(species, Te);
	double fp = Getfp(Enu, Te);
	double fm = Getfm(Enu, Te);
	double fpm = Getfpm(Enu, Te);

	if(nu==0){}
	else // swap gL and gR
	{
		double gtemp;
		gtemp = gL; gL = gR; gR = gtemp;
	}

	if(Te <= TeMax(Enu))
	{
		double z = Te/Enu;
		double q = Enu;
		ds = s0/m * (gL*gL*(1+alpha/pi * fm) + gR*gR*pow(1-z,2)*(1+alpha/pi * fp) - gR*gL*m*z/q * (1+alpha/pi * fpm));
	}

	return ds;
}

void NuElasticCS::VisElecSpec(TH1D* NuSpec, TH1D* ElecSpec, int species, int nu)
{
  // Nue Spectrum * Cross section => true electron spectrum
  double NuEmin = NuSpec->GetXaxis()->GetXmin();
  double NuEmax = NuSpec->GetXaxis()->GetXmax();
  double NuNbin = NuSpec->GetXaxis()->GetNbins();

  double binWidth = 5e-3;   // 5 keV

	int    TeNbin = ElecSpec->GetXaxis()->GetNbins();


  for( int bin=1; bin<=NuNbin; bin++ )  
	{
    double c = NuSpec->GetBinContent( bin );
    double Enu = NuSpec->GetBinCenter( bin );

    for( int el=1; el<=TeNbin; el++ )  
		{
      double Te = ElecSpec->GetBinCenter( el );
      double curr = ElecSpec->GetBinContent( el );

      // Neutrino scatters on electron
      //            N/cm2/s/MeV * MeV                    *  [cm2/MeV]  = [1/(s MeV)]
      double Addi = c           * NuSpec->GetBinWidth(1) * dsdTe( Enu, Te, species, nu );
      ElecSpec->SetBinContent( el, curr+Addi );
    }
  }
}

double NuElasticCS::GetgL(int species, double T)
{
	double gL = 0;
	if(species==0) gL = rhoNC * (1./2. - Getk(species, T) * s2tw) -1;
	else gL = rhoNC * (1./2. - Getk(species, T) * s2tw);
	return gL;
}

double NuElasticCS::GetgR(int species, double T)
{
	return -1 * rhoNC * Getk(species, T) * s2tw;
}

double NuElasticCS::Getx(double T)
{
	return sqrt(1 + 2*m/T);
}

double NuElasticCS::GetI(double T)
{
	return 1./6. * (1./3. + (3- pow(Getx(T),2))*(1./2. * Getx(T) * log((Getx(T)+1)/(Getx(T)-1)) - 1));
}

double NuElasticCS::Getk(int species, double T)
{
	double kk=0;
	if(species==0) kk = 0.9791 + 0.0097*GetI(T);
	else kk = 0.9970 - 0.00037*GetI(T);
	return kk;
}

double NuElasticCS::Getl(double E)
{
	return sqrt(E*E-m*m);
}

double NuElasticCS::GetE(double T)
{
	return T+m;
}

double NuElasticCS::Getbeta(double T)
{
	return Getl(GetE(T))/GetE(T);
}

double NuElasticCS::L(double x)
{
	return -1 * TMath::DiLog(x);
}

double NuElasticCS::Getfp(double Enu, double T)
{
	double z = T/Enu;
	double beta = Getbeta(T);
	double E = GetE(T);
	double l = Getl(E);
	return pow(1-z,-2)*(
			(E/l * log((E+l)/m)-1.) * (pow(1-z,2)*(2*log(1-z-m/(E+l))-log(1-z)-log(z)/2. -2./3.)-(z*z*log(z)+1-z)/2)
			- pow(1-z,2)/2. * (pow(log(1-z),2)+beta*(L(1-z)-log(z)*log(1-z)))
			+ log(1-z) * (z*z/2. * log(z) + (1-z)/3. * (2*z-1./2.))
			- z*z/2. * L(1-z) - z*(1-2*z)/3 * log(z) - z*(1-z)/6.
			- beta/12. * (log(z) + (1-z) * ((115-109*z)/6))
			);
}

double NuElasticCS::Getfm(double Enu, double T)
{
	double beta = Getbeta(T);
	double E = GetE(T);
	double q = Enu;
	double l = Getl(E);
	double z = T/Enu;
	return (E/l * log((E+l)/m)-1)*(2*log(1-z-m/(E+l))-log(1-z)-1./2. * log(z) - 5./12.)
		+ 1./2. *(L(z)-L(beta))-1./2. * pow(log(1-z),2) - (11./12. + z/2.)*log(1-z)
		+ z*(log(z)+1./2. * log(2*q/m)) - (31./18. + 1./12. * log(z))*beta - 11./12. * z + z*z/24.;
}

double NuElasticCS::Getfpm(double Enu, double T)
{
	double E = GetE(T);
	double l = Getl(E);
	double z = T/Enu;
	return (E/l * log((E+l)/m)-1) * 2 * log(1-z-m/(E+l));
}

double NuElasticCS::TeMax(double Enu)
{
	return 2 * Enu * Enu / (m + 2 * Enu);
}
