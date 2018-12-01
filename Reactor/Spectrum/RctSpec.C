/*
  Standard reactor neutrino spectrum

  Zhe Wang, Oct. 29, 2014 Prague
*/

#include "RctSpec.h"

RctSpec::RctSpec()
{
  fFractionU235=0.58;
  fFractionU238=0.07;
  fFractionPu239=0.30;
  fFractionPu241=0.05;

  fHSFluxParU235[0]=3.519;
  fHSFluxParU235[1]=-3.517;
  fHSFluxParU235[2]=1.595;
  fHSFluxParU235[3]=-4.171e-1;
  fHSFluxParU235[4]=5.004e-2;
  fHSFluxParU235[5]=-2.303e-3;
  fHSFluxParU238[0]=0.976;
  fHSFluxParU238[1]=-0.162;
  fHSFluxParU238[2]=-0.0790;
  fHSFluxParU238[3]=0;
  fHSFluxParU238[4]=0;
  fHSFluxParU238[5]=0;
  fHSFluxParPu239[0]=2.560;
  fHSFluxParPu239[1]=-2.654;
  fHSFluxParPu239[2]=1.256;
  fHSFluxParPu239[3]=-3.617e-1;
  fHSFluxParPu239[4]=4.547e-2;
  fHSFluxParPu239[5]=-2.143e-3;
  fHSFluxParPu241[0]=1.487;
  fHSFluxParPu241[1]=-1.038;
  fHSFluxParPu241[2]=4.130e-1;
  fHSFluxParPu241[3]=-1.423e-1;
  fHSFluxParPu241[4]=1.866e-2;
  fHSFluxParPu241[5]=-9.229e-4;

  fFissionsPerSecondGW = 3.1e19; // fissions/(MeV*second*GW)

  dNdE_Expected = new TF1("dNdE_Expected","([0]*(exp(pol5(1))) + [7]*(exp(pol5(8))) + [14]*(exp(pol5(15))) + [21]*(exp(pol5(22))))*[28]",0.0,10.0);

  dNdE_Expected->SetParameter(0,fFractionU235);
  dNdE_Expected->SetParameter(1,fHSFluxParU235[0]);
  dNdE_Expected->SetParameter(2,fHSFluxParU235[1]);
  dNdE_Expected->SetParameter(3,fHSFluxParU235[2]);
  dNdE_Expected->SetParameter(4,fHSFluxParU235[3]);
  dNdE_Expected->SetParameter(5,fHSFluxParU235[4]);
  dNdE_Expected->SetParameter(6,fHSFluxParU235[5]);

  dNdE_Expected->SetParameter(7,fFractionU238);
  dNdE_Expected->SetParameter(8,fHSFluxParU238[0]);
  dNdE_Expected->SetParameter(9,fHSFluxParU238[1]);
  dNdE_Expected->SetParameter(10,fHSFluxParU238[2]);
  dNdE_Expected->SetParameter(11,fHSFluxParU238[3]);
  dNdE_Expected->SetParameter(12,fHSFluxParU238[4]);
  dNdE_Expected->SetParameter(13,fHSFluxParU238[5]);

  dNdE_Expected->SetParameter(14,fFractionPu239);
  dNdE_Expected->SetParameter(15,fHSFluxParPu239[0]);
  dNdE_Expected->SetParameter(16,fHSFluxParPu239[1]);
  dNdE_Expected->SetParameter(17,fHSFluxParPu239[2]);
  dNdE_Expected->SetParameter(18,fHSFluxParPu239[3]);
  dNdE_Expected->SetParameter(19,fHSFluxParPu239[4]);
  dNdE_Expected->SetParameter(20,fHSFluxParPu239[5]);

  dNdE_Expected->SetParameter(21,fFractionPu241);
  dNdE_Expected->SetParameter(22,fHSFluxParPu241[0]);
  dNdE_Expected->SetParameter(23,fHSFluxParPu241[1]);
  dNdE_Expected->SetParameter(24,fHSFluxParPu241[2]);
  dNdE_Expected->SetParameter(25,fHSFluxParPu241[3]);
  dNdE_Expected->SetParameter(26,fHSFluxParPu241[4]);
  dNdE_Expected->SetParameter(27,fHSFluxParPu241[5]);

  dNdE_Expected->SetParameter(28,fFissionsPerSecondGW);

  // Get a histogram for easy use later
  dNdE_Histo = new TH1D("RctSpec","Reactor spectrum",10000,0,10);  // 1 keV/bin
  dNdE_Histo->GetXaxis()->SetTitle("Neutrino energy [MeV]");
  dNdE_Histo->GetYaxis()->SetTitle("Fissions/MeV/second/GW");

  for(int bin=1; bin<=10000; bin++)  {
    double energy = dNdE_Histo->GetBinCenter( bin );
    double flux = dNdE_Expected->Eval( energy );
    dNdE_Histo->SetBinContent( bin, flux );
  }
}

RctSpec::~RctSpec()
{}  

TF1* RctSpec::Get()
{
  return dNdE_Expected;
}

TH1D* RctSpec::GetSpec()
{
  return dNdE_Histo;
}

