#include "RctSum.h"
#include <stdlib.h>
#include "../../MSW-Num/OsciNum.h"
using namespace std;

RctSum::RctSum()
{
  // Get $Jinping root path
  char* jp = getenv("Jinping");
  if( jp==0 ) {
    cout<<"Environment variable 'Jinping' undefined"<<endl;
    exit(0);
  }
  string Jinping = jp;
  // Then reach each data file
  //DataFile.push_back(Jinping+"/Sensitivity/Reactor/RctData/Reactor.dat");  
  DataFile.push_back(Jinping+"/Sensitivity/Reactor/RctData/ReactorsCanada.txt");
  DataFile.push_back(Jinping+"/Sensitivity/Reactor/RctData/ReactorsChina.txt");
  DataFile.push_back(Jinping+"/Sensitivity/Reactor/RctData/ReactorsEurope.txt");
  DataFile.push_back(Jinping+"/Sensitivity/Reactor/RctData/ReactorsUSA.txt");
  DataFile.push_back(Jinping+"/Sensitivity/Reactor/RctData/ReactorsOthers.txt");


  // Earth radius
  EarthRadius = 6371e3; // Meter now. 6,371,000 m = 6,371 km
  
  // Prepare the output spectra
  ExpName[0] = "Jinping";
  ExpName[1] = "Borexino";
  ExpName[2] = "SNO";
  ExpName[3] = "DayaBay";
  ExpName[4] = "JUNO";

  // Each experiment's location
  // Jinping  
  // 28¡ã11¡ä07¡åN 101¡ã37¡ä42¡åE from wiki
  ExpLatitude[0] = 28.18;
  ExpLongitude[0] = 101.62;
  // Borexino 
  // 42.46¡ãN 13.5E from wiki
  ExpLatitude[1] = 42.46;
  ExpLongitude[1] = 13.57;
  // SNO+
  // 46¡ã28¡ä00¡åN 81'W from wiki
  ExpLatitude[2] = 46.47;
  ExpLongitude[2] = -81.17;
  // DayaBay
  // 22.5953¡ãN 114.5431¡ãEfrom wiki
  ExpLatitude[3] = 22.5953;
  ExpLongitude[3] = 114.5431;
  // JUNO
  // 22¡ã34¡äN 113¡ã04for Jiangmen from wiki
  ExpLatitude[4] = 22.57;
  ExpLongitude[4] = 113.07;

  // Read in all reactor information
  ReadDataBase();

  // Set up unit reactor spectrum
  RctSpec RSpec;
  TH1D OsciSpec("OsciSpec","OsciSpec",10000,0,10);  // A temporary spec for addition of electron neutrino energy
  TH1D MuTauSpec("MuTauSpec","MuTauSpec",10000,0,10);  // A temporary spec for addition of NuMu and NuTau neutrino energy

  // Set up oscillation matrix
  OsciNum Osci;
  Osci.setup( 7.63e-5*3.0*4./1000, 4./1000, -1 ); // Set up oscillation with mean energy 0.004 GeV and Earth surface density 3.0 g/cm3.


  for(unsigned int exp=0; exp<NExp_RctSum; exp++)  {
    string title;
    title = ExpName[exp]+"Nue";
    NueSpec[exp]     = new TH1D( title.c_str(), title.c_str(), 10000,0,10);
    title = ExpName[exp]+"NonNue";
    NonNueSpec[exp]  = new TH1D( title.c_str(), title.c_str(), 10000,0,10);

    title = ExpName[exp]+"IBD";
    VisableIBD[exp]  = new TH1D( title.c_str(), title.c_str(), 10000,0,10);
    VisableIBD[exp]->GetXaxis()->SetTitle("Visible energy [MeV]");
    VisableIBD[exp]->GetYaxis()->SetTitle("Number of IBD/proton/MeV/second");

    title = ExpName[exp]+"NueElas";
    VisableNueElas[exp]    = new TH1D( title.c_str(), title.c_str(), 10000,0,10);
    VisableNueElas[exp]->GetXaxis()->SetTitle("Visible energy of Nue [MeV]");
    VisableNueElas[exp]->GetYaxis()->SetTitle("(Number of ElasScat)/electron/MeV/second");

    title = ExpName[exp]+"NonNueElas";
    VisableNonNueElas[exp] = new TH1D( title.c_str(), title.c_str(), 10000,0,10);
    VisableNonNueElas[exp]->GetXaxis()->SetTitle("Visible energy of Non-Nue [MeV]");
    VisableNonNueElas[exp]->GetYaxis()->SetTitle("(Number of ElasScat)/electron/MeV/second");

    title = ExpName[exp]+"Elas";
    VisableElas[exp]       = new TH1D( title.c_str(), title.c_str(), 10000,0,10);
    VisableElas[exp]->GetXaxis()->SetTitle("Visible energy [MeV]");
    VisableElas[exp]->GetYaxis()->SetTitle("(Number of ElasScat)/electron/MeV/second");
  }

  // Sum of all reactor background spectra
  for(unsigned int exp=0; exp<NExp_RctSum; exp++)  {
    for(unsigned int rct=0; rct<RctName.size(); rct++)  {

      double d = Distance(ExpLatitude[exp], ExpLongitude[exp], RctLatitude[rct], RctLongitude[rct], EarthRadius);
      if( d<500 ) d = 500; // Make sure all Exp are 500 meter away from reactors to avoid the input uncertainty of latitude and longitude

      for(int bin=1; bin<=10000; bin++)  {
	double energy = RSpec.GetSpec()->GetBinCenter( bin );  // MeV
	Osci.E = energy/1000.; // In principle, Osci.setup(...) should be called again. However, it is close to vacuum oscillation.
                               // Its matter effect dependence is very weak. So skip the setup now for high speed.
                               // Osci requires [GeV] for input.
	double prob = Osci.prob(0,0,d/1000.);  // Osci requires [km] for input.
	
	double RawRate = RSpec.GetSpec()->GetBinContent( bin ); // fissions/MeV/second/GW
	// 100: 1/m2 -> 1/cm2. Output: fissions/MeV/second/cm2
	OsciSpec. SetBinContent( bin, RawRate* prob   *1/(4*3.1415926*d*100*d*100)*RctPower[rct] );
	// 100: 1/m2 -> 1/cm2. Output: fissions/MeV/second/cm2
	MuTauSpec.SetBinContent( bin, RawRate*(1-prob)*1/(4*3.1415926*d*100*d*100)*RctPower[rct] );
      }
      NueSpec[exp]->Add( &OsciSpec );
      NonNueSpec[exp]->Add( &MuTauSpec );
    }
    
    // Visable spectrum vis IBD process
    InverseBetaCS IbdCs;  // cm2 is the crosssection used inside.
    // Output: N/keV/second/hydrogen
    IbdCs.VisIbdSpec( NueSpec[exp], VisableIBD[exp] );

    //VisableElas = Sum * ElasCs;
    NuElasticCS   ElasCs; // cm2 is the crosssection used inside.
    ElasCs.VisElecSpec( NueSpec[exp],    VisableNueElas[exp],    0 /* Nue */,     1 /* anti */);
    ElasCs.VisElecSpec( NonNueSpec[exp], VisableNonNueElas[exp], 1 /* NuMuTau */, 1 /* anti */);
    // Output: N/keV/second/electron
    VisableElas[exp]->Add( VisableNueElas[exp], VisableNonNueElas[exp] );
  }
}

RctSum::~RctSum()
{}

int RctSum::ReadDataBase()
{
  for(unsigned int ct=0; ct<DataFile.size(); ct++)  {
    cout<<"Processing "<<DataFile[ct]<<endl;
    Table RctTable(DataFile[ct]);
    for( int i=0; i<RctTable.NRows; i++ ) {
      RctName.push_back(      RctTable.Columns["Name"][i]);
      RctLatitude.push_back(  atof( RctTable.Columns["Latitude"][i].c_str() ) );
      RctLongitude.push_back( atof( RctTable.Columns["Longitude"][i].c_str()) );
      RctPower.push_back(     atof( RctTable.Columns["Power"][i].c_str())     );
    }
  }
  
  return 1;
}

double RctSum::Distance(double LatitudeE, double LongitudeE, double LatitudeR, double LongitudeR, double Radius)
{
  double d=0;
  double xE,yE,zE;
  double xR,yR,zR;

  // Convert all degree to radian
  LatitudeE  = LatitudeE  * 3.1415926/180;
  LongitudeE = LongitudeE * 3.1415926/180;
  LatitudeR  = LatitudeR  * 3.1415926/180;
  LongitudeR = LongitudeR * 3.1415926/180;

  xE = Radius*cos(LatitudeE)*cos(LongitudeE);
  yE = Radius*cos(LatitudeE)*sin(LongitudeE);
  zE = Radius*sin(LatitudeE);

  xR = Radius*cos(LatitudeR)*cos(LongitudeR);
  yR = Radius*cos(LatitudeR)*sin(LongitudeR);
  zR = Radius*sin(LatitudeR);

  d = sqrt((xE-xR)*(xE-xR)+(yE-yR)*(yE-yR)+(zE-zR)*(zE-zR));

  return d;
}
  
TH1D* RctSum::GetVisableIBD(int Exp)
{
  return VisableIBD[Exp];
}

TH1D* RctSum::GetVisableElas(int Exp)
{
  return VisableElas[Exp];
}
