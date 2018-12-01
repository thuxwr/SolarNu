/*
  Return the visable energy spectrum of reactor neutrinos
  for a given site. 
  
  Exp:
   0, Jinping
   1, Borexino
   2, SNO (SNO+)
   3, DayaBay
   4, JUNO

  Output unit:
   Spectrum got by GetVisableIBD has a unit of (N of IBD/second/hydrogen/MeV) with 20000 binning from 0 to 20 MeV (0.001 MeV/bin).
   Spectrum got by GetVisableElas has a unit of (N of IBD/second/electron/MeV) with 20000 binning from 0 to 20 MeV (0.001 MeV/bin). 

   Zhe Wang, Oct. 29, 2014 Prague
*/

#ifndef _RCT_SUM_H_
#define _RCT_SUM_H_

#include "TH1D.h"
#include "RctSpec.h"
#include "Table.h"
#include <vector>
#include "../../CrossSection/InverseBetaCS.h"
#include "../../CrossSection/NuElasticCS.h"

#define NExp_RctSum 5

class RctSum
{
 public:
  RctSum();
  ~RctSum();
  
  // All with output unit Events/second/hydrogen/MeV
  TH1D* GetVisableIBD(int Exp);
  // All with output unit Events/second/electron/MeV
  TH1D* GetVisableElas(int Exp);

 private:
  // Calculate the distance between (LatitueE, LongitudeE) and (LatitudeR, LongitudeR) for a sphere with radius.
  double Distance(double LatitudeE, double LongitudeE, double LatitudeR, double LongitudeR, double Radius);
  
  int ReadDataBase();

 public:
  std::string ExpName[ NExp_RctSum ];

  // Raw anti nue spectrum
  TH1D*  NueSpec[ NExp_RctSum ];
  TH1D*  NonNueSpec[ NExp_RctSum ];
  // Inverse beta decay process
  TH1D*  VisableIBD[ NExp_RctSum ];
  // Elastic scattering  
  TH1D*  VisableNueElas[ NExp_RctSum ];
  TH1D*  VisableNonNueElas[ NExp_RctSum ];
  TH1D*  VisableElas[ NExp_RctSum ];

 private:
  std::vector<std::string> DataFile;
  double ExpLatitude[ NExp_RctSum ], ExpLongitude[ NExp_RctSum ]; // Latitude and Longitude for three experimental sites
  double EarthRadius; // Earth radius 
  std::vector<std::string> RctName;  // Reactor name
  std::vector<double> RctLatitude, RctLongitude; // Latitude and Longitude for three experimental sites
  std::vector<double> RctPower;    // Reactor thermal power in GW
};

#endif // _RCT_SUM_H_
