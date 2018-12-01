#include "Response.h"
#include <iostream>
#include "TMath.h"
#include "Math/ProbFuncMathCore.h"

using namespace std;

Response::Response()
{
	for(int i=0; i<nGauss; i++)
	{
		double x = -5 + 10. / nGauss * i;
		gauss_array[i] = ROOT::Math::gaussian_cdf(x);
	}
}

inline double Response::gauss_cdf(double x)
{
	if(x < -5) return 0;
	if(x > 5) return 1;
	return gauss_array[(int)(nGauss*(x+5)/10.)];
}

void Response::Convert(TH1D*& DetS, TH1D* TrueS, string ExpName)
{
	/* Create Detected spectrum. */
	double EMin = 0;
	double TrueEMax = TrueS->GetXaxis()->GetXmax();
	double DetEMax = TrueEMax + 5 * Resolution(TrueEMax, ExpName);
	double binWidth = TrueS->GetXaxis()->GetBinWidth(1);
	int binNum = (int)(floor(DetEMax/binWidth))+1;
	double EMax = binNum * binWidth + EMin;
	string name = TrueS->GetName();
	name = name + "Det";

	DetS = new TH1D(name.c_str(), name.c_str(), binNum, EMin, EMax);
	for(int Tbin=1; Tbin<=TrueS->GetNbinsX(); Tbin++)
	{
		/* Convolution from TrueS to DetS. */
		double Energy = TrueS->GetXaxis()->GetBinCenter(Tbin); 
		double TrueContent = TrueS->GetBinContent(Tbin);
		double res = Resolution(Energy, ExpName);
		for(int Dbin=DetS->FindBin(Energy-5*res); Dbin<=DetS->FindBin(Energy+5*res); Dbin++)
		{
			double cdf_low = gauss_cdf((DetS->GetXaxis()->GetBinLowEdge(Dbin) - Energy)/res);
			double cdf_up  = gauss_cdf((DetS->GetXaxis()->GetBinUpEdge(Dbin)  - Energy)/res);
			double content = (cdf_up - cdf_low) * TrueContent;
			DetS->Fill(DetS->GetBinCenter(Dbin), content);
		}
	}
}

double Response::Resolution(double Energy, string ExpName)
{
	if(ExpName=="SKI")        return 0.2468  + 0.1492*sqrt(Energy) + 0.0690*Energy;
	else if(ExpName=="SKII")  return 0.0536  + 0.5200*sqrt(Energy) + 0.0458*Energy;
	else if(ExpName=="SKIII") return -0.123  + 0.376*sqrt(Energy)  + 0.0349*Energy;
	else if(ExpName=="SKIV")  return -0.0839 + 0.349*sqrt(Energy)  + 0.0397*Energy;
	else cout << "The experiment " << ExpName << " has not been supported." << endl;
	return -1e8; //To make an error deliberatedly.
}
