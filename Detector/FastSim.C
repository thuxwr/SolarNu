#include "FastSim.h"
#include "TMath.h"
#include <string>

using namespace std;

FastSim::FastSim()
{
	strcpy(ModeName[0], "1000PE");
	strcpy(ModeName[1], "500PE");
	strcpy(ModeName[2], "200PE");
	strcpy(ModeName[3], "50PE");
	strcpy(ModeName[4], "10PE");

	Yield[0] = 1000;
	Yield[1] = 500; /* PhysRevD.89.112007 (2014) */
	Yield[2] = 200;
	Yield[3] = 50;
	Yield[4] = 10;

	/* Calculate gaussian_cdf in +-5sigma and save to array. */
	for(int i=0; i<nGauss; i++)
	{
		double x = -5 + 10. / nGauss * i;
		gauss_array[i] = ROOT::Math::gaussian_cdf(x);
	}
}

inline double FastSim::gauss_cdf(double x)
{
	if(x < -5) return 0;
	if(x > 5) return 1;
	return gauss_array[(int)(nGauss*(x+5)/10.)];
}

void FastSim::Convert(TH1D*& DetS, TH1D* TrueS, int mode /* 1-5 for differend yields */)
{
	/* Create Detected spectrum. */
	double EMin = 0;
	double TrueEMax = TrueS->GetXaxis()->GetXmax();
	double DetEMax = TrueEMax + 5 * Resolution(TrueEMax, mode);
	double binWidth = TrueS->GetXaxis()->GetBinWidth(1);
	int binNum = (int)(floor(DetEMax/binWidth))+1;
	double EMax = binNum * binWidth + EMin;
	string name = TrueS->GetName();
	name = name + "Det";

	DetS = new TH1D("", "", binNum, EMin, EMax);
	for(int Tbin=1; Tbin<=TrueS->GetNbinsX(); Tbin++)
	{
		/* Convolution from TrueS to DetS. */
		double Energy = TrueS->GetXaxis()->GetBinCenter(Tbin); 
		double TrueContent = TrueS->GetBinContent(Tbin);
		double res = Resolution(Energy, mode);
		for(int Dbin=DetS->FindBin(Energy-5*res); Dbin<=DetS->FindBin(Energy+5*res); Dbin++)
		{
			double cdf_low = gauss_cdf((DetS->GetXaxis()->GetBinLowEdge(Dbin) - Energy)/res);
			double cdf_up  = gauss_cdf((DetS->GetXaxis()->GetBinUpEdge(Dbin)  - Energy)/res);
			double content = (cdf_up - cdf_low) * TrueContent;
			DetS->Fill(DetS->GetBinCenter(Dbin), content);
		}
	}
}

double FastSim::Resolution(double Energy, int mode)
{
	return Energy/sqrt(Energy*Yield[mode-1]);
}
