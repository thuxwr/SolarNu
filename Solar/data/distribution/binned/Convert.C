/*
	 Convert flux distribution to binned distribution. This is used in numeric calculation.
	 Total bin number is about 120.

	 Weiran, Feb 26, 2018
*/

#include "../../../SolarNu.h"
#include "TH1D.h"
#include "TFile.h"
#include <string>

using namespace std;

string Model[2] = {"GS", "AGS"};

int main()
{
	SolarNu solar;
	for(int model=1; model<=2; model++)
	{
		TFile* file = new TFile(("Dist-" + Model[model-1] + ".root").c_str(),"RECREATE");
		for(int comp=1; comp<=9; comp++)
		{
			TH1D* hist = solar.GetFluxDist(model, comp);
			hist->Rebin(10);
			hist->Scale(1/(hist->Integral("width")));
			hist->Write();
		}
		file->Close();
	}
	return 0;
}
