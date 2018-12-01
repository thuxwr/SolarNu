#include "../SolarElec.h"
#include "TH1D.h"
#include "TFile.h"

using namespace std;

int main()
{
	SolarElec elec(3);
	TH1D* day = elec.GetElecSpec(0.307,7.5e-5,0,1,5.25e-4);
	TH1D* night = elec.GetElecSpec(0.307,7.5e-5,1,1,5.25e-4);
	TFile* file = new TFile("test.root","RECREATE");
	day->SetName("day");
	night->SetName("night");
	day->Write();
	night->Write();
	file->Close();
	return 0;
}
