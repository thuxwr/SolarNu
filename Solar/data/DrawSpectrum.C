/*
	 Draw solar neutrino spectra according to their shape and flux.
	
	 Weiran, Feb 25, 2018
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include "TH1D.h"
#include <string>
#include "TFile.h"

using namespace std;

const double Be7_384Ratio = 0.1052;
string ModelName[2]={"GS98","AGS09"};
string CompName[9]={"pp","pep","hep","Be7_384","Be7_862","B8","N13","O15","F17"};

void ReadFile(TH1D*& hist, const char* file);

int main()
{
	ifstream infile("./flux/Flux.dat");
	string line;
	double Flux[5][9];
	double Error[5][9];
	int mode=0;
	while(getline(infile, line))
	{
		char f = line[0];
		if(f=='#') continue;
		istringstream iss(line);
		iss >> Flux[mode][0] >> Error[mode][0] >> Flux[mode][1] >> Error[mode][1] >> Flux[mode][2] >> Error[mode][2] >> Flux[mode][3] >> Error[mode][3] >> Flux[mode][5] >> Error[mode][5] >> Flux[mode][6] >> Error[mode][6] >> Flux[mode][7] >> Error[mode][7] >> Flux[mode][8] >> Error[mode][8];

		/* Be7 384 and 862 lines */
		Flux[mode][4] = Flux[mode][3] * (1-Be7_384Ratio);
		Error[mode][4] = Error[mode][3];
		Flux[mode][3] = Flux[mode][3] * Be7_384Ratio;

		mode++;
	}

	TH1D* Spectrum[2][9]; //Only GS98 and AGS09 is drawn.
	ReadFile(Spectrum[0][0], "./spectrum/pp.dat");
	ReadFile(Spectrum[0][1], "./spectrum/pep.dat");
	ReadFile(Spectrum[0][2], "./spectrum/hep.dat");
	ReadFile(Spectrum[0][3], "./spectrum/Be7-384.dat");
	ReadFile(Spectrum[0][4], "./spectrum/Be7-862.dat");
	ReadFile(Spectrum[0][5], "./spectrum/B8.dat");
	ReadFile(Spectrum[0][6], "./spectrum/N13.dat");
	ReadFile(Spectrum[0][7], "./spectrum/O15.dat");
	ReadFile(Spectrum[0][8], "./spectrum/F17.dat");

	for(int comp=0; comp<9; comp++) Spectrum[1][comp] = new TH1D(*Spectrum[0][comp]);

	/* Rescale */
	for(int mod=3; mod<5; mod++)
	{
		for(int comp=0; comp<9; comp++)
		{
			double scale = 1/Spectrum[mod-3][comp]->GetXaxis()->GetBinWidth(1);
			Spectrum[mod-3][comp]->Scale(scale);
			Spectrum[mod-3][comp]->Scale(Flux[mod][comp]);
			string name = ModelName[mod-3] + CompName[comp];
			Spectrum[mod-3][comp]->SetTitle(name.c_str());
			Spectrum[mod-3][comp]->SetName(name.c_str());

			Spectrum[mod-3][comp]->GetXaxis()->SetTitle("Neutrino Energy [MeV]");
			Spectrum[mod-3][comp]->GetYaxis()->SetTitle("Flux [#times 10^{10} /(cm^{2} s MeV)]");
		}
	}
			
	TFile* file = new TFile("./spectrum/NuSpec.root","RECREATE");
	for(int mod=3; mod<5; mod++)
	{
		for(int comp=0; comp<9; comp++)
		{
			Spectrum[mod-3][comp]->Write();
		}
	}
	file->Close();

	return 0;
}

/* Read spectrum and normalize. */
void ReadFile(TH1D*& hist, const char* file)
{
	const int nline = 5000;
	double E[nline], Spec[nline];

	ifstream infile(file);
	string line;
	int n=0;

	while(getline(infile, line))
	{
		istringstream iss(line);
		iss >> E[n] >> Spec[n];
		n++;
	}

	int nbin = n;
	double min = E[0];
	double max = E[nbin-1];
	double width = (max-min)/(nbin-1);

	min = min - width/2.;
	max = max + width/2.;

	hist = new TH1D(file, file, nbin, min, max);
	for(int i=1; i<=nbin; i++)
	{
		hist->SetBinContent(i, Spec[i-1]);
	}

	hist->Scale(1/hist->Integral());
}
