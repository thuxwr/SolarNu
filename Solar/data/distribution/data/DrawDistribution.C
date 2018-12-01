/*
	 Draw neutrino flux distribution for GS and AGS model.
	 Since this data only provides distribution in the center half of the Sun, both electron density and neutron density cannot be drawn.

	 Weiran, Feb 25, 2018
*/

#include "../../../SolarNu.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "TFile.h"

using namespace std;
string Model[2] = {"BS05_GS", "BS05_AGS"};
string ModelName[2] = {"GS", "AGS"};

const int nline = 5000;
double radius[nline], temp[nline], edens[nline], mass[nline], be7_frac[nline], pp[nline],
			 b8[nline], n13[nline], o15[nline], f17[nline], be7[nline], pep[nline], hep[nline];

SolarNu solar;

int main()
{
	for(int model=0; model<2; model++)
	{
		ifstream infile(("./" + Model[model] + ".dat").c_str());
		string line;
		int n=0;
		while(getline(infile, line))
		{
			istringstream iss(line);
			iss >> radius[n] >> temp[n] >> edens[n] >> mass[n] >> be7_frac[n] >> pp[n] >> b8[n] >> n13[n] >> o15[n] >> f17[n] >> be7[n] >> pep[n] >> hep[n];
			n++;
		}
		int nbins = n;

		TH1D* dist[9];
		for(int comp=1; comp<=9; comp++)
		{
			string name = ModelName[model] + solar.GetCompName(comp) + "Dist";
			dist[comp-1] = new TH1D(name.c_str(), name.c_str(), nbins, 0, radius[n-1]);
			dist[comp-1]->GetXaxis()->SetTitle("radius/R_{sun}");
			dist[comp-1]->GetYaxis()->SetTitle("d(Flux)/d(radius/R_{sun})");
		}

		for(int bin=1; bin<=nbins; bin++)
		{
			dist[0]->SetBinContent(bin, pp[bin-1]);
			dist[1]->SetBinContent(bin, pep[bin-1]);
			dist[2]->SetBinContent(bin, hep[bin-1]);
			dist[3]->SetBinContent(bin, be7[bin-1]);
			dist[4]->SetBinContent(bin, be7[bin-1]);
			dist[5]->SetBinContent(bin, b8[bin-1]);
			dist[6]->SetBinContent(bin, n13[bin-1]);
			dist[7]->SetBinContent(bin, o15[bin-1]);
			dist[8]->SetBinContent(bin, f17[bin-1]);
		}

		TFile* file = new TFile(("../Dist-" + ModelName[model] + ".root").c_str(),"RECREATE");
		for(int comp=1; comp<=9; comp++)
		{
			dist[comp-1]->Scale(1/(dist[comp-1]->Integral("width")));
			dist[comp-1]->Write();
		}
		file->Close();
	}
	return 0;
}




		

