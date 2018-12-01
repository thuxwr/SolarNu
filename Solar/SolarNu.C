#include "SolarNu.h"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include "TGraph.h"
#include "TMath.h"

using namespace std;

SolarNu::SolarNu()
{
	string nu = getenv("neutrino");
	if(nu=="") 
	{
		cout << "Environment variable 'neutrino' undefined." << endl;
		exit(0);
	}

	/* Set Model Name and Component Name. */
	ModelName[0] = "GS98";
	ModelName[1] = "AGS09";
	Model[0] = "GS";
	Model[1] = "AGS";

	CompName[0] = "pp";
	CompName[1] = "pep";
	CompName[2] = "hep";
	CompName[3] = "Be7_384";
	CompName[4] = "Be7_862";
	CompName[5] = "B8";
	CompName[6] = "N13";
	CompName[7] = "O15";
	CompName[8] = "F17";

	/* Get Flux and Error for each model. */
	ifstream infile((nu+"/Solar/data/flux/Flux.dat").c_str());
	string line;
	int mode=0;
	while(getline(infile, line))
	{
		if(line[0]=='#') continue;
		istringstream iss(line);
		iss >> Flux[mode][0] >> Error[mode][0] >> Flux[mode][1] >> Error[mode][1] >> Flux[mode][2] >> Error[mode][2] >> Flux[mode][3] >> Error[mode][3] >> Flux[mode][5] >> Error[mode][5] >> Flux[mode][6] >> Error[mode][6] >> Flux[mode][7] >> Error[mode][7] >> Flux[mode][8] >> Error[mode][8];

		Flux[mode][4] = Flux[mode][3] * (1-Be7_384Ratio);
		Error[mode][4] = Error[mode][3];
		Flux[mode][3] = Flux[mode][3] * Be7_384Ratio;

		mode++;
	}
	
	/* Prepare neutrino spectrum, flux distribution and number density distribution graph. */
	NuSpec = new TFile((nu + "/Solar/data/spectrum/NuSpec.root").c_str(),"READ");
	Density = new TFile((nu + "/Solar/data/distribution/Density.root").c_str(),"READ");

	EDensity = (TGraph*)Density->Get("EDensity");
	NDensity = (TGraph*)Density->Get("NDensity");
	PDensity = (TGraph*)Density->Get("PDensity");

	Dist[0] = new TFile((nu + "/Solar/data/distribution/Dist-GS.root").c_str(), "READ");
	Dist[1] = new TFile((nu + "/Solar/data/distribution/Dist-AGS.root").c_str(), "READ");

	for(int model=1; model<=2; model++)
	{
		for(int comp=1; comp<=9; comp++)
		{
			TH1D* temp = (TH1D*)Dist[model-1]->Get((Model[model-1]+CompName[comp-1]+"Dist").c_str());
			Distribution[model-1][comp-1] = new TH1D(*temp); 
		}
	}

	BinnedDist[0] = new TFile((nu + "/Solar/data/distribution/binned/Dist-GS.root").c_str(), "READ");
	BinnedDist[1] = new TFile((nu + "/Solar/data/distribution/binned/Dist-AGS.root").c_str(), "READ");

	for(int model=1; model<=2; model++)
	{
		for(int comp=1; comp<=9; comp++)
		{
			TH1D* temp = (TH1D*)BinnedDist[model-1]->Get((Model[model-1] + CompName[comp-1] + "Dist").c_str());
			BinnedDistribution[model-1][comp-1] = new TH1D(*temp); 
		}
	}

	Be7_384Ratio = 0.1052;
	NA = 6.022141e23;
}

SolarNu::~SolarNu()
{
	for(int model=1; model<=2; model++)
	{
		Dist[model-1]->Close();
		BinnedDist[model-1]->Close();
		delete Dist[model-1]; delete BinnedDist[model-1];
	}
	NuSpec->Close(); Density->Close(); 
	delete NuSpec; delete EDensity; delete NDensity;
}

	

double SolarNu::GetFlux(int model, int comp)
{
	return Flux[model+2][comp-1];
}

double SolarNu::GetError(int model, int comp)
{
	return Error[model+2][comp-1] * Flux[model+2][comp-1];
}

double SolarNu::GetModelSeparation(int comp)
{
	return abs(Flux[3][comp-1]-Flux[4][comp-1]);
}

TH1D* SolarNu::GetSpec(int model, int comp)
{
	TH1D* spec = (TH1D*) NuSpec->Get((ModelName[model-1] + CompName[comp-1]).c_str());
	return spec;
}

TH1D* SolarNu::GetSpec(int comp, double flux, int option)
{
	if(option==310) cout << "Happy birthday!" << endl;
	TH1D* spec = GetSpec(1, comp);
	spec->Scale(1/spec->Integral("width"));
	spec->Scale(flux);
	return spec;
}


string SolarNu::GetModelName(int model)
{
	return ModelName[model-1];
}

string SolarNu::GetCompName(int comp)
{
	return CompName[comp-1];
}

double SolarNu::GetEDensity(double radius)
{
	return pow(10,EDensity->Eval(radius)) * NA;
}

double SolarNu::GetNDensity(double radius)
{
	return pow(10,NDensity->Eval(radius)) * NA;
}

double SolarNu::GetPDensity(double radius)
{
	return pow(10,PDensity->Eval(radius)) * NA;
}

TH1D* SolarNu::GetFluxDist(int model, int comp)
{
	return Distribution[model-1][comp-1];
}

TH1D* SolarNu::GetBinnedFluxDist(int model, int comp)
{
	return BinnedDistribution[model-1][comp-1];
}

double SolarNu::GetResonantDensity(double Energy, double sin_2_theta, double ms)
{
	double a = 2.5348e-31 /* 2sqrt(2)G_{F} */;
	double cos_2theta = 1 - 2*sin_2_theta;
	double ne = ms * cos_2theta / a / Energy;
	return ne;
}

double SolarNu::ModelConstrain(int model, double* nuflux)
{
	double chi2 = 0;
	for(int comp=1; comp<=9; comp++)
	{
		if(comp==4) continue; //Be7_384
		if(comp==5) continue; //Be7_862
		double ModelFlux = GetFlux(model, comp);
		double Error = TMath::Abs(GetFlux(1, comp)-GetFlux(2, comp));
		chi2 += pow((nuflux[comp-1]-ModelFlux)/Error, 2);
	}
	/* Deal with Be7. */
	double Be7Flux = GetFlux(model, 4) + GetFlux(model, 5);
	double Be7Error = TMath::Abs(GetFlux(1,4)+GetFlux(1,5)-GetFlux(2,4)-GetFlux(2,5));
	double Be7Det = nuflux[3]+nuflux[4];
	chi2 += pow((Be7Flux-Be7Det)/Be7Error, 2);
	return chi2;
}


