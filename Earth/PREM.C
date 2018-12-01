#include "PREM.h"
#include <iostream>
#include <string>
#include "TMath.h"

using namespace std;

PREM::PREM(int Land)
{
	Rearth = 6371;
	land = Land;
	layers = 100;
	Boundary[0] = 1221.5; Boundary[1] = 3480.; Boundary[2] = 5701.; Boundary[3] = 5771.; Boundary[4] = 5971.;
	Boundary[5] = 6151.; Boundary[6] = 6346.6; Boundary[7] = 6356.; Boundary[8] = 6368.; Boundary[9] = 6371.;

	string nu = getenv("neutrino");
	if(nu=="")
	{
		cout << "Environment variable 'neutrino' undefined." << endl;
		exit(0);
	}

	string jppath = nu + "/Experiment/Jinping/Simulation/result/livetime.root";
	TFile* jpfile = new TFile(jppath.c_str(),"READ");
	livetime[0] = (TH1D*)jpfile->Get("livetime");

	string snopath = nu + "/Experiment/SNO/reproduce/result/livetime.root";
	TFile* snofile = new TFile(snopath.c_str(),"READ");
	livetime[1] = (TH1D*)snofile->Get("livetime");

	string skpath = nu + "/Experiment/SuperK/data/Zenith/livetime.root";
	TFile* skfile = new TFile(skpath.c_str(),"READ");
	livetime[2] = (TH1D*)skfile->Get("livetime");
}

double PREM::GetNe(double radius)
{
	double density = GetDensity(radius);
	double ProtonRatio = 0; //Ratio of protons to nucleons, Z/A.
	if(radius<=3480) ProtonRatio = 0.468;
	else ProtonRatio = 0.497;
	return density*ProtonRatio*NA;
}

double PREM::GetNn(double radius)
{
	double density = GetDensity(radius);
	double ProtonRatio = 0;
	if(radius<=3480) ProtonRatio = 0.468;
	else ProtonRatio = 0.497;
	return density*(1-ProtonRatio)*NA;
}

double PREM::GetDensity(double radius)
{
	double density = 0;
	double x = radius/Rearth;
	if(radius<=Boundary[0]) //Inner core
		density = 13.0885 - 8.8381*x*x;
	else if(radius<=Boundary[1]) //Outer core
		density = 12.5815 - 1.2638*x - 3.6426*x*x - 5.5281*x*x*x;
	else if(radius<=Boundary[2]) //Lower mantle
		density = 7.9565 - 6.4761*x + 5.5283*x*x - 3.0807*x*x*x;
	else if(radius<=Boundary[3]) //Transition zone 1
		density = 5.3197 - 1.4836*x;
	else if(radius<=Boundary[4]) //Transition zone 2
		density = 11.2494 - 8.0298*x;
	else if(radius<=Boundary[5]) //Transition zone 3
		density = 7.1089 - 3.8045*x;
	else if(radius<=Boundary[6]) //LVZ&LID
		density = 2.6910 + 0.6924*x;
	else if(radius<=Boundary[7]) //Crust 1
		density = 2.9;
	else if(radius<=Boundary[8]) //Crust 2
		density = 2.6;
	else if(radius<=Boundary[9]) //Ocean
	{
		if(land==1) density = 2.6;
		else density = 1.02;
	}
	return density;
}

void PREM::Intersect(double cos_zenith, int &nseg, double* L, double* edens, double* ndens)
{
	/* Decide slabs in each layer. */
	int nlayers[7];
	int TotalLayers = 0;
	for(int i=0; i<7; i++)
	{
		if(i==0) nlayers[i] = (int)(Boundary[i]/Boundary[9] * layers);
		else nlayers[i] = (int)((Boundary[i]-Boundary[i-1])/Boundary[9] * layers);
		if(nlayers[i]==0) nlayers[i] = 1;
		TotalLayers += nlayers[i];
	}

	/* Crust is divided into two layers and Ocean is one layer. */
	TotalLayers += 3;

	/* Decide all boundaries. */
	double r[TotalLayers];
	int InnerLayers = 0;
	for(int i=0; i<7; i++) 
	{
		for(int lyr=0; lyr<nlayers[i]; lyr++)
		{
			if(i==0) r[InnerLayers+lyr] = Boundary[i]/nlayers[i] * (lyr+1);
			else r[InnerLayers+lyr] = Boundary[i-1] + (Boundary[i]-Boundary[i-1])/nlayers[i] * (lyr+1);
		}
		InnerLayers += nlayers[i];
	}
	r[InnerLayers] = Boundary[7]; r[InnerLayers+1] = Boundary[8]; r[InnerLayers+2] = Boundary[9];

	/* Start calculating intersect regions. */

	/* Going downwards */
	if(cos_zenith >= 0)
	{
		nseg = 1;
		L[0] = 0;
		edens[0] = 0;
		ndens[0] = 0;
		return;
	}

	/* Going upwards, intersect with earth shell. */
	double fullLength[TotalLayers];
	double ccd = Rearth * sqrt(1-cos_zenith*cos_zenith); //chord-center distance
	for(int lyr=0; lyr<TotalLayers; lyr++)
	{
		fullLength[lyr] = 0;
		if(ccd<=r[lyr]) fullLength[lyr] = 2 * sqrt(r[lyr]*r[lyr] - ccd*ccd);
	}

	nseg = 0;
	for(int lyr=0; lyr<TotalLayers-1; lyr++)
	{
		double Length = (fullLength[TotalLayers-1-lyr]-fullLength[TotalLayers-2-lyr])/2.;
		if(Length>0)
		{
			L[nseg] = Length;
			edens[nseg] = (GetNe(r[TotalLayers-1-lyr])+GetNe(r[TotalLayers-2-lyr]+1e-10))/2.;
			ndens[nseg] = (GetNn(r[TotalLayers-1-lyr])+GetNn(r[TotalLayers-2-lyr]+1e-10))/2.;
			nseg++;
		}
	}

	if(fullLength[0]<=0) L[nseg-1] *= 2;
	else
	{
		L[nseg] = fullLength[0];
		edens[nseg] = (GetNe(r[0])+GetNe(0))/2.;
		ndens[nseg] = (GetNn(r[0])+GetNn(0))/2.;
		nseg++;
	}

	int Totalseg = 2 * nseg - 1;
	while(nseg < Totalseg)
	{
		L[nseg] = L[Totalseg-nseg-1];
		edens[nseg] = edens[Totalseg-nseg-1];
		ndens[nseg] = ndens[Totalseg-nseg-1];
		nseg++;
	}
}

TH1D* PREM::GetLivetime(string experiment)
{
	if(experiment == "Jinping") return livetime[0];
	else if(experiment == "SNO") return livetime[1];
	else if(experiment == "SK") return livetime[2];
	else 
	{
		cout << "Livetime distribution for " << experiment << " is not available." << endl;
		TH1D* nullh = new TH1D;
		return nullh;
	}
}



