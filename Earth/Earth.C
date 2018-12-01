#include "Earth.h"
#include <iostream>
#include <string>

using namespace std;

Earth::Earth(int Land)
{
	/* Six regions in km. */
	r[0] = 1221;
	r[1] = 3480;
	r[2] = 5701;
	r[3] = 5971;
	r[4] = 6363;
	r[5] = 6371;
	r[6] = 6371+1e-8; //Since in allowed parameter region, oscillation in vacuum can be seen degenerated, atmosphere layer should have no effect.

	/* Corresponding Ne and Nn. */
	Ne[0] = 6.15;
	Ne[1] = 5.36;
	Ne[2] = 2.47;
	Ne[3] = 1.93;
	Ne[4] = 1.50;
	if(Land==1) Ne[5] = 1.50;
	else Ne[5] = 0.50;
	Ne[6] = 0;

	Nn[0] = 6.5;
	Nn[1] = 5.65;
	Nn[2] = 2.5;
	Nn[3] = 1.95;
	Nn[4] = 1.5;
	if(Land==1) Nn[5] = 1.5;
	else Nn[5] = 0.40;
	Nn[6] = 0;

	/* Change units for Ne and Nn to cm^-3. */
	for(int i=0; i<7; i++)
	{
		Ne[i] *= NA;
		Nn[i] *= NA;
	}

	rho[0] = 13.0;
	rho[1] = 11.3;
	rho[2] = 5.0;
	rho[3] = 3.9;
	rho[4] = 3.0;
	if(Land==1) rho[5] = 3.0;
	else rho[5] = 1.0;
	rho[6] = 0;

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

double Earth::GetNe(double radius)
{
	if(radius<r[0]) return Ne[0];
	else if(radius<r[1]) return Ne[1];
	else if(radius<r[2]) return Ne[2];
	else if(radius<r[3]) return Ne[3];
	else if(radius<r[4]) return Ne[4];
	else if(radius<r[5]) return Ne[5];
	else return Ne[6];
}

double Earth::GetNn(double radius)
{
	if(radius<r[0]) return Nn[0];
	else if(radius<r[1]) return Nn[1];
	else if(radius<r[2]) return Nn[2];
	else if(radius<r[3]) return Nn[3];
	else if(radius<r[4]) return Nn[4];
	else if(radius<r[5]) return Nn[5];
	else return Nn[6];
}

double Earth::GetDensity(double radius)
{
	if(radius<r[0]) return rho[0];
	else if(radius<r[1]) return rho[1];
	else if(radius<r[2]) return rho[2];
	else if(radius<r[3]) return rho[3];
	else if(radius<r[4]) return rho[4];
	else if(radius<r[5]) return rho[5];
	else return rho[6];
}

void Earth::Intersect(double cos_zenith, int &nseg, double* L, double* edens, double* ndens)
{
	/* Going downwards */
	if(cos_zenith >=0)
	{
		nseg = 1;
		L[0] = 0;
		edens[0] = 0;
		ndens[0] = 0;
		return;
	}

	/* Going upwards, intersect with earth shell. */
	int layer = 6;
	double Radius = r[5];
	double fullLength[layer];
	double ccd = Radius * sqrt(1-cos_zenith*cos_zenith); //chord-center distance
	for(int lyr=0; lyr<layer; lyr++)
	{
		fullLength[lyr] = 0;
		if(ccd<=r[lyr]) fullLength[lyr] = 2 * sqrt(r[lyr]*r[lyr] - ccd*ccd);
	}
	nseg = 0;
	for(int lyr=0; lyr<layer-1; lyr++)
	{
		double Length = (fullLength[layer-1-lyr]-fullLength[layer-2-lyr])/2.;
		if(Length>1e-9) //Length>0, take 1e-9.
		{
			L[nseg] = Length;
			edens[nseg] = Ne[layer-lyr-1];
			ndens[nseg] = Nn[layer-lyr-1];
			nseg++;
		}
	}
	if(fullLength[0]<=0) L[nseg-1] *= 2;
	if(fullLength[0]>0)
	{
		L[nseg] = fullLength[0];
		edens[nseg] = Ne[0];
		ndens[nseg] = Nn[0];
		nseg++;
	}

	int Totalseg = 2 * nseg - 1;
	while(nseg<Totalseg)
	{
		L[nseg] = L[Totalseg-nseg-1];
		edens[nseg] = edens[Totalseg-nseg-1];
		ndens[nseg] = ndens[Totalseg-nseg-1];
		nseg++;
	}
}

TH1D* Earth::GetLivetime(string experiment)
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

	



