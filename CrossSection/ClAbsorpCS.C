/*
	Use "root ClAbsorpCS.C" to create "ClAbsorpCS.root" in data directory.
	 
	This is the best-estimate nue absorption cross section of chlorine, here ignore the uncertainty, as what Homestake did.
	From Bahcall et al., Phys. Rev. C 54, 411 (1996)

	Weiran, Aug 30, 2017

*/
#include <iostream>

using namespace std;

void ClAbsorpCS()
{
	double x[200], y[200],z[200];
	x[0]=1;x[1]=2;x[2]=3;x[3]=4;x[4]=5;x[5]=6;x[6]=7;x[7]=8;x[8]=9;x[9]=10;x[10]=11;x[11]=12;x[12]=13;x[13]=14;x[14]=15;
	x[15]=16;x[16]=18;x[17]=20;

	//Best-estimated cross section
	y[0]=5.21;y[1]=3.7e1;y[2]=1.02e2;y[3]=2.23e2;y[4]=5.38e2;y[5]=1.44e3;y[6]=4.62e3;y[7]=1.01e4;y[8]=1.85e4;
	y[9]=3e4;y[10]=4.45e4;y[11]=6.21e4;y[12]=8.27e4;y[13]=1.06e5;y[14]=1.33e5;y[15]=1.62e5;y[16]=2.28e5;y[17]=3.05e5;

	//Error of cross section
	z[0]=0;z[1]=0;z[2]=0.13e2;z[3]=0.4e2;z[4]=0.25e2;z[5]=0.08e3;z[6]=0.14e3;z[7]=0.01e4;z[8]=0.06e4;z[9]=0.23e4;
	z[10]=0.48e4;z[11]=0.83e4;z[12]=1.27e4;z[13]=1.8e4;z[14]=0.24e5;z[15]=0.31e5;z[16]=0.47e5;z[17]=0.67e5;
	int n = 18;
	TFile* file = new TFile("./data/ClAbsorpCS.root","RECREATE");
	TGraph* cs = new TGraph(n,x,y);
	cs->SetTitle("Chlorine Absorption Cross Section");
	cs->GetXaxis()->SetTitle("Energy/[MeV]");
	cs->GetYaxis()->SetTitle("Cross Section /[1e-46 cm^2]");
	cs->SetName("CrossSection");

	TGraph* gr = new TGraph(n,x,z);
	gr->SetTitle("Chlorine Absorption Cross Section Error");
	gr->GetXaxis()->SetTitle("Energy/[MeV]");
	gr->GetYaxis()->SetTitle("Cross Section /[1e-46 cm^2]");
	gr->SetName("Error");

	cs->Write();
	gr->Write();

	gROOT->ProcessLine(".q");
}
