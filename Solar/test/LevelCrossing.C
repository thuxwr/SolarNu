/*
	 Draw level crossing scheme inside the Sun.

	 Weiran, May 21, 2018.

*/

#include "../SolarNu.h"
#include "../../Oscillation/Osci.h"
#include "Eigen/Dense"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLatex.h"

using namespace Eigen;

Oscillation osci(4); //3 generations
SolarNu solar;
ComplexEigenSolver<MatrixXcd> ces;

int main()
{
	int npoints = 1000;
	double sin_2_t12 = 0.307;
	double ms12 = 7.5e-5;
	double sin_2_t13 = 0.0241;
	double sin_2_t23 = 0.5;
	double sin_2_t14 = 0.0001;
	double ms13 = 2.5e-3;
	double ms14 = 1.2e-5;
	double energy = 1; //1MeV
	//double energy = 4;
	osci.SetupParameters(sin_2_t12, sin_2_t13, sin_2_t23, sin_2_t14, ms12, ms13, ms14);
	double radius[npoints], nu1[npoints], nu2[npoints], nu3[npoints];
	double levels = -1000;
	for(int i=0; i<npoints; i++) 
	{
		radius[i] = ((double)i)/npoints;
		double ne = solar.GetEDensity(radius[i]);
		double nn = solar.GetNDensity(radius[i]);
		MatrixXcd Ham = osci.Hamiltonian(energy, ne, nn);
		ces.compute(Ham);

		nu1[i] = ces.eigenvalues()[0].real() * 1e-6; //Change MeV to eV.
		nu2[i] = ces.eigenvalues()[1].real() * 1e-6; //Change MeV to eV.
		nu3[i] = ces.eigenvalues()[2].real() * 1e-6; //Change MeV to eV.

		/* Determine whether to switch nu1 and nu2. */
		if(energy==4 && nu2[i]>1.51e-12)
		{
			double tmp;
			tmp = nu1[i];
			nu1[i] = nu2[i];
			nu2[i] = tmp;
		}
	}

	TGraph* g1 = new TGraph(npoints, radius, nu1);
	TGraph* g2 = new TGraph(npoints, radius, nu2);
	TGraph* g3 = new TGraph(npoints, radius, nu3);

	double rangeYmax = 0;
	if(rangeYmax < g1->GetYaxis()->GetXmax()) rangeYmax = g1->GetYaxis()->GetXmax();
	if(rangeYmax < g3->GetYaxis()->GetXmax()) rangeYmax = g3->GetYaxis()->GetXmax();
	rangeYmax = 45e-12;

	TCanvas* c1 = new TCanvas;
	g1->GetYaxis()->SetRangeUser(0, rangeYmax);
	g1->SetLineColor(kRed);
	g1->SetTitle("");
	g1->GetXaxis()->SetTitle("R/R_{#odot}");
	g1->GetXaxis()->SetTitleSize(0.045);
	g1->GetYaxis()->SetTitle("Mass Eigenvalue [eV]");
	g1->GetYaxis()->SetTitleSize(0.045);
	g1->GetXaxis()->SetRangeUser(0,1);
	g3->SetLineColor(kRed);
	g1->SetLineWidth(2);
	g2->SetLineWidth(2);
	g3->SetLineWidth(2);
	g1->Draw("al");
	g2->Draw("lsame");
	g3->Draw("lsame");
	
	TLatex latex;
	latex.DrawLatex(0.6,40.5e-12,"Energy = 1MeV");

	c1->SaveAs("test.pdf");

	//TFile* file = new TFile("test.root", "RECREATE");
	//g1->SetName("g1");
	//g1->Write();
	//g2->SetName("g2");
	//g2->Write();
	//file->Close();

	return 0;
}

		

