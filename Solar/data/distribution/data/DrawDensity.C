/*
	 I found this data, which includes all radius inside the Sun.
	 Neutron number density is calculated using fractions of all elements and their neutron numbers in one atom.

	 Weiran, Feb 25, 2018
*/

#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TGraph.h"
#include "TMath.h"
#include "TAxis.h"

using namespace std;

int main()
{
	TFile* file = new TFile("../Density.root", "UPDATE");
	for(int times=1; times<2; times++)
	{
		if(times==0)
		{
			/* Draw electron number density.*/
			int nline = 1763;
			double R[nline], Dens[nline];
			ifstream fin("BP2000.dat");
			for(int i=0; i<nline; i++)
			{
				fin >> R[i] >> Dens[i];
			}

			TGraph* dens = new TGraph(nline, R, Dens);
			dens->SetName("EDensity");
			dens->GetXaxis()->SetTitle("radius/R_{sun}");
			dens->GetYaxis()->SetTitle("Log(n_{e}/N_{A})/cm^{3}");
			dens->Write();
		}
		else
		{
			/* Draw neutron/proton number density. */
			int nline = 1268;
			double smass[nline],radius[nline],temp[nline],rho[nline],pressure[nline],lum[nline],frac_H[nline],
						 frac_He4[nline],frac_He3[nline],frac_C[nline],frac_N[nline],frac_O[nline];
			double n_dens[nline];
			double p_dens[nline];
			ifstream fin("BS05_GS_elements.dat");
			for(int i=0; i<nline; i++)
			{
				fin >> smass[i] >> radius[i] >> temp[i] >> rho[i] >> pressure[i] >> lum[i] >> frac_H[i] >> frac_He4[i] >> frac_He3[i] >> frac_C[i] >> frac_N[i] >> frac_O[i];
				//n_dens[i] = rho[i]*(frac_He4[i]/4.002602*2 + frac_He3[i]/3.0160293 + frac_C[i]/12.*6 + frac_N[i]/14.003074*7 + frac_O[i]/15.9949146*8);
				//n_dens[i] = log(n_dens[i])/log(10);
				p_dens[i] = rho[i]*(frac_H[i]/1.008 + frac_He4[i]/4.002602*2 + frac_He3[i]/3.0160293*2 + frac_C[i]/12.*6 + frac_N[i]/14.003074*7 + frac_O[i]/15.9949146*8);
				p_dens[i] = log(p_dens[i])/log(10);
			}

			TGraph* dens = new TGraph(nline, radius, p_dens);
			dens->SetName("PDensity");
			dens->GetXaxis()->SetTitle("radius/R_{sun}");
			dens->GetYaxis()->SetTitle("Log(n_{p}/N_{A})/cm^{3}");
			dens->Write();
		}
	}
	file->Close();
	return 0;
}
			

