#include "SurvProb.h"
#include <iostream>
#include "TMath.h"

using namespace Eigen;

SurvProb::SurvProb(int generation)
{
	gene = generation;
	DN[0] = "day";
	DN[1] = "night";
	osci = new Oscillation(gene);
	SunEarthDistance = 7.58122e11;
	DistPeriAp = 0.01672;
	RSun = 3.52461e9; //radius_sun/(hbar*c) = 6.955e8[m]/(hbar*c)
}

MatrixXcd SurvProb::CalculateRotateMatrix(double Energy, double radius, double L)
{
	std::complex<double> img(0,-1);
	double ne = solar.GetEDensity(radius);
	double nn = solar.GetNDensity(radius);
	double np = solar.GetPDensity(radius);
	MatrixXcd Ham = osci->Hamiltonian(Energy, ne, nn, np);
	return (img * RSun * Ham * L).exp();
}

MatrixXcd SurvProb::CalculateEarthRotate(double Energy, double ne, double nn, double np, double L)
{
	/* Since L outputs are of units km, the zero-unit Schrodinger's equation needs a coefficient. */
	double coeff = 5067.73094; /* from km to MeV/eV^2 */
	std::complex<double> img(0,-1);
	MatrixXcd Ham = osci->Hamiltonian(Energy, ne, nn, np);
	return (img * coeff * Ham * L).exp();
}

MatrixXcd SurvProb::CalculateVacuumRotate(double Energy, double L)
{
	/* L has been normalized. */
	std::complex<double> img(0,-1);
	MatrixXcd Ham = osci->Hamiltonian(Energy, 0, 0);
	return (img * SunEarthDistance * Ham * L).exp();
}

void SurvProb::SetupRotateMatrix(double Energy, int &nsteps, double* radius, MatrixXcd* Rotate)
{
	/* First step is not important since step size is adjusted automatically. */
	/* The only thing is, if first step is too large, numeric calculation will get unprecise when level crossing happens at solar core coincidently. */

	double firststep = 1e-5; 
	nsteps = 0;
	double Radius = 0;
	double tm12_prev, tm13_prev, tm14_prev;
	tm12_prev = tm13_prev = tm14_prev = 0; // Theta at previous step. Save outside loop to compare with latter step.
	radius[0] = firststep;
	while(Radius<=1)
	{
		double ne = solar.GetEDensity(Radius);
		double nn = solar.GetNDensity(Radius);
		double np = solar.GetPDensity(Radius);
		Rotate[nsteps] = CalculateRotateMatrix(Energy, Radius, radius[nsteps]);
		Radius += radius[nsteps];
		nsteps++;

		/* Change step automatically */
		double HighBound = 2e-4;
		double LowBound = 1e-4;
		if(gene==4)
		{
			HighBound = 6e-4;
			LowBound = 3e-4;
		}
		if(gene==2)
		{
			double tm12 = osci->GetTheta12(Energy, ne, nn);
			if(abs(tm12 - tm12_prev) > HighBound) radius[nsteps] = radius[nsteps-1] * 0.8;
			else if(abs(tm12 - tm12_prev) < LowBound)  radius[nsteps] = radius[nsteps-1] * 1.2; 
			else radius[nsteps] = radius[nsteps-1];
			tm12_prev = tm12;
		}
		if(gene==3)
		{
			double tm12 = osci->GetTheta12(Energy, ne, nn, np);
			double tm13 = osci->GetTheta13(Energy, ne, nn, np);
			if(abs(tm12 - tm12_prev) > HighBound || abs(tm13 - tm13_prev) > HighBound) radius[nsteps] = radius[nsteps-1] * 0.8;
			else if(abs(tm12 - tm12_prev) < LowBound  && abs(tm13 - tm13_prev) < LowBound)  radius[nsteps] = radius[nsteps-1] * 1.2;
			else radius[nsteps] = radius[nsteps-1];
			tm12_prev = tm12;
			tm13_prev = tm13;
		}
		if(gene==4)
		{
			double tm12 = osci->GetTheta12(Energy, ne, nn);
			double tm13 = osci->GetTheta13(Energy, ne, nn);
			double tm14 = osci->GetTheta14(Energy, ne, nn);
			if(abs(tm12 - tm12_prev) > HighBound || abs(tm13 - tm13_prev) > HighBound || abs(tm14 - tm14_prev) > HighBound) radius[nsteps] = radius[nsteps-1] * 0.8;
			else if(abs(tm12 - tm12_prev) < LowBound  && abs(tm13 - tm13_prev) < LowBound  && abs(tm14 - tm14_prev) < LowBound)  radius[nsteps] = radius[nsteps-1] * 1.2;
			else radius[nsteps] = radius[nsteps-1];
			tm12_prev = tm12;
			tm13_prev = tm13;
			tm14_prev = tm14;
		}
	}
	
	/* Deal with last rotation. */
	Rotate[nsteps-1] = CalculateRotateMatrix(Energy, Radius-radius[nsteps-1], 1-(Radius-radius[nsteps-1]));
	radius[nsteps-1] = 1 - (Radius - radius[nsteps-1]);
}

void SurvProb::SetupProb(double Energy, double eeel, double eeer, double etel, double eter, double eeu, double eed, double euu, double eud, double etu, double etd, double uuu, double uud, double utu, double utd, double ttu, double ttd)
{
	//osci->SetupParameters(sin_2_t12, sin_2_t13, 0.5, sin_2_t14, ms12, ms13, ms14, 0); //No CP violation and maximal mixing between 2 and 3.
	osci->SetupNSI(eeel, eeer, etel, eter, eeu, eed, euu, eud, etu, etd, uuu, uud, utu, utd, ttu, ttd);

	/* Setup mass eigenstates in flavor basis. */
	VectorXcd Mass[gene];
	MatrixXcd Uv = osci->GetVacuumMixingMatrix();
	for(int i=0; i<gene; i++)
	{
		Mass[i] = VectorXcd::Zero(gene);
		for(int j=0; j<gene; j++)
		{
			Mass[i](j) = Uv.transpose().conjugate()(i,j);
		}
	}

	/* Get flux distribution for GS model. */
	for(int comp=1; comp<=9; comp++)
	{
		Dist[comp-1] = solar.GetBinnedFluxDist(model, comp);
	}

	/* Initiate state vector. */
	int nDistBins = Dist[0]->GetNbinsX(); // All components share the same binning.
	double DistRadius[nDistBins]; // Center position for distribution bins.
	VectorXcd Psi[nDistBins]; // State vector for each distribution bin.
	VectorXcd Flavor[gene]; // Flavor eigen vector in flavor basis.
	for(int i=0; i<gene; i++)
	{
		Flavor[i] = VectorXcd::Zero(gene);
		for(int j=0; j<gene; j++) 
		{
			Flavor[i](j) = 0;
			Flavor[i](i) = 1;
		}
	}
	
	
	for(int nbins=0; nbins<nDistBins; nbins++)
	{
		Psi[nbins] = Flavor[0];
		DistRadius[nbins] = Dist[0]->GetBinCenter(nbins+1);
	}

	/* Determine whether solar numeric calculation should be done, by comparing L with Sun Radius. */
	/* Since 3 is degenerated from 1 and 2 in the Sun, only ms12 contributes. */
	bool IsSolar = true;
	double ms12 = 7.54e-5;
	double OsciLength12 = 4 * TMath::Pi() * Energy / ms12;
	double neMax = solar.GetEDensity(0);
	double nnMax = solar.GetNDensity(0);
	double npMax = solar.GetPDensity(0);
	double MassSplittingMax = osci->GetMass12(Energy, neMax, nnMax, npMax);
	double OsciLengthSun12 = 4 * TMath::Pi() * Energy / MassSplittingMax; 
	if((OsciLengthSun12 > RSun * 50) && (OsciLength12 > RSun * 50)) IsSolar = false;

	
	if(IsSolar)
	{
		/* Setup rotate matrix for each small slab. */
		int nsteps;
		double radius[200000];
		MatrixXcd Rotate[200000];
	
		SetupRotateMatrix(Energy, nsteps, radius, Rotate);
	
		//std::cout << "Total steps for numeric calculation:  " << nsteps << std::endl;
	
		/* Calculate state vector after Sun MSW. */
		MatrixXcd RotTotal = MatrixXcd::Identity(gene,gene);
		double Radius = 1;
		int DistBinNum = nDistBins-1;
		for(int slab=nsteps-1; slab>0; slab--)
		{
			RotTotal = RotTotal * Rotate[slab];
			Radius -= radius[slab];
			if(slab>0)
			{
				if(DistRadius[DistBinNum] > Radius - radius[slab-1])
				{
					MatrixXcd FirstRot = CalculateRotateMatrix(Energy, DistRadius[DistBinNum], Radius-DistRadius[DistBinNum]); // The first rotate matrix for state vector at bin center.
					Psi[DistBinNum] = RotTotal * FirstRot * Psi[DistBinNum];
					DistBinNum--;
				}
			}
			if(DistBinNum<0) break; // All distribution bins complete calculation.
		}
	}
	else{;} // If no solar numeric calculation, then all state vectors should be the same as original ones.


	/* Sun to Earth. For regions other than VAC, this part can be regarded incoherent, and there is no need to consider changes in distance between Sun and Earth. */
	/* See Note.txt for detailed analysis of degeneracy. */
	bool IsVac = false;

	/* Initialize operator for VAC. Detect <psi|Operator|psi>. */
	MatrixXcd Operator[2][gene];
	for(int daynight=0; daynight<2; daynight++)
	{
		for(int flavor=0; flavor<gene; flavor++)
		{
			Operator[daynight][flavor] = MatrixXcd::Zero(gene, gene);
			Operator[0][flavor](flavor,flavor) = 1; 
		}
	}

	/* Setup livetime. */
	TH1D* jplivetime = earth.GetLivetime("Jinping");
	TH1D* snolivetime = earth.GetLivetime("SNO");
	TH1D* sklivetime = earth.GetLivetime("SK");
	int nseg;
	double L[300], edens[300], ndens[300], pdens[300]; //100 shells for Earth density distribution.
	double jpNightLivetimeRatio = 0;
	double snoNightLivetimeRatio = 0;
	double skNightLivetimeRatio = 0;
	int zenbins = jplivetime->GetNbinsX(); //Zenith angle distribution bins.

	std::cout << "Solar finish." << std::endl;
	/* Save abs(amplitudes) from psi to mass eigenstate, and mass eigenstate to flavor eigenstate. For degeneracy. */
	double psi2mass[nDistBins][gene], mass2flv[gene][gene], mass2Ham2flv[gene][zenbins][gene][3];
	for(int bin=0; bin<nDistBins; bin++)
	{
		for(int mass=0; mass<gene; mass++) 
		{
			psi2mass[bin][mass] = std::abs((Mass[mass].transpose().conjugate() * Psi[bin])(0,0));
		}
	}
	for(int mass=0; mass<gene; mass++) for(int flv=0; flv<gene; flv++)
	{
		mass2flv[mass][flv] = std::abs((Flavor[flv].transpose().conjugate() * Mass[mass])(0,0));
	}

	/* Calculate Earth MSW . */
	for(int bin=1; bin<=zenbins; bin++)
	{
		/* Initialize mass2Ham2flv. */
		for(int mass=0; mass<gene; mass++)	for(int flavor=0; flavor<gene; flavor++) for(int exp=0; exp<3; exp++)
		{
			mass2Ham2flv[mass][bin-1][flavor][exp] = 0;
		}

		double zenith = jplivetime->GetBinCenter(bin);
		if(zenith >= 0) continue; //Downwards.
		double jpdlivetime = jplivetime->GetBinContent(bin); //Differential livetime.
		jpNightLivetimeRatio += jpdlivetime;
		double snodlivetime = snolivetime->GetBinContent(bin);
		snoNightLivetimeRatio += snodlivetime;
		double skdlivetime = sklivetime->GetBinContent(bin);
		skNightLivetimeRatio += skdlivetime;

		earth.Intersect(zenith, nseg, L, edens, ndens);
		MatrixXcd RotateForGivenZenith = MatrixXcd::Identity(gene, gene);
		MatrixXcd RotateForGivenShell[nseg];

		for(int EarthSeg=0; EarthSeg<nseg; EarthSeg++)
		{
			pdens[EarthSeg] = ndens[EarthSeg];
			if(2*EarthSeg < nseg) RotateForGivenShell[EarthSeg] = CalculateEarthRotate(Energy, edens[EarthSeg], ndens[EarthSeg], pdens[EarthSeg], L[EarthSeg]);
			else RotateForGivenShell[EarthSeg] = RotateForGivenShell[nseg-EarthSeg-1]; // Shorten calculation time using symmetry.
			RotateForGivenZenith = RotateForGivenZenith * RotateForGivenShell[EarthSeg];
		}

		if(IsVac){;}
		else // Do degenerated calculation
		{
			for(int mass=0; mass<gene; mass++)	for(int flavor=0; flavor<gene; flavor++)
			{
				mass2Ham2flv[mass][bin-1][flavor][0] = std::sqrt(jpdlivetime) * std::abs((Flavor[flavor].transpose().conjugate() * RotateForGivenZenith * Mass[mass])(0,0));
				mass2Ham2flv[mass][bin-1][flavor][1] = std::sqrt(snodlivetime) * std::abs((Flavor[flavor].transpose().conjugate() * RotateForGivenZenith * Mass[mass])(0,0));
				mass2Ham2flv[mass][bin-1][flavor][2] = std::sqrt(skdlivetime) * std::abs((Flavor[flavor].transpose().conjugate() * RotateForGivenZenith * Mass[mass])(0,0));
			}
		}
	}

	if(IsVac){;}
	else 
	{
		/* Energy eigenstates from Sun to Earth are degenerated. */
		for(int comp=0; comp<9; comp++) for(int flavor=0; flavor<gene; flavor++)
		{
			prob[0][flavor][comp][0] = 0;
			prob[1][flavor][comp][0] = 0;
			prob[0][flavor][comp][1] = 0;
			prob[1][flavor][comp][1] = 0;
			prob[0][flavor][comp][2] = 0;
			prob[1][flavor][comp][2] = 0;
			for(int mass=0; mass<gene; mass++) for(int distbin=1; distbin<=nDistBins; distbin++)
			{
				double dflux = Dist[comp]->GetBinContent(distbin) * Dist[comp]->GetBinWidth(distbin);
				prob[0][flavor][comp][0] += std::pow(psi2mass[distbin-1][mass] * mass2flv[mass][flavor], 2) * dflux;
				for(int bin=1; bin<zenbins; bin++)
				{
					prob[1][flavor][comp][0] += std::pow(psi2mass[distbin-1][mass] * mass2Ham2flv[mass][bin][flavor][0], 2) * dflux;
					prob[1][flavor][comp][1] += std::pow(psi2mass[distbin-1][mass] * mass2Ham2flv[mass][bin][flavor][1], 2) * dflux;
					prob[1][flavor][comp][2] += std::pow(psi2mass[distbin-1][mass] * mass2Ham2flv[mass][bin][flavor][2], 2) * dflux;
				}
			}
			prob[1][flavor][comp][0] /= jpNightLivetimeRatio;
			prob[1][flavor][comp][1] /= snoNightLivetimeRatio;
			prob[1][flavor][comp][2] /= skNightLivetimeRatio;
			prob[0][flavor][comp][1] = prob[0][flavor][comp][0];
			prob[0][flavor][comp][2] = prob[0][flavor][comp][0];
		}
	}
	std::cout << "Finish preparing all survival probabilities. " << std::endl;
}

double SurvProb::GetProbFromCalculation(int daynight, int flavor, int comp, int exp)
{
	return prob[daynight][flavor-1][comp-1][exp];
}

