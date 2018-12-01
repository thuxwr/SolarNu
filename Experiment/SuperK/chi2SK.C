#include "chi2SK.h"
#include "TMath.h"
#include <iostream>
#include "Eigen/Dense"

using namespace std;

chi2SK* chi2SK::mSK;

chi2SK::chi2SK(int generation)
{
	mSK = this;
	gene = generation;
	for(int phase=1; phase<=4; phase++)
	{
		SKspec[phase-1] = SuperK.GetSpec(phase);
		SKspec[phase-1]->SetName(ExpName(phase).c_str());
	}
	ElecOsci = new SolarElec(gene);
	ElecNoOsci = new SolarElec(gene, 1);
	binwidth = 5e-3;
	nbins = int(floor(20./binwidth))+1;

	B8fluxMC = 5.25e6;
	hepfluxMC = 7.88e3;

	/* Setup unoscillated spectrum. */
	for(int phase=1; phase<=4; phase++)
	{
		TH1D* UnOsciB8 = GetDetSpec(0.3, 7.5e-5, 6, B8fluxMC, phase, 1);
		TH1D* UnOscihep = GetDetSpec(0.3, 7.5e-5, 3, hepfluxMC, phase, 1);
		UnOscispec[phase-1] = new TH1D("","",SuperK.GetNbinsX(phase), SuperK.GetXbins(phase));
		UnOscispec[phase-1]->Add(UnOsciB8, UnOscihep);
		delete UnOsciB8; delete UnOscihep;

		for(int daynight=0; daynight<3; daynight++)
		{
			B8spec[phase-1][daynight] = new TH1D();
			hepspec[phase-1][daynight] = new TH1D();
		}
	}

	minuit = new TMinuit(9); //4*2 for deltaS and deltaR, and 1 for correlated deltaB.
	minuit->SetPrintLevel(-1);

	string nu = getenv("neutrino");
	if(nu=="")
	{
		cout << "Environment variable 'neutrino' undefined." << endl;
		exit(0);
	}

	//file = new TFile((nu+"/Earth/test/result/SKEarthAffectRegion.root").c_str(), "READ");
	//earth = (TH2D*)file->Get("earth");

	/* Calculate uncertainty for measuring total flux. */
	for(int phase=1; phase<=4; phase++)
	{
		double Sigma0_2_Inverse = 0;
		double average = 0; //Average detected/MC rate.
		for(int bin=1; bin<=SuperK.GetNbinsX(phase); bin++)
		{
			Sigma0_2_Inverse += 1./pow(SKspec[phase-1]->GetBinError(bin), 2);
			average += SKspec[phase-1]->GetBinContent(bin)/pow(SKspec[phase-1]->GetBinError(bin), 2);
		}
		average /= Sigma0_2_Inverse;

		sigma0_2[phase-1] = 1./Sigma0_2_Inverse;
		sigmar_2[phase-1] = pow(average * SuperK.GetFluxSysError(phase), 2);
		alpha[phase-1] = sigma0_2[phase-1]/(sigma0_2[phase-1] + sigmar_2[phase-1]);
	}

	//file->Close();
	//delete file;
}

chi2SK::~chi2SK()
{
	delete ElecOsci;
	delete ElecNoOsci;
	for(int phase=1; phase<=4; phase++)
	{
		delete SKspec[phase-1];
		delete UnOscispec[phase-1];
	}
}

void chi2SK::SetupParameter(double sin_2_theta, double ms)
{
	/* Pass parameters to other fcns. */
	msin_2_theta = sin_2_theta;
	mms = ms;

	for(int phase=1; phase<=4; phase++)
	{
		for(int daynight=0; daynight<3; daynight++)
		{
			B8spec[phase-1][daynight]->Reset();
			B8spec[phase-1][daynight]->SetBins(SuperK.GetNbinsX(phase), SuperK.GetXbins(phase));
			hepspec[phase-1][daynight]->Reset();
			hepspec[phase-1][daynight]->SetBins(SuperK.GetNbinsX(phase), SuperK.GetXbins(phase));

			TH1D* PredB8 = GetPredB8(sin_2_theta, ms, phase, daynight);
			TH1D* Predhep = GetPredhep(sin_2_theta, ms, phase, daynight);

			for(int bin=1; bin<=B8spec[phase-1][daynight]->GetNbinsX(); bin++)
			{
				B8spec[phase-1][daynight]->SetBinContent(bin, PredB8->GetBinContent(bin)/UnOscispec[phase-1]->GetBinContent(bin));
				hepspec[phase-1][daynight]->SetBinContent(bin, Predhep->GetBinContent(bin)/UnOscispec[phase-1]->GetBinContent(bin));
			}

			delete PredB8; delete Predhep;
		}
	}
}

double chi2SK::chi2spec(double B8flux, double hepflux)
{
	/* Set parameters. */
	pB8flux = B8flux;
	phepflux = hepflux;

	/* delta_x is a gaussian distribution with center value 0 and sigma 1. */
	int ierflg;
	minuit->mnparm(0,"deltaB",0,1,0,0,ierflg);
	minuit->mnparm(1,"deltaS_phase1",0,1,0,0,ierflg);
	minuit->mnparm(2,"deltaR_phase1",0,1,0,0,ierflg);
	minuit->mnparm(3,"deltaS_phase2",0,1,0,0,ierflg);
	minuit->mnparm(4,"deltaR_phase2",0,1,0,0,ierflg);
	minuit->mnparm(5,"deltaS_phase3",0,1,0,0,ierflg);
	minuit->mnparm(6,"deltaR_phase3",0,1,0,0,ierflg);
	minuit->mnparm(7,"deltaS_phase4",0,1,0,0,ierflg);
	minuit->mnparm(8,"deltaR_phase4",0,1,0,0,ierflg);

	minuit->SetFCN(skfcn);
	minuit->SetErrorDef(1);

	double arglist[10]; //Argument list.
	arglist[0] = 1e6;

	minuit->mnexcm("MIGRAD",arglist,0,ierflg);

	/* Output chi2. */
	double chi2 = 0;
	double edm, errdef;
	int nvpar, nparx, icstat;
	minuit->mnstat(chi2, edm, errdef, nvpar, nparx, icstat);
	mchi2 = chi2;

	/* Save parameters from spectrum fit and pass to time variation calculation. */
	double err;
	minuit->GetParameter(0, mdeltaB, err);
	minuit->GetParameter(1, mdeltaS[0], err);
	minuit->GetParameter(2, mdeltaR[0], err);
	minuit->GetParameter(3, mdeltaS[1], err);
	minuit->GetParameter(4, mdeltaR[1], err);
	minuit->GetParameter(5, mdeltaS[2], err);
	minuit->GetParameter(6, mdeltaR[2], err);
	minuit->GetParameter(7, mdeltaS[3], err);
	minuit->GetParameter(8, mdeltaR[3], err);

	/* Calculate fluxes for each phase while minimum. */
	double Par[9]; 
	Par[0] = mdeltaB; Par[1] = mdeltaS[0]; Par[2] = mdeltaR[0]; Par[3] = mdeltaS[1]; Par[4] = mdeltaR[1];
	Par[5] = mdeltaS[2]; Par[6] = mdeltaR[2]; Par[7] = mdeltaS[3]; Par[8] = mdeltaR[3];

	for(int phase=1; phase<=4; phase++)
	{
		/* Determine chi2_p_min. */
		TMinuit* Min = new TMinuit(2); //For flux minimization
		Min->SetPrintLevel(-1);
		int ierflg;
		mphase = phase;
		mpar = Par;
		Min->mnparm(0, "B8flux", 5.25e6, 0.2e6, 0, 0, ierflg);
		Min->mnparm(1, "hepflux", 8e3, 16e3, 0, 0, ierflg);
		Min->SetFCN(phasefcn);

		double arglist[10];
		arglist[0] = 1e5;
		Min->mnexcm("MIGRAD", arglist, 0, ierflg);

		double err1, err2;
		Min->GetParameter(0, B8fluxmin[phase-1], err1);
		Min->GetParameter(1, hepfluxmin[phase-1], err2);

		delete Min;
	}

	return chi2;
}
			
//double chi2SK::chi2tv()
//{
//	double xearth = msin_2_theta/(1-msin_2_theta);
//	double yearth = mms;
//	double ADN;
//	if(yearth>-9) ADN = earth->Interpolate(xearth, yearth);
//	else ADN = 0;
//	/* Best-fit day night asymmetry from SK is 3.3 \pm 1.118 percent. */
//	return pow((3.3*0.01-ADN)/(1.118*0.01), 2);
//}

double chi2SK::chi2tv()
{
	double chi2 = 0;
	for(int phase=1; phase<=4; phase++)
	{
		double pchi2, pchi2min;
		pchi2 = GetPhasetvChi2(phase, pB8flux, phepflux);
		pchi2min = GetPhasetvChi2(phase, B8fluxmin[phase-1], hepfluxmin[phase-1]);

		/* Add systematic uncertainty on total flux. */
		pchi2 = pchi2min + alpha[phase-1] * (pchi2 - pchi2min);

		pchi2 += pow(mdeltaS[phase-1], 2);
		pchi2 += pow(mdeltaR[phase-1], 2);

		chi2 += pchi2;
	}

	chi2 += pow(mdeltaB, 2);
	chi2 -= mchi2;

	return chi2;
}

double chi2SK::GetPhasetvChi2(int phase, double B8flux, double hepflux)
{
	double chi2 = 0;

	/* Get SK day spec and night spec. */
	TH1D* day = SuperK.GetSpec(phase, 0);
	TH1D* night = SuperK.GetSpec(phase, 1);

	for(int bin=1; bin<=SuperK.GetNbinsX(phase); bin++)
	{
		double energy = SKspec[phase-1]->GetBinCenter(bin);

		/* Construct cov matrix. For each energy bin, chi2_daynight is defined by cov(day, night). */
		Eigen::Matrix2d cov;
		cov(0,0) = pow(day->GetBinError(bin), 2);
		cov(1,1) = pow(night->GetBinError(bin), 2);

		double s1sys = SuperK.GetUncorreSysError(energy, phase) * 0.01 * day->GetBinContent(bin);
		double s2sys = SuperK.GetUncorreSysError(energy, phase) * 0.01 * night->GetBinContent(bin);
		cov(0,1) = cov(1,0) = s1sys*s2sys;

		/* Calculate chi2. */
		Eigen::Vector2d v;
		double shape_factor = 1;

		/* Distortion due to deltaB. */
		if(mdeltaB>=0) shape_factor /= (1 + mdeltaB * SuperK.GetCorreSysError(energy, 1, phase, 0));
		else shape_factor /= (1 - mdeltaB * SuperK.GetCorreSysError(energy, 1, phase, 1));

		/* Distortion due to deltaS. */
		if(mdeltaS[phase-1]>=0) shape_factor /= (1 + mdeltaS[phase-1] * SuperK.GetCorreSysError(energy, 2, phase, 0));
		else shape_factor /= (1 - mdeltaS[phase-1] * SuperK.GetCorreSysError(energy, 2, phase, 1));

		/* Distortion due to deltaR. */
		if(mdeltaR[phase-1]>=0) shape_factor /= (1 + mdeltaR[phase-1] * SuperK.GetCorreSysError(energy, 3, phase, 0));
		else shape_factor /= (1 - mdeltaR[phase-1] * SuperK.GetCorreSysError(energy, 3, phase, 1));

		double betab = B8flux / B8fluxMC * B8spec[phase-1][0]->GetBinContent(bin);
		double etah = hepflux / hepfluxMC * hepspec[phase-1][0]->GetBinContent(bin);
		v(0) = day->GetBinContent(bin) - (betab+etah) * shape_factor; //Day, Det-Pred

		betab = B8flux / B8fluxMC * B8spec[phase-1][1]->GetBinContent(bin);
		etah = hepflux / hepfluxMC * hepspec[phase-1][1]->GetBinContent(bin);
		v(1) = night->GetBinContent(bin) - (betab+etah) * shape_factor; //Night, Det-Pred

		/* The first bin in SK-II for day-night asymmetry is empty. */
		if(phase==2 && bin==1)
		{
			double sigma = SKspec[phase-1]->GetBinError(bin);
			double d = SKspec[phase-1]->GetBinContent(bin);
			betab = B8flux / B8fluxMC * B8spec[phase-1][2]->GetBinContent(bin);
			etah = hepflux / hepfluxMC * hepspec[phase-1][2]->GetBinContent(bin);

			chi2 += pow((d - (betab + etah) * shape_factor) / sigma, 2);
			continue;
		}

		chi2 += v.transpose() * cov.inverse() * v;
	}

	delete day; delete night;
	return chi2;
}

string chi2SK::ExpName(int phase)
{
	string name;
	if(phase==1) name = "SKI";
	else if(phase==2) name = "SKII";
	else if(phase==3) name = "SKIII";
	else if(phase==4) name = "SKIV";
	else name = "";
	return name;
}

TH1D* chi2SK::GetDetSpec(double sin_2_theta, double ms, int comp, double flux, int phase, int IsOsci, int daynight)
{
	TH1D* ElecSpec;
	if(daynight==2) //Average of day and night det spec.
	{
		if(IsOsci==0) ElecSpec = ElecOsci->GetElecSpec(sin_2_theta, ms, comp, flux, "SK");
		else ElecSpec = ElecNoOsci->GetElecSpec(sin_2_theta, ms, comp, flux, "SK");
	}
	else
	{
		if(IsOsci==0) ElecSpec = ElecOsci->GetElecSpec(sin_2_theta, ms, daynight, comp, flux, "SK");
		else ElecSpec = ElecNoOsci->GetElecSpec(sin_2_theta, ms, daynight, comp, flux, "SK");
	}

	TH1D* DetSpecTmp;
	res.Convert(DetSpecTmp, ElecSpec, ExpName(phase));

	/* Rebin according to phase. */
	TH1D* DetSpec = new TH1D("","",SuperK.GetNbinsX(phase), SuperK.GetXbins(phase));
	for(int bin=1; bin<=DetSpec->GetNbinsX(); bin++)
	{
		double InteLow = DetSpec->GetBinLowEdge(bin);
		double InteUp = DetSpec->GetBinWidth(bin) + InteLow;
		int LowBin = DetSpecTmp->FindBin(InteLow);
		int UpBin = DetSpecTmp->FindBin(InteUp)-1;

		double DetRate = DetSpecTmp->Integral(LowBin, UpBin);
		DetSpec->SetBinContent(bin, DetRate);
	}
	delete ElecSpec; delete DetSpecTmp;
	return DetSpec;
}

TH1D* chi2SK::GetPredB8(double sin_2_theta, double ms, int phase, int daynight)
{
	return GetDetSpec(sin_2_theta, ms, 6, B8fluxMC, phase, 0, daynight);
}

TH1D* chi2SK::GetPredhep(double sin_2_theta, double ms, int phase, int daynight)
{
	return GetDetSpec(sin_2_theta, ms, 3, hepfluxMC, phase, 0, daynight);
}

void chi2SK::skfcn(int &npar, double* gin, double &f, double* par, int iflag)
{
	double chi2 = 0;

	/* See ./test/Minimum.dat for more info. */
	for(int phase=1; phase<=4; phase++)
	{
		/* Determine chi2_p_min. */
		TMinuit* Min = new TMinuit(2); //For flux minimization
		Min->SetPrintLevel(-1);
		int ierflg;
		mSK->mphase = phase;
		mSK->mpar = par;
		Min->mnparm(0, "B8flux", 5.25e6, 0.2e6, 0, 0, ierflg);
		Min->mnparm(1, "hepflux", 8e3, 16e3, 0, 0, ierflg);
		Min->SetFCN(phasefcn);

		double arglist[10];
		arglist[0] = 1e5;
		Min->mnexcm("MIGRAD", arglist, 0, ierflg);

		double chi2pmin, pedm, perrdef;
		int pnvpar, pnparx, picstat;
		Min->mnstat(chi2pmin, pedm, perrdef, pnvpar, pnparx, picstat);

		/* Calculate chi2 for each phase. */
		double chi2p = 0;

		for(int bin=1; bin<= (mSK->SuperK).GetNbinsX(phase); bin++)
		{
			chi2p += mSK->GetBinChi2(phase, bin, par, mSK->pB8flux, mSK->phepflux);
		}

		/* Revise chi2p according to alpha_p. */
		chi2p = chi2pmin + mSK->alpha[phase-1] * (chi2p-chi2pmin);


		/* Add chi2 from nuisance parameters except for B8 spec. */
		chi2p += pow(par[phase*2-1], 2);
		chi2p += pow(par[phase*2], 2);

		chi2 += chi2p;
		delete Min;
	}

	/* Add chi2 from B8 spec theoretical uncertainty. */
	chi2 += pow(par[0], 2);

	f = chi2;
}

void chi2SK::phasefcn(int &npar, double* gin, double &f, double* par, int iflag)
{
	double chi2 = 0;
	for(int bin=1; bin<= (mSK->SuperK).GetNbinsX(mSK->mphase); bin++)
	{
		chi2 += mSK->GetBinChi2(mSK->mphase, bin, mSK->mpar, par[0], par[1]);
	}
	f = chi2;
}



double chi2SK::GetBinChi2(int phase, int bin, double* par, double B8flux, double hepflux)
{
	double chi2 = 0;

	double energy = SKspec[phase-1]->GetBinCenter(bin);
	double d = SKspec[phase-1]->GetBinContent(bin);
	double betab = B8flux / B8fluxMC * B8spec[phase-1][2]->GetBinContent(bin);
	double etah = hepflux / hepfluxMC * hepspec[phase-1][2]->GetBinContent(bin);
	double sigma = SKspec[phase-1]->GetBinError(bin);
	double shape_factor = 1;

	/* Distortion due to deltaB. */
	if(par[0]>=0) shape_factor /= (1 + par[0] * SuperK.GetCorreSysError(energy, 1, phase, 0));
	else shape_factor /= (1 - par[0] * SuperK.GetCorreSysError(energy, 1, phase, 1));

	/* Distortion due to deltaS. */
	if(par[phase*2-1]>=0) shape_factor /= (1 + par[phase*2-1] * SuperK.GetCorreSysError(energy, 2, phase, 0));
	else shape_factor /= (1 - par[phase*2-1] * SuperK.GetCorreSysError(energy, 2, phase, 1));

	/* Distortion due to deltaR. */
	if(par[phase*2]>=0) shape_factor /= (1 + par[phase*2] * SuperK.GetCorreSysError(energy, 3, phase, 0));
	else shape_factor /= (1 - par[phase*2] * SuperK.GetCorreSysError(energy, 3, phase, 1));

	/* Add chi2 from each bin. */
	chi2 = pow((d - (betab + etah) * shape_factor) / sigma, 2);
	return chi2;
}
