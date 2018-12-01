#include <string>
#include "TMath.h"
#include "SolarElec.h"

using namespace std;

SolarElec::SolarElec(int generation, int Option)
{
	gene = generation;
	option = Option;
	surv = new SurvProb(gene);
	for(int comp=1; comp<=9; comp++)
	{
		spec[comp-1] = solar.GetSpec(1, comp);
	}
	IsNSI = false;
	binwidth = 5e-3;
}

TH1D* SolarElec::GetElecSpec(double sin_2_theta, double ms, int daynight, int comp, double flux, string experiment)
{
	TGraph* SurvP = surv->GetProb(sin_2_theta, ms, comp, daynight, "e", experiment);
	spec[comp-1]->Scale(1/spec[comp-1]->Integral("width"));
	spec[comp-1]->Scale(flux);

	double NuEmin = spec[comp-1]->GetXaxis()->GetXmin();
	double NuEmax = spec[comp-1]->GetXaxis()->GetXmax();
	int NuNbin = spec[comp-1]->GetXaxis()->GetNbins();

	TH1D* Nue = new TH1D("","",NuNbin,NuEmin,NuEmax);
	TH1D* Numu = new TH1D("","",NuNbin,NuEmin,NuEmax);
	TH1D* Nutau = new TH1D("","",NuNbin,NuEmin,NuEmax);
	TH1D* Nuo = new TH1D("","",NuNbin,NuEmin,NuEmax);

	TGraph* ProbMu = surv->GetProb(sin_2_theta, ms, comp, daynight, "mu", experiment);
	TGraph* ProbTau = surv->GetProb(sin_2_theta, ms, comp, daynight, "tau", experiment);

	for(int bin=1; bin<=NuNbin; bin++)
	{
		double dFlux = spec[comp-1]->GetBinContent(bin);
		double Energy = spec[comp-1]->GetBinCenter(bin);
		double Sprob = SurvP->Eval(Energy);
		if(option==1) Sprob = 1;

		/* Convert output to unit: /cm^2 /s /MeV. */
		Nue->SetBinContent(bin, dFlux*Sprob*1e10);
		if(gene<=3) Nuo->SetBinContent(bin, dFlux*(1-Sprob)*1e10);
		else //Deal with sterile neutrinos.
		{
	    TGraph* ProbS = surv->GetProb(sin_2_theta, ms, comp, daynight, "s");
			double probS = ProbS->Eval(Energy);
			if(option==1) probS = 0;
			Nuo->SetBinContent(bin, dFlux*(1-Sprob-probS)*1e10);
			delete ProbS;
		}
		if(IsNSI)
		{
			double probmu = ProbMu->Eval(Energy);
			Numu->SetBinContent(bin, dFlux*probmu * 1e10);
			double probtau = ProbTau->Eval(Energy);
			Nutau->SetBinContent(bin, dFlux*probtau * 1e10);
		}
	}
	delete ProbMu; delete ProbTau;

	double TeEmax = NuEmax;
	int TeNbin = int(floor(TeEmax/binwidth)) +1;

	/* Convert neutrino spectra to electron energy spectra. */
	TH1D* Tee = new TH1D("","",TeNbin, 0, TeNbin * binwidth);
	TH1D* Teo = new TH1D("","",TeNbin, 0, TeNbin * binwidth);
	scat.VisElecSpec(Nue, Tee, 0, 0);
	if(IsNSI)
	{
		TH1D* Temu = new TH1D("","",TeNbin, 0, TeNbin * binwidth);
		TH1D* Tetau = new TH1D("","",TeNbin, 0, TeNbin * binwidth);
		scat.VisElecSpec(Numu, Temu, 1, 0);
		scat.VisElecSpec(Nutau, Tetau, 2, 0);
		Teo->Add(Temu, Tetau);
		delete Temu; delete Tetau;
	}
	else scat.VisElecSpec(Nuo, Teo, 1, 0);

	TH1D* ElecSpec = new TH1D(*Tee);
	ElecSpec->Reset();
	ElecSpec->Add(Tee, Teo);

	ElecSpec->SetName("ElecSpec");
	ElecSpec->GetXaxis()->SetTitle("Electron kinetic energy [MeV]");
	ElecSpec->GetYaxis()->SetTitle("Event rate [(s MeV)^{-1}]");

	delete Nue; delete Nuo; delete Tee; delete Teo;
	delete Numu; delete Nutau;
	return ElecSpec;
}

TH1D* SolarElec::GetElecSpec(int daynight, int model, int comp)
{
	/* Get flux according to model. */
	double flux = solar.GetFlux(model, comp);

	double sin_2_theta = 0.327, ms = 4.8e-5;
	return GetElecSpec(sin_2_theta, ms, daynight, comp, flux);
}

TH1D* SolarElec::GetElecSpec(int model, int comp)
{
	TH1D* dayspec = GetElecSpec(0, model, comp);
	TH1D* nightspec = GetElecSpec(1, model, comp);
	TH1D* totalspec = (TH1D*)dayspec->Clone();
	totalspec->Add(dayspec, nightspec);
	totalspec->Scale(0.5);
	delete dayspec; delete nightspec;
	return totalspec;
}

TH1D* SolarElec::GetElecTotal(double sin_2_theta, double ms, int daynight, double* flux)
{
	TH1D* ElecSpec[9];
	double Emax = 0;
	for(int comp=1; comp<=9; comp++)
	{
		ElecSpec[comp-1] = GetElecSpec(sin_2_theta, ms, daynight, comp, flux[comp-1]);
		double emax = ElecSpec[comp-1]->GetXaxis()->GetXmax();
		if(Emax<=emax) Emax = emax;
	}
	int nbins = int(floor(Emax/binwidth))+1;
	TH1D* ElecTotal = new TH1D("","",nbins,0,nbins*binwidth);
	
	for(int comp=1; comp<=9; comp++)
	{
		for(int bin=1; bin<=ElecSpec[comp-1]->GetNbinsX(); bin++)
		{
			double add = ElecSpec[comp-1]->GetBinContent(bin);
			double cur = ElecTotal->GetBinContent(bin);
			ElecTotal->SetBinContent(bin, cur+add);
		}
	}

	ElecTotal->GetXaxis()->SetTitle("Electron kinetic energy [MeV]");
	ElecTotal->GetYaxis()->SetTitle("Event rate [(s MeV electron)^{-1}]");

	for(int i=0; i<9; i++) delete ElecSpec[i];
	return ElecTotal;
}

TH1D* SolarElec::GetElecTotal(int daynight, int model, double sin_2_theta, double ms)
{
	double flux[9];
	for(int comp=1; comp<=9; comp++) flux[comp-1] = solar.GetFlux(model, comp);
	return GetElecTotal(sin_2_theta, ms, daynight, flux);
}

TH1D* SolarElec::GetElecSpec(double sin_2_theta, double ms, int comp, double flux, string experiment)
{
	TH1D* dayspec = GetElecSpec(sin_2_theta, ms, 0, comp, flux, experiment);
	TH1D* nightspec = GetElecSpec(sin_2_theta, ms, 1, comp, flux, experiment);
	TH1D* totalspec = new TH1D(*dayspec);
	for(int bin=1; bin<=totalspec->GetNbinsX(); bin++)
	{
		double daycontent = dayspec->GetBinContent(bin);
		double nightcontent = nightspec->GetBinContent(bin);
		double content = (daycontent+nightcontent)/2.;
		totalspec->SetBinContent(bin, content);
	}
	delete dayspec; delete nightspec;
	return totalspec;
}

TH1D* SolarElec::GetElecSpec(TF1* SurvProb, int comp, double flux)
{
	spec[comp-1]->Scale(1/spec[comp-1]->Integral("width"));
	spec[comp-1]->Scale(flux);

	double NuEmin = spec[comp-1]->GetXaxis()->GetXmin();
	double NuEmax = spec[comp-1]->GetXaxis()->GetXmax();
	int NuNbin = spec[comp-1]->GetXaxis()->GetNbins();

	TH1D* Nue = new TH1D("","",NuNbin,NuEmin,NuEmax);
	TH1D* Nuo = new TH1D("","",NuNbin,NuEmin,NuEmax);

	for(int bin=1; bin<=NuNbin; bin++)
	{
		double dFlux = spec[comp-1]->GetBinContent(bin);
		double Energy = spec[comp-1]->GetBinCenter(bin);
		double Sprob = SurvProb->Eval(Energy);
		if(option==1) Sprob = 1;

		/* Convert output to unit: /cm^2 /s /MeV. */
		Nue->SetBinContent(bin, dFlux*Sprob*1e10);
		Nuo->SetBinContent(bin, dFlux*(1-Sprob)*1e10);
	}

	double TeEmax = NuEmax;
	int TeNbin = int(floor(TeEmax/binwidth)) +1;

	/* Convert neutrino spectra to electron energy spectra. */
	TH1D* Tee = new TH1D("","",TeNbin, 0, TeNbin * binwidth);
	TH1D* Teo = new TH1D("","",TeNbin, 0, TeNbin * binwidth);
	scat.VisElecSpec(Nue, Tee, 0, 0);
	scat.VisElecSpec(Nuo, Teo, 1, 0);

	TH1D* ElecSpec = new TH1D(*Tee);
	ElecSpec->Reset();
	ElecSpec->Add(Tee, Teo);

	ElecSpec->SetName("ElecSpec");
	ElecSpec->GetXaxis()->SetTitle("Electron kinetic energy [MeV]");
	ElecSpec->GetYaxis()->SetTitle("Event rate [(s MeV)^{-1}]");

	delete Nue, Nuo, Tee, Teo;
	return ElecSpec;
}
	
TH1D* SolarElec::GetElecSpec(TGraph* SurvProb, int comp, double flux)
{
	spec[comp-1]->Scale(1/spec[comp-1]->Integral("width"));
	spec[comp-1]->Scale(flux);

	double NuEmin = spec[comp-1]->GetXaxis()->GetXmin();
	double NuEmax = spec[comp-1]->GetXaxis()->GetXmax();
	int NuNbin = spec[comp-1]->GetXaxis()->GetNbins();

	TH1D* Nue = new TH1D("","",NuNbin,NuEmin,NuEmax);
	TH1D* Nuo = new TH1D("","",NuNbin,NuEmin,NuEmax);

	for(int bin=1; bin<=NuNbin; bin++)
	{
		double dFlux = spec[comp-1]->GetBinContent(bin);
		double Energy = spec[comp-1]->GetBinCenter(bin);
		double Sprob = SurvProb->Eval(Energy);
		if(option==1) Sprob = 1;

		/* Convert output to unit: /cm^2 /s /MeV. */
		Nue->SetBinContent(bin, dFlux*Sprob*1e10);
		Nuo->SetBinContent(bin, dFlux*(1-Sprob)*1e10);
	}

	double TeEmax = NuEmax;
	int TeNbin = int(floor(TeEmax/binwidth)) +1;

	/* Convert neutrino spectra to electron energy spectra. */
	TH1D* Tee = new TH1D("","",TeNbin, 0, TeNbin * binwidth);
	TH1D* Teo = new TH1D("","",TeNbin, 0, TeNbin * binwidth);
	scat.VisElecSpec(Nue, Tee, 0, 0);
	scat.VisElecSpec(Nuo, Teo, 1, 0);

	TH1D* ElecSpec = new TH1D(*Tee);
	ElecSpec->Reset();
	ElecSpec->Add(Tee, Teo);

	ElecSpec->SetName("ElecSpec");
	ElecSpec->GetXaxis()->SetTitle("Electron kinetic energy [MeV]");
	ElecSpec->GetYaxis()->SetTitle("Event rate [(s MeV)^{-1}]");

	delete Nue, Nuo, Tee, Teo;
	return ElecSpec;
}

void SolarElec::SetupNSI(double el, double er, double tl, double tr)
{
	if(gene!=3)
	{
		cout << "NSI analysis should only be used for 3nu paradigm. " << endl;
		exit(0);
	}
	scat.SetupNSI(el, er, tl, tr);
	IsNSI = true;
}
