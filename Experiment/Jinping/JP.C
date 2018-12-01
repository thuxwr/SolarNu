#include "JP.h"
#include "TMath.h"
#include <sstream>

using namespace std;

JP* JP::mJP;

JP::JP(int generation, int model, int GlobalOrSolar, int option, int path)
{
	mJP = this;
	gene = generation;
	mmodel = model;
	elec = new SolarElec(gene);

	string nu = getenv("neutrino");
	if(nu=="")
	{
		cout << "Environment variable 'neutrino' undefined." << endl;
		exit(0);
	}

	string filepath = nu + "/Experiment/Jinping/Simulation/";
	string ExpName;

	if(option==1) ExpName = "Jinping200t";
	else ExpName = "Jinping2k";

	if(option==1) filepath = filepath + "200t/";
	if(option==2) filepath = filepath + "Bi210scan/";

	filepath = filepath + "result/";

	if(option==0) //golbal or solar
	{
		if(GlobalOrSolar==1) filepath = filepath+"global/";
		else filepath = filepath + "solar/";
	}
	else //for Bi210 scan.
	{
		filepath = filepath + "Bi210_1/";
	}

	if(path>=0)
	{
		stringstream ss;
		ss << path;
		string pathnum;
		ss >> pathnum;
		filepath = "/work/wangzhe9_work/neutrino/Experiment/Jinping/Sensitivity/Flux/result/tmp/" + pathnum + "/";
	}

	if(model==1) filepath = filepath + "GS98";
	else if(model==2) filepath = filepath + "AGS09";
	else
	{
		cout << "This model is not supported. Please choose model 1 or model 2." << endl;
		exit(0);
	}

	filepath = filepath + "-500PE.root";
	cout << "Open file from: " << filepath << endl;

	file = new TFile(filepath.c_str(), "READ");

	bkgname[0] = "Kr85"; bkgname[1] = "Bi210"; bkgname[2] = "C11"; bkgname[3] = "C14"; 
	bkgname[4] = "C10"; bkgname[5] = "Tl208"; bkgname[6] = "Be11"; bkgname[7] = "ExtTl208";

	for(int bkg=0; bkg<8; bkg++)
	{
		hbkg[bkg] = (TH1D*)file->Get((ExpName+bkgname[bkg]).c_str());
		bkgflux[bkg] = hbkg[bkg]->Integral();
		hbkg[bkg]->Scale(1/bkgflux[bkg]);
	}

	hsim[0] = (TH1D*)file->Get((ExpName+"Simday").c_str());
	hsim[1] = (TH1D*)file->Get((ExpName+"Simnight").c_str());

	hsimtot = (TH1D*)hsim[0]->Clone();
	hsimtot->Add(hsim[0], hsim[1]);

		/*            binwidth * day100ton         * Mass * Runtime */
	if(option==1) //200t
		ScaleFactor = 5e-3     * 24*60*60*3.307e31 * 2    * 750;
	else
		ScaleFactor = 5e-3     * 24*60*60*3.307e31 * 20   * 365/2. * 5;


	threshold = 0.05; //Default threshold: 50keV.

	EngUp = 20;	EngLow = 3.5; //For JP, this is reasonable.
	//EngUp = 15; EngLow = 3.5; //test JP at SK's threshold.
}

void JP::SetupParameter(double sin_2_theta, double ms)
{
	msin_2_theta = sin_2_theta;
	mms = ms;

	/* Get convoluted spectrum. */
	for(int comp=1; comp<=9; comp++) 
	{
		for(int daynight=0; daynight<2; daynight++)
		{
			delete hsig[daynight][comp-1];
			//hsig[daynight][comp-1]->Reset();
			TH1D* ElecSpec = elec->GetElecSpec(sin_2_theta, ms, daynight, comp, 1, "JP");
			TH1D* DetSpec;
			sim.Convert(DetSpec, ElecSpec, 2);
			DetSpec->Scale(ScaleFactor);
			hsig[daynight][comp-1] = DetSpec;
			delete ElecSpec;
		}
		hsigtot[comp-1] = (TH1D*)hsig[0][comp-1]->Clone();
		hsigtot[comp-1]->Add(hsig[0][comp-1], hsig[1][comp-1]);
		//hsigtot[comp-1]->Scale(0.5);
	}
}

void JP::SetupParameter(TF1* SurvProb)
{
	for(int comp=1; comp<=9; comp++) 
	{
		TH1D* ElecSpec = elec->GetElecSpec(SurvProb, comp, 2); 
		TH1D* DetSpec;
		sim.Convert(DetSpec, ElecSpec, 2);
		DetSpec->Scale(ScaleFactor);
		hsigtot[comp-1] = DetSpec;
	}
}

void JP::SetupParameter(TGraph* SurvProb)
{
	for(int comp=1; comp<=9; comp++) 
	{
		TH1D* ElecSpec = elec->GetElecSpec(SurvProb, comp, 2); 
		TH1D* DetSpec;
		sim.Convert(DetSpec, ElecSpec, 2);
		DetSpec->Scale(ScaleFactor);
		hsigtot[comp-1] = DetSpec;
	}
}

//double JP::chi2(double** bkgflux, double* nuflux, double penalty)
//{
//	/* chi2 = -2 log likelihood. */
//	double likelihood = 0;
//
//	/* All predicted and simulated spectra should share the same xmin(0) and same binwidth(5keV). */
//	for(int daynight=0; daynight<2; daynight++) for(int bin=1; bin<hsim[daynight]->GetNbinsX(); bin++)
//	{
//		double mea = hsim[daynight]->GetBinContent(bin);
//		mea /= (1 + penalty*0.015); //Systematic uncertainty: 1.5%
//		double pre = 0;
//
//		/* Add all background. */
//		for(int bkg=0; bkg<8; bkg++)
//		{
//			if(hbkg[bkg]->GetNbinsX()<bin) continue; 
//			pre += hbkg[bkg]->GetBinContent(bin) * bkgflux[daynight][bkg];
//		}
//
//		/* Add all signals. */
//		for(int sig=0; sig<9; sig++)
//		{
//			if(hsig[daynight][sig]->GetNbinsX()<bin) continue;
//			pre += hsig[daynight][sig]->GetBinContent(bin) * nuflux[sig];
//		}
//
//		double dlikelihood = 0;
//		if(mea==0) dlikelihood = pre;
//		else if(mea>=50) dlikelihood = (mea+0.5)*log(mea) - mea + 0.5*log(2*TMath::Pi()) + 1./12./mea - mea*log(pre) + pre;
//		else dlikelihood = log(TMath::Factorial((int)mea)) - mea*log(pre) + pre;
//
//		likelihood += dlikelihood;
//	}
//
//	double chi2 = 2 * likelihood;
//	chi2 += pow(penalty, 2);
//	return chi2;
//}

double JP::chi2(double** bkgflux, double* nuflux, double penalty)
{
	/* chi2 = -2 log likelihood. */
	double likelihood = 0;

	double BinEdge[444];
	for(int i=0; i<=330; i++) BinEdge[i] = 0.1 + 0.005 * i; //0.1~1.75MeV, binwidth = 5keV;
	for(int i=331; i<=415; i++) BinEdge[i] = 1.75 + 0.05 * (i-330); //1.75~6MeV, binwidth = 50keV;
	for(int i=416; i<=435; i++) BinEdge[i] = 6 + 0.2 * (i-415); //6~10MeV, binwidth = 200keV;
	for(int i=436; i<=441; i++) BinEdge[i] = 10 + 0.5 * (i-435); //10~13MeV, binwidth = 500keV;
	BinEdge[442] = 14; BinEdge[443] = 15; BinEdge[444] = 20; //Others.

	/* All predicted and simulated spectra should share the same xmin(0) and same binwidth(5keV). */
	for(int bin=1; bin<=443; bin++) for(int daynight=0; daynight<2; daynight++)
	{
		int lowbin = hsim[daynight]->FindBin(BinEdge[bin-1]);
		int highbin = hsim[daynight]->FindBin(BinEdge[bin])-1;
		double mea = hsim[daynight]->Integral(lowbin, highbin);
		mea /= (1 + penalty*0.015); //Systematic uncertainty: 1.5%
		double pre = 0;

		/* Add all background. */
		for(int bkg=0; bkg<8; bkg++)
		{
			if(hbkg[bkg]->GetXaxis()->GetXmax() < BinEdge[bin-1]) continue; 
			pre += hbkg[bkg]->Integral(lowbin, highbin) * bkgflux[daynight][bkg];
		}

		/* Add all signals. */
		for(int sig=0; sig<9; sig++)
		{
			if(hsig[daynight][sig]->GetXaxis()->GetXmax() < BinEdge[bin-1]) continue;
			pre += hsig[daynight][sig]->Integral(lowbin, highbin) * nuflux[sig];
		}

		double dlikelihood = 0;
		if(mea==0) dlikelihood = pre;
		else if(mea < 1) dlikelihood = 0;
		//else if(mea>=50) dlikelihood = (mea+0.5)*log(mea) - mea + 0.5*log(2*TMath::Pi()) + 1./12./mea - mea*log(pre) + pre;
		//else dlikelihood = log(TMath::Factorial((int)mea)) - mea*log(pre) + pre;
		else dlikelihood = mea * log(mea/pre) + pre - mea;

		likelihood += dlikelihood;
	}

	double chi2 = 2 * likelihood;
	chi2 += pow(penalty, 2);
	return chi2;
}

double JP::chi2le(double* bkgflux, double* nuflux, double penalty)
{
	/* chi2 = -2 log likelihood. */
	double likelihood = 0;

	/* All predicted and simulated spectra should share the same xmin(0) and same binwidth(5keV). */
	for(int bin=hsimtot->FindBin(0.1); bin<=hsimtot->FindBin(1.75); bin++)
	{
		double mea = hsimtot->GetBinContent(bin);
		mea /= (1 + penalty*0.015); //Systematic uncertainty: 1.5%
		double pre = 0;

		/* Add all background. */
		for(int bkg=0; bkg<8; bkg++)
		{
			if(hbkg[bkg]->GetXaxis()->GetXmax()<hsimtot->GetBinCenter(bin)) continue; 
			pre += hbkg[bkg]->GetBinContent(bin) * bkgflux[bkg];
		}

		/* Add all signals. */
		for(int sig=0; sig<9; sig++)
		{
			if(hsigtot[sig]->GetXaxis()->GetXmax() < hsimtot->GetBinCenter(bin)) continue;
			pre += hsigtot[sig]->GetBinContent(bin) * nuflux[sig];
		}

		double dlikelihood = 0;
		if(mea==0) dlikelihood = pre;
		else if(mea < 1) dlikelihood = 0;
		//else if(mea>=50) dlikelihood = (mea+0.5)*log(mea) - mea + 0.5*log(2*TMath::Pi()) + 1./12./mea - mea*log(pre) + pre;
		//else dlikelihood = log(TMath::Factorial((int)mea)) - mea*log(pre) + pre;
		else dlikelihood = mea * log(mea/pre) + pre - mea;

		likelihood += dlikelihood;
	}

	double chi2 = 2 * likelihood;
	chi2 += pow(penalty, 2);
	return chi2;
}

double JP::chi2(double* nuflux)
{
	/* Pass NuFlux to mNuFlux. */
	mNuFlux = nuflux;

	TMinuit* Min = new TMinuit;
	Min->SetPrintLevel(-1);
	int ierflg;
	for(int bkg=0; bkg<8; bkg++)
	{
		Min->mnparm(bkg, (GetBkgName(bkg) + "day").c_str(), GetBkgFlux(bkg), sqrt(GetBkgFlux(bkg)),0,0,ierflg);
		Min->mnparm(bkg+8, (GetBkgName(bkg) + "night").c_str(), GetBkgFlux(bkg), sqrt(GetBkgFlux(bkg)),0,0,ierflg);
	}
	Min->mnparm(16, "penalty", 0, 1, 0, 0, ierflg);
	Min->FixParameter(6);
	Min->FixParameter(14);

	Min->SetFCN(jpfcn);
	Min->SetErrorDef(1);
	Min->SetMaxIterations(50000);
	Min->Migrad();

	double chi2 = 0;
	double edm, errdef;
	int a, b, c;
	Min->mnstat(chi2, edm, errdef, a, b, c);
	delete Min;

	return chi2;
}

double JP::chi2he(double* bkgflux, double* nuflux, double penalty)
{
	/* chi2 = -2 log likelihood. */
	double likelihood = 0;

	/* Rebin */
	//double BinEdge[10] = {1.75, 2.5, 3.5, 5.0, 6.5, 8.0, 9.5, 12.0 ,15.0, 20.0};
	double BinEdge[85];
	for(int i=0; i<80; i++) BinEdge[i] = 1.75+ i * 0.05;
	BinEdge[80] = 7.0;
	BinEdge[81] = 8.5;
	BinEdge[82] = 10.0;
	BinEdge[83] = 12.0;
	BinEdge[84] = 20.0;

	/* All predicted and simulated spectra should share the same xmin(0) and same binwidth(5keV). */
	for(int bin=1; bin<=84; bin++)
	{
		int lowbin = hsimtot->FindBin(BinEdge[bin-1]);
		int highbin = hsimtot->FindBin(BinEdge[bin])-1;
		double mea = hsimtot->Integral(lowbin, highbin);
		mea /= (1 + penalty * 0.015);
		double pre = 0;

		/* Add all background. */
		int Bkg[5] = {2,4,5,6,7};
		for(int bkg=0; bkg<5; bkg++)
		{
			if(hbkg[Bkg[bkg]]->GetXaxis()->GetXmax()<BinEdge[bin-1]) continue; 
			pre += hbkg[Bkg[bkg]]->Integral(lowbin, highbin) * bkgflux[bkg];
		}

		/* Add all signals. */
		int Sig[2] = {5,2};
		for(int sig=0; sig<2; sig++)
		{
			if(hsigtot[Sig[sig]]->GetXaxis()->GetXmax()<BinEdge[bin-1]) continue;
			pre += hsigtot[Sig[sig]]->Integral(lowbin, highbin) * nuflux[sig];
		}

		double dlikelihood = 0;
		if(mea==0) dlikelihood = pre;
		else if(mea<1) dlikelihood = 0;
		else if(mea>=60) dlikelihood = (mea+0.5)*log(mea) - mea + 0.5*log(2*TMath::Pi()) + 1./12./mea - mea*log(pre) + pre;
		else dlikelihood = log(TMath::Factorial((int)mea)) - mea*log(pre) + pre;
		//else dlikelihood = mea * log(mea/pre) + pre - mea;

		likelihood += dlikelihood;
	}

	double chi2 = 2 * likelihood;
	chi2 += pow(penalty, 2);
	return chi2;
}


//double JP::chi2he(double* bkgflux, double* nuflux, double penalty)
//{
//	/* chi2 = -2 log likelihood. */
//	double likelihood = 0;
//
//	/* Rebin */
//	//double BinEdge[10] = {1.75, 2.5, 3.5, 5.0, 6.5, 8.0, 9.5, 12.0 ,15.0, 20.0};
//
//	/* All predicted and simulated spectra should share the same xmin(0) and same binwidth(5keV). */
//	for(int bin=hsimtot->FindBin(EngLow); bin<hsimtot->FindBin(EngUp); bin++)
//	{
//		double mea = hsimtot->GetBinContent(bin);
//		mea /= (1 + penalty * 0.015);
//		double pre = 0;
//
//		/* Add all background. */
//		int Bkg[5] = {2,4,5,6,7};
//		for(int bkg=0; bkg<5; bkg++)
//		{
//			if(hbkg[Bkg[bkg]]->GetNbinsX()<bin) continue; 
//			pre += hbkg[Bkg[bkg]]->GetBinContent(bin) * bkgflux[bkg];
//		}
//
//		/* Add all signals. */
//		int Sig[2] = {5,2};
//		for(int sig=0; sig<2; sig++)
//		{
//			if(hsigtot[Sig[sig]]->GetNbinsX()<bin) continue;
//			pre += hsigtot[Sig[sig]]->GetBinContent(bin) * nuflux[sig];
//		}
//
//		double dlikelihood = 0;
//		if(mea==0) dlikelihood = pre;
//		else if(mea<1) dlikelihood = 0;
//		else if(mea>=60) dlikelihood = (mea+0.5)*log(mea) - mea + 0.5*log(2*TMath::Pi()) + 1./12./mea - mea*log(pre) + pre;
//		else dlikelihood = log(TMath::Factorial((int)mea)) - mea*log(pre) + pre;
//
//		likelihood += dlikelihood;
//	}
//
//	double chi2 = 2 * likelihood;
//	chi2 += pow(penalty, 2);
//	return chi2;
//}


double JP::chi2he(double** bkgflux, double* nuflux, double penalty)
{
	/* chi2 = -2 log likelihood. */
	double likelihood = 0;

	/* All predicted and simulated spectra should share the same xmin(0) and same binwidth(5keV). */
	for(int daynight=0; daynight<2; daynight++) for(int bin=hsim[daynight]->FindBin(EngLow); bin<hsim[daynight]->FindBin(EngUp); bin++)
	{
		double mea = hsim[daynight]->GetBinContent(bin);
		mea /= (1 + penalty * 0.015);
		double pre = 0;

		/* Add all background. */
		int Bkg[5] = {2,4,5,6,7};
		for(int bkg=0; bkg<5; bkg++)
		{
			if(hbkg[Bkg[bkg]]->GetNbinsX()<bin) continue; 
			pre += hbkg[Bkg[bkg]]->GetBinContent(bin) * bkgflux[daynight][bkg];
		}

		/* Add all signals. */
		int Sig[2] = {5,2};
		for(int sig=0; sig<2; sig++)
		{
			if(hsig[daynight][Sig[sig]]->GetNbinsX()<bin) continue;
			pre += hsig[daynight][Sig[sig]]->GetBinContent(bin) * nuflux[sig];
		}

		double dlikelihood = 0;
		if(mea==0) dlikelihood = pre;
		else if(mea>=60) dlikelihood = (mea+0.5)*log(mea) - mea + 0.5*log(2*TMath::Pi()) + 1./12./mea - mea*log(pre) + pre;
		else dlikelihood = log(TMath::Factorial((int)mea)) - mea*log(pre) + pre;

		likelihood += dlikelihood;
	}

	double chi2 = 2 * likelihood;
	chi2 += pow(penalty, 2);
	return chi2;
}

double JP::GetBkgFlux(int bkg)
{
	return bkgflux[bkg];
}

string JP::GetBkgName(int bkg)
{
	return bkgname[bkg];
}

void JP::SetupThreshold(double thres)
{
	threshold = thres;
}

double JP::B8Constrain(double B8flux)
{
	return pow((B8flux-5.25e-4)/0.2e-4, 2);
}

double JP::B8ModelConstrain(double B8flux)
{
	double ModelFlux = solar.GetFlux(mmodel, 6);
	double Error = TMath::Abs(solar.GetFlux(1, 6) - solar.GetFlux(2, 6));
	double chi2 = pow((B8flux - ModelFlux)/Error, 2);
	return chi2;
}


double JP::ModelConstrain(double* nuflux)
{
	double chi2 = 0;
	for(int comp=1; comp<=9; comp++)
	{
		if(comp==6) continue; //B8
		if(comp==9) continue; //F17
		if(comp==5) continue; //Be7_862
		double ModelFlux, Error;
		if(comp==4 || comp==8)
		{
			ModelFlux = solar.GetFlux(mmodel, comp) + solar.GetFlux(mmodel, comp+1);
			Error = TMath::Abs(solar.GetFlux(1, comp)+solar.GetFlux(1, comp+1) - solar.GetFlux(2, comp)-solar.GetFlux(2, comp+1));
			chi2 += pow((nuflux[comp-1] + nuflux[comp] - ModelFlux)/Error, 2);
		}
		else 
		{
			ModelFlux = solar.GetFlux(mmodel, comp);
			Error = TMath::Abs(solar.GetFlux(1, comp) - solar.GetFlux(2, comp));
			chi2 += pow((nuflux[comp-1] - ModelFlux)/Error, 2);
		}
	}
	return chi2;
}

void JP::jpfcn(int &npar, double* gin, double &f, double* par, int iflag)
{
	double* BkgFlux[2];
	double BkgFluxday[8];
	double BkgFluxnight[8];
	for(int bkg=0; bkg<8; bkg++)
	{
		BkgFluxday[bkg] = par[bkg];
		BkgFluxnight[bkg] = par[bkg]+8;
	}
	BkgFlux[0] = BkgFluxday; BkgFlux[1] = BkgFluxnight;

	f = mJP->chi2(BkgFlux, mJP->mNuFlux, par[16]);
}

void JP::SetupNSI(double el, double er, double tl, double tr)
{
	elec->SetupNSI(el,er,tl,tr);
}
