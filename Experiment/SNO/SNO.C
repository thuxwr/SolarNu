#include "SNO.h"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include "TGraph.h"
#include "TMath.h"
#include "Eigen/Dense"

using namespace std;

SNO* SNO::mSNO;

SNO::SNO(int generation)
{
	mSNO = this;
	gene = generation;
	string nu = getenv("neutrino");
	if(nu=="")
	{
		cout << "Environment variable 'neutrino' undefined." << endl; 
		exit(0);
	}

	string filepath = nu + "/Experiment/SNO/data/MC/PredSpec.dat";
	string line;

	nbins = 0;
	Energy = new double[150];
	ifstream fin(filepath.c_str());

	while(getline(fin, line))
	{
		if(line[0]=='#') continue;
		istringstream iss(line);
		iss >> Energy[nbins] >> CC[0][nbins] >> CC[1][nbins] >> ESe[0][nbins] 
												 >> ESe[1][nbins] >> ESo[0][nbins] >> ESo[1][nbins];
		nbins++;
	}

	Energy[nbins] = Energy[nbins-1] + 0.1; //binwidth is 0.1MeV in SNO's MC.

	/* Create Monte Carlo spectrum. */
	for(int daynight=0; daynight<2; daynight++)
	{
		MCspec[daynight] = new TH1D("", "", nbins, Energy);
		for(int bin=1; bin<=nbins; bin++)
		{
			MCspec[daynight]->SetBinContent(bin, CC[daynight][bin-1] + ESe[daynight][bin-1] + ESo[daynight][bin-1]);
		}
	}

	surv = new SurvProb(gene);
	minuit = new TMinuit(5);
	minuit->SetPrintLevel(-1);
}

void SNO::SetupParameter(double sin_2_theta, double ms)
{
	/* Prepare for survival prob and MC distorted spectra. */
	TGraph* probgraph[2];
	TH1D* MCdistspectmp[2];
	TGraph* probsgraph[2]; //For e->s oscillation.

	for(int daynight=0; daynight<2; daynight++)
	{
		/* In order for memory release. */
		MCdistspectmp[daynight] = (TH1D*)mSNO->MCspec[daynight]->Clone("");
		MCdistspec[daynight] = MCdistspectmp[daynight];

		probgraph[daynight] = surv->GetProb(sin_2_theta, ms, 6, daynight, "e", "SNO");
		if(gene==4) probsgraph[daynight] = surv->GetProb(sin_2_theta, ms, 6, daynight, "s", "SNO");

		/* This can be used in sterile neutrino study directly. */
		for(int bin=1; bin<=nbins; bin++)
		{
			double energy = MCdistspec[daynight]->GetBinCenter(bin);
			double probe = probgraph[daynight]->Eval(energy);
			double econt = MCspec[daynight]->GetBinContent(bin);
			if(gene<4)
				MCdistspec[daynight]->SetBinContent(bin, econt * probe);
			else if(gene==4)
			{
				double probs = probsgraph[daynight]->Eval(energy);
				MCdistspec[daynight]->SetBinContent(bin, econt * probe/(1-probs));
			}
		}

		delete probgraph[daynight];
	}

	/* Set parameters according to their best-fit values. */
	minuit->mnparm(0, "c0", 0.317, 0.018, 0, 0, ierflg);
	minuit->mnparm(1, "c1", 0.0039, 0.008, 0, 0, ierflg);
	minuit->mnparm(2, "c2", -0.001, 0.0033, 0, 0, ierflg);
	minuit->mnparm(3, "a0", 0.046, 0.034, 0, 0, ierflg);
	minuit->mnparm(4, "a1", -0.016, 0.027, 0, 0, ierflg);

	minuit->SetFCN(snofcn);
	minuit->SetErrorDef(1); 

	double arglist[10];
	arglist[0] = 1e5; //This is enough.

	minuit->mnexcm("MIGRAD",arglist,0,ierflg);
	
	/* Get fit result. */
	double c0err, c1err, c2err, a0err, a1err;
	minuit->GetParameter(0,c0,c0err);
	minuit->GetParameter(1,c1,c1err);
	minuit->GetParameter(2,c2,c2err);
	minuit->GetParameter(3,a0,a0err);
	minuit->GetParameter(4,a1,a1err);
}

double SNO::chi2()
{
	/* Set correlation matrix and calculate chi square. */
	Eigen::MatrixXd R = Eigen::MatrixXd::Zero(5,5);
	             /*  c0       c1       c2       a0       a1  */
	R << /* c0  */  1.000,  -0.299,  -0.366,  -0.376,   0.129,
			 /* c1  */ -0.299,   1.000,  -0.206,   0.219,  -0.677,
			 /* c2  */ -0.366,  -0.206,   1.000,   0.008,  -0.035,
			 /* a0  */ -0.376,   0.219,   0.008,   1.000,  -0.297,
			 /* a1  */  0.129,  -0.677,  -0.035,  -0.297,   1.000;
	double err[5] = { sqrt(0.016*0.016+0.009*0.009), sqrt(0.0066*0.0066+0.0045*0.0045),
										sqrt(0.0029*0.0029+0.0015*0.0015), sqrt(0.031*0.031+0.014*0.014), sqrt(0.025*0.025+0.011*0.011)};

	/* Covariance matrix. */
	Eigen::MatrixXd Cov = Eigen::MatrixXd::Zero(5,5);
	for(int i=0; i<5; i++) for(int j=0; j<5; j++)
	{
		Cov(i,j) = R(i,j) * err[i] * err[j];
	}

	Eigen::VectorXd vSigEx(5);
	vSigEx << 0.317, 0.0039, -0.0010, 0.046, -0.016;
	Eigen::VectorXd vModel(5);
	vModel << c0, c1, c2, a0, a1;

	double chi2 = 0;
	chi2 = (vSigEx-vModel).transpose() * Cov.inverse() * (vSigEx-vModel);

	/* If minimization failed, set chi2 to be very large. */
	if(ierflg!=0) chi2 = 10000;
	return chi2;
}

double SNO::chi2(double B8flux)
{
	/* Set correlation matrix and calculate chi square. */
	Eigen::MatrixXd R = Eigen::MatrixXd::Zero(6,6);
	             /*  phi      c0       c1       c2       a0       a1  */
	R << /* phi */  1.000,  -0.723,   0.302,  -0.168,   0.028,  -0.012,
			 /* c0  */ -0.723,   1.000,  -0.299,  -0.366,  -0.376,   0.129,
			 /* c1  */  0.302,  -0.299,   1.000,  -0.206,   0.219,  -0.677,
			 /* c2  */ -0.168,  -0.366,  -0.206,   1.000,   0.008,  -0.035,
			 /* a0  */  0.028,  -0.376,   0.219,   0.008,   1.000,  -0.297,
			 /* a1  */ -0.012,   0.129,  -0.677,  -0.035,  -0.297,   1.000;
	double err[6] = { sqrt(0.16*0.16+0.12*0.12), sqrt(0.016*0.016+0.009*0.009), sqrt(0.0066*0.0066+0.0045*0.0045),
										sqrt(0.0029*0.0029+0.0015*0.0015), sqrt(0.031*0.031+0.014*0.014), sqrt(0.025*0.025+0.011*0.011)};

	/* Covariance matrix. */
	Eigen::MatrixXd Cov = Eigen::MatrixXd::Zero(6,6);
	for(int i=0; i<6; i++) for(int j=0; j<6; j++)
	{
		Cov(i,j) = R(i,j) * err[i] * err[j];
	}

	Eigen::VectorXd vSigEx(6);
	vSigEx << 5.25, 0.317, 0.0039, -0.0010, 0.046, -0.016;
	Eigen::VectorXd vModel(6);
	vModel << B8flux*1e-6, c0, c1, c2, a0, a1;

	double chi2 = 0;
	chi2 = (vSigEx-vModel).transpose() * Cov.inverse() * (vSigEx-vModel);

	/* If minimization failed, set chi2 to be very large. */
	if(ierflg!=0) chi2 = 10000;
	return chi2;
}

void SNO::snofcn(int &npar, double* gin, double &f, double* par, int iflag)
{
	double chi2 = 0;
	/* Prepare for predicted event number in each bin. */
	for(int bin=1; bin<=mSNO->nbins; bin++)
	{
		for(int daynight=0; daynight<2; daynight++)
		{
			double pred = 0;
			double theo = 0;
			
			/* Pee_day = c0 + c1(Enu-10MeV) + c2(Enu-10MeV)^2 */
			/* Pee_night = Pee_day * (1+A/2)/(1-A/2)*/
			/* A = a0 + a1(Enu-10MeV) */
			double energy = mSNO->MCdistspec[daynight]->GetBinCenter(bin); 
			double probe = par[0] + par[1]*(energy-10) + par[2]*pow(energy-10,2);
			double Asym = par[3] + par[4]*(energy-10);
			if(daynight==1) probe *= (1+Asym/2.)/(1-Asym/2.);

			/* For sterile neutrinos, see PHYSICAL REVIEW C 88, 025501 (2013), Appendix A. */
			double econt = mSNO->MCspec[daynight]->GetBinContent(bin);

			pred = econt * probe;
			theo = mSNO->MCdistspec[daynight]->GetBinContent(bin);

			/* Get chi2 for each bin. */
			double dchi2 = 0;
			if(pred==0) dchi2 = theo;
			else dchi2 = pred * log(pred/theo) + theo - pred;
			chi2 += dchi2;
		}
	}
	f = chi2 * 2;
}









	




