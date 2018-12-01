#include <iostream>
#include "SurvProb.h"
#include "TFile.h"
#include "TMath.h"
#include "TH1D.h"

using namespace std;

SurvProb::SurvProb(int generation)
{
	gene = generation;
	string nu = getenv("neutrino");
	if(nu=="")
	{
		cout << "Environment variable 'neutrino' undefined." << endl;
		exit(0);
	}

	Daynight[0] = "day";
	Daynight[1] = "night";
	path = nu + "/SurvProb";
	if(gene==2) path = path + "/2nu";
	if(gene==3) path = path + "/3nu";
	if(gene==4) path = path + "/4nu/solar";
	Path = path;
	LMApath = nu + "/SurvProb/LMA";
	LMA2nupath = nu + "/SurvProb/LMA2nu";
	IHpath = nu + "/SurvProb/LMA_inv";
}

TGraph* SurvProb::GetProb(double sin_2_theta, double ms, int comp, int daynight, string flavorname, string experiment, int whichtheta)
{
	path = Path;
	if(gene==3)
	{
		if(whichtheta==0) path = Path + "/12";
		else path = Path + "/13";
	}
	/* Extrapolate value at any point, using linear interpolation. */
	double thetabin = GetThetaBin(sin_2_theta);
	double massbin = GetMassBin(ms);

	bool IsLMA = false;
	/* For LMA region and 3 generations, use refined calculation. */
	if(gene==3 && sin_2_theta<=0.898/1.898 && sin_2_theta>=0.102/1.102 && ms<=14.965e-5 && ms>=1.035e-5)
	{
		IsLMA = true;
		path = LMApath;
		thetabin = GetThetaBinLMA(sin_2_theta);
		massbin = GetMassBinLMA(ms);
	}

	if(gene==2 && sin_2_theta<=0.898/1.898 && sin_2_theta>=0.102/1.102 && ms<=14.965e-5 && ms>=1.035e-5)
	{
		IsLMA = true;
		path = LMA2nupath;
		thetabin = GetThetaBinLMA(sin_2_theta);
		massbin = GetMassBinLMA(ms);
	}

	if(gene==3 && sin_2_theta<=0.898/1.898 && sin_2_theta>=0.102/1.102 && ms>=-14.965e-5 && ms<=-1.035e-5) //IH
	{
		IsLMA = true;
		path = IHpath;
		thetabin = GetThetaBinLMA(sin_2_theta);
		massbin = GetMassBinLMA(0.-ms);
	}

	char theta1[5], theta2[5];
	sprintf(theta1,"%d",(int)thetabin);
	sprintf(theta2,"%d",(int)(thetabin+1));
	char mass1[5], mass2[5];
	sprintf(mass1,"%d",(int)massbin);
	sprintf(mass2,"%d",(int)(massbin+1));

	TFile* f1 = new TFile((path+"/theta"+theta1+"/mass"+mass1+"/SurvProb.root").c_str(),"READ");
	TFile* f2 = new TFile((path+"/theta"+theta1+"/mass"+mass2+"/SurvProb.root").c_str(),"READ");
	TFile* f3 = new TFile((path+"/theta"+theta2+"/mass"+mass1+"/SurvProb.root").c_str(),"READ");
	TFile* f4 = new TFile((path+"/theta"+theta2+"/mass"+mass2+"/SurvProb.root").c_str(),"READ");

	string graphname = solar.GetCompName(comp) + Daynight[daynight] + flavorname;
	if(IsLMA && daynight==1 && gene!=4) graphname = solar.GetCompName(comp) + Daynight[daynight] + experiment + flavorname;

	TGraph* g1 = (TGraph*)f1->Get(graphname.c_str());
	TGraph* g2 = (TGraph*)f2->Get(graphname.c_str());
	TGraph* g3 = (TGraph*)f3->Get(graphname.c_str());
	TGraph* g4 = (TGraph*)f4->Get(graphname.c_str());

	/* linear interpolation */
	double a = thetabin - (int)thetabin;
	double b = 1. - a;
	double c = massbin - (int)massbin;
	double d = 1. - c;

	int Nbins = 200;
	if(g1->GetN()>Nbins) Nbins = g1->GetN();
	if(g2->GetN()>Nbins) Nbins = g2->GetN();
	if(g3->GetN()>Nbins) Nbins = g3->GetN();
	if(g4->GetN()>Nbins) Nbins = g4->GetN();

	double eng[Nbins], prob[Nbins];
	for(int bin=0; bin<Nbins; bin++)
	{
		eng[bin] = 0.05 + 20. * bin / Nbins;
		prob[bin] = g1->Eval(eng[bin])*b*d + g2->Eval(eng[bin])*b*c + g3->Eval(eng[bin])*a*d + g4->Eval(eng[bin])*a*c;
	}

	TGraph* ProbGraph = new TGraph(Nbins, eng, prob);
	ProbGraph->GetXaxis()->SetTitle(g1->GetXaxis()->GetTitle());
	ProbGraph->GetYaxis()->SetTitle(g1->GetYaxis()->GetTitle());
	ProbGraph->SetTitle(g1->GetTitle());

	f1->Close();
	f2->Close();
	f3->Close();
	f4->Close();
	delete g1; delete g2; delete g3; delete g4;
	delete f1; delete f2; delete f3; delete f4;
	return ProbGraph;
}

double SurvProb::GetProbAnaly(double Energy, double sin_2_t12, double ms12, double sin_2_t13, double ms13, int comp)
{
	TH1D* Dist = solar.GetBinnedFluxDist(1, comp);
	int nDistBins = Dist->GetNbinsX();
	double prob = 0;
	for(int bin=1; bin<=nDistBins; bin++)
	{
		double dflux = Dist->GetBinContent(bin) * Dist->GetBinWidth(bin);
		double radius = Dist->GetBinCenter(bin);
		double ne = solar.GetEDensity(radius);
		prob += dflux * Prob(ne, Energy, sin_2_t12, ms12, sin_2_t13, ms13);
	}
	return prob;
}

TGraph* SurvProb::GetProbAnaly(double sin_2_t12, double ms12, double sin_2_t13, double ms13, int comp)
{
	int n=200;
	double x[n], y[n];
	for(int i=0; i<n; i++)
	{
		x[i] = 0.05 + i * 20. / n;
		y[i] = GetProbAnaly(x[i], sin_2_t12, ms12, sin_2_t13, ms13, comp);
	}

	TGraph* graph = new TGraph(n, x, y);
	graph->SetTitle("Survival Probability");
	graph->GetXaxis()->SetTitle("Energy [MeV]");
	graph->GetYaxis()->SetTitle("Probability");

	return graph;
}
	




double SurvProb::GetThetaBin(double sin_2_theta)
{
	double tan_2_theta = sin_2_theta/(1-sin_2_theta);
	double log_theta = log(tan_2_theta)/log(10.);
	double thetabin = 0; 
	if(gene<=3) thetabin = (log_theta + 4) * 199. / 5.;
	if(gene==4) thetabin = (log_theta + 6) * 199. / 7.;
	return thetabin;
}

double SurvProb::GetMassBin(double ms)
{
	double log_ms = log(ms)/log(10.);
	double msbin = 0;
	if(gene<=3) msbin = (log_ms + 12) * 199. / 9.;
	if(gene==4) msbin = (log_ms + 12) * 199. / 9.;
	return msbin;
}

double SurvProb::Prob(double ne, double Energy, double sin_2_t12, double ms12, double sin_2_t13, double ms13)
{
	double cos_2t12 = 1 - 2 * sin_2_t12;
	double sin_2_2t12 = 1 - pow(cos_2t12, 2);

	double Gamma_r = ms12 * sin_2_2t12 * 7.7978e8 /(2. * Energy * cos_2t12 );
	double F = 1 - sin_2_t12/(1-sin_2_t12);
	double p = -TMath::Pi()/2. * Gamma_r * F;
	double Pc = (exp(p)-exp(p/sin_2_t12))/(1-exp(p/sin_2_t12));

	double cos_2_t13 = 1 - sin_2_t13;
	double cos_4_t13 = pow(cos_2_t13, 2);
	double sin_4_t13 = pow(sin_2_t13, 2);
	double beta = 2.5348e-31 * cos_2_t13 * (ne * cos_2_t13) * Energy / ms12;
	double cos_2tm12 = (cos_2t12 - beta )/ sqrt(pow((cos_2t12 - beta), 2) + sin_2_2t12);

	double prob = 0;
	prob = cos_4_t13 * (1./2. + (1./2. - Pc ) * cos_2tm12 * cos_2t12) + sin_4_t13;

	return prob;
}

TGraph* SurvProb::GetProb(double sin_2_theta, double ms, int comp, string flavorname, string experiment, int whichtheta)
{
	TGraph* dayprob = GetProb(sin_2_theta, ms, comp, 0, flavorname, experiment, whichtheta);
	TGraph* nightprob = GetProb(sin_2_theta, ms, comp, 1, flavorname, experiment,  whichtheta);

	int npoints = dayprob->GetN();
	TGraph* avgprob = new TGraph(npoints);
	for(int point = 0; point<npoints; point++)
	{
		double x1, y1, x2, y2;
		dayprob->GetPoint(point, x1, y1);
		nightprob->GetPoint(point, x2, y2);
		if(x1!=x2) cout << "An error occured while averaging daytime and nighttime survival probability." << endl;
		avgprob->SetPoint(point, x1, (y1+y2)/2.);
	}
	delete dayprob; delete nightprob;
	return avgprob;
}

double SurvProb::GetThetaBinLMA(double sin_2_theta)
{
	double tan_2_theta = sin_2_theta/(1-sin_2_theta);
	double thetabin = (tan_2_theta-0.1) * 200. / 0.8 - 0.5;
	return thetabin;
}

double SurvProb::GetMassBinLMA(double ms)
{
	double massbin = (ms-1e-5) * 200. / 14e-5 - 0.5;
	return massbin;
}

