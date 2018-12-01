/*
	 Reproduce SNO livetime for all three phases.

	 Weiran, Apr. 8, 2018.

*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "TH1D.h"
#include "TGraph.h"
#include "TFile.h"
#include "TMath.h"

using namespace std;
int npoints = 1000;

double binwidth = 2./npoints; // -1~1
int livetime()
{
	/* Input all files and graphs. */
	TFile* fD2O = new TFile("./result/livetimeD2O.root","READ");
	TH1D* D2O = (TH1D*)fD2O->Get("livetime");

	TFile* fSalt = new TFile("./result/livetimeSalt.root","READ");
	TH1D* Salt = (TH1D*)fSalt->Get("livetime");

	TFile* fNCD = new TFile("./result/livetimeNCD.root","READ");
	TH1D* NCD = (TH1D*)fNCD->Get("livetime");
	
	/* Change x axis to cos(theta_z), and normalize all spectra. */
	//NCD->Scale(1./NCDorg->Integral("width"));
	double NCDmin = NCD->GetXaxis()->GetXmin();
	double NCDmax = NCD->GetXaxis()->GetXmax();
	TGraph* gNCDorg = new TGraph(NCD);
	TGraph* gNCD = new TGraph(npoints);

	TGraph* gSaltorg = new TGraph(Salt);
	TGraph* gD2Oorg = new TGraph(D2O);

	TGraph* gSalt = new TGraph(npoints);
	TGraph* gD2O = new TGraph(npoints);

	double x = -1 + binwidth/2.;
	int point = 0;

	while(x<1)
	{
		double theta = acos(x); //Zenith angle, in unit rad
		theta = theta * 180. / TMath::Pi(); //Change unit to degree
		theta = 180 - theta; //Change to nadir angle in unit degree

		double y = 0;
		if(theta <= NCDmin) y = 0;
		else if(theta >= NCDmax) y = 0;
		else y = gNCDorg->Eval(theta);

		gNCD->SetPoint(point, x, y);
		gSalt->SetPoint(point, x, gSaltorg->Eval(x));
		gD2O->SetPoint(point, x, gD2Oorg->Eval(x));
		point++;
		x += binwidth;
	}

	/* Normalize according to livetime. D2O:279.27, Salt:389.85, NCD:385.17. */
	double timeD2O, timeSalt, timeNCD;
	timeD2O = gD2O->Integral();
	timeSalt = gSalt->Integral();
	timeNCD = gNCD->Integral();

	TGraph* total = new TGraph(npoints);
	x = -1 + binwidth/2.;
	point = 0;

	while(x<1)
	{
		double xx, yD2O, ySalt, yNCD;
		gD2O->GetPoint(point, xx, yD2O);
		gSalt->GetPoint(point, xx, ySalt);
		gNCD->GetPoint(point, xx, yNCD);

		double y = yD2O/timeD2O * 279.27 + ySalt/timeSalt * 389.85 + yNCD/timeNCD * 385.17;

		total->SetPoint(point, x, y);
		point ++;
		x += binwidth;
	}

	total->SaveAs("./result/livetimeTotal.root");
	delete D2O; delete Salt; delete NCD;

	/* Create a histogram to save the graph and convert x axis to theta and do normalization. */
	/* The output should be the form of JP's livetime. */
	TH1D* SNOlivetime = new TH1D("livetime","SNO livetime for full data set",360,-1,1);
	for(int bin=1; bin<=360; bin++)
	{
		double center = SNOlivetime->GetBinCenter(bin);
		SNOlivetime->SetBinContent(bin, total->Eval(center));
	}
	SNOlivetime->Scale(1./SNOlivetime->Integral("width"));
	SNOlivetime->GetXaxis()->SetTitle("cos(#theta_{z})");
	SNOlivetime->GetYaxis()->SetTitle("livetime(normalized)");

	SNOlivetime->SaveAs("./result/livetime.root");




	return 0;
}

	
