/*
	 Draw chi2 contour from KamLAND.

	 From http://www.awa.tohoku.ac.jp/KamLAND/datarelease/2ndresult.html

	 Weiran, Mar. 28, 2018.

*/

#include <iostream>
#include <fstream>
#include <sstream>
#include "TH2D.h"
#include <string>

using namespace std;

double tan_2_theta[23000], ms[23000], chi2[23000];
double xbin[1000], ybin[1000];
double Xbin[1000], Ybin[1000];

int Draw()
{
	string filepath = "./data/chi2.dat";
	ifstream fin(filepath.c_str());
	string line;
	int nline=0;
	while(getline(fin, line))
	{
		istringstream iss(line);
		iss >> tan_2_theta[nline] >> ms[nline] >> chi2[nline];
		nline++;
	}

	int nbinx = 1;
	xbin[0] = tan_2_theta[0]; ybin[0] = 100;
	for(int i=1; i<nline; i++)
	{
		if(tan_2_theta[i]!=tan_2_theta[i-1]) 
		{
			xbin[nbinx] = tan_2_theta[i];
			nbinx++;
		}
	}

	int nbiny = 201;
	for(int i=0; i<=nbiny; i++) ybin[i] = ms[i];

	Xbin[0] = xbin[0] - (xbin[1]-xbin[0])/2.;
	for(int i=0; i<nbinx-1; i++)
	{
		Xbin[i+1] = (xbin[i]+xbin[i+1])/2.;
	}
	Xbin[nbinx] = xbin[nbinx-1] + (xbin[nbinx-1]-xbin[nbinx-2])/2.;

	Ybin[0] = ybin[0] - (ybin[1]-ybin[0])/2.;
	for(int i=0; i<nbiny-1; i++)
	{
		Ybin[i+1] = (ybin[i]+ybin[i+1])/2.;
	}
	Ybin[nbiny] = ybin[nbiny-1] + (ybin[nbiny-1]-ybin[nbiny-2])/2.;

	TH2D* th = new TH2D("KamContour","",nbinx,Xbin, nbiny, Ybin);
	for(int i=0; i<nline; i++)
	{
		int binxx = th->GetXaxis()->FindBin(tan_2_theta[i]);
		int binyy = th->GetYaxis()->FindBin(ms[i]);
		th->SetBinContent(binxx, binyy, chi2[i]);
	}

	double min = th->GetMinimum();
	th->GetYaxis()->SetRangeUser(4e-5, 1e-4);
	th->SetContour(3);
	th->SetContourLevel(0,min+2.29582);
	th->SetContourLevel(1,min+5.99146);
	th->SetContourLevel(2,min+11.8290);
	th->SetStats(kFALSE);
	th->GetXaxis()->SetTitle("tan^{2}#theta_{12}");
	th->GetYaxis()->SetTitle("#Delta m^{2}_{12} [eV^{2}]");
	th->GetXaxis()->SetTitleSize(0.045);
	th->GetYaxis()->SetTitleSize(0.045);

	TCanvas* c1 = new TCanvas;

	th->Draw("cont3");
	c1->SetLogx();
	c1->SetLogy();
	th->SaveAs("test.root");



	return 0;
}





