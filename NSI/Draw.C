/*
	 Draw out SK's contour in LMA region.

	 Weiran, Mar. 27, 2018.

*/

#include <iostream>
#include <fstream>
#include "TH2D.h"
#include <string>

int nbinstheta = 300;
int nbinsmass = 300;

using namespace std;
int Draw()
{
	TH2D* contour = new TH2D("JPcontour", "", nbinstheta, -1.5, 1.5, nbinsmass, -1.5, 1.5);
	double min = 10000;
	int minx, miny;
	for(int xbin = 1; xbin <= nbinstheta; xbin++)
	{
		string path = "./result/tmp2/";
		char corenum[5];
		sprintf(corenum, "%d", xbin-1);
		string core = corenum;
		ifstream fin((path+"core"+core+".dat").c_str());
		for(int ybin=1; ybin<=nbinsmass; ybin++)
		{
			int XBin, YBin;
			double chi2;
			fin >> XBin >> YBin >> chi2;
			int binnum = contour->GetBin(XBin+1, YBin+1);
			if(chi2<=min) {
				min = chi2;
				minx = XBin+1;
				miny = YBin+1;
			}
			contour->SetBinContent(binnum, chi2);
		}
		fin.close();
		//cout << "Finish " << xbin << "  thetas." << endl;

	}
	/* Check if any bin failed. */
	for(int xbin = 2; xbin < nbinstheta; xbin++) for(int ybin=2; ybin<nbinsmass; ybin++)
	{
		double center = contour->GetBinContent(xbin, ybin);
		double a1 = contour->GetBinContent(xbin-1, ybin);
		double a2 = contour->GetBinContent(xbin+1, ybin);
		double b1 = contour->GetBinContent(xbin, ybin-1);
		double b2 = contour->GetBinContent(xbin, ybin+1);

		//if(TMath::Abs(center-a1)>0.2 && TMath::Abs(center-a2)>0.2 && TMath::Abs(center-b1)>0.2 && TMath::Abs(center-b2)>0.2)
		//	contour->SetBinContent(xbin, ybin, (a1+a2+b1+b2)/4.);
		if(TMath::Abs(center-a1)>2 || TMath::Abs(center-a2)>2 || TMath::Abs(center-b1)>2 || TMath::Abs(center-b2)>2)
		{
			double mini;
			mini = a1;
			if(a2<a1) mini = a2;
			if(b1<mini) mini = b1;
			if(b2<mini) mini = b2;
		//	if(TMath::Abs(mini-center)>1)
		//		contour->SetBinContent(xbin, ybin, mini);
		}
	}

	cout << "min  " << min << endl;
	cout << "tan^2 theta min:  " << contour->GetXaxis()->GetBinCenter(minx) << endl;
	cout << "ms min:  " << contour->GetYaxis()->GetBinCenter(miny) << endl;
	contour->GetZaxis()->SetRangeUser(min, min+15);
	contour->GetXaxis()->SetTitle("ed");
	contour->GetYaxis()->SetTitle("en");
	contour->SetTitle("JP contour in LMA region");
	contour->SetStats(kFALSE);
	contour->SetContour(4);
	contour->SetContourLevel(1, min+2.29582);
	contour->SetContourLevel(2, min+5.99146);
	contour->SetContourLevel(3, min+11.8290);
	//TCanvas* c1 = new TCanvas;
	contour->SaveAs("test.root");
	//contour->Draw("colz");
	//c1->SaveAs("test.C");
	return 0;
}
