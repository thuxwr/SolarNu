/*
	 Draw out SK's contour in LMA region.

	 Weiran, Mar. 27, 2018.

*/

#include <iostream>
#include <fstream>
#include "TH2D.h"
#include <string>

const int nbinstheta = 300;
const int nbinsmass = 300;

using namespace std;
int Drawfull()
{
	double xx[nbinstheta+1];
	double yy[nbinsmass+1];
	for(int xbin=0; xbin<=nbinstheta; xbin++)
		xx[xbin] = pow(10,-4+5.* xbin / nbinstheta);
	for(int ybin=0; ybin<=nbinsmass; ybin++)
		yy[ybin] = pow(10, -12 + 9. * ybin / nbinsmass);

	TH2D* contour = new TH2D("JPcontour", "", nbinstheta, xx, nbinsmass, yy);
	double min = 3e4;
	int minx, miny;
	for(int xbin = 1; xbin <= nbinstheta; xbin++)
	{
		string path = "./tmp2/";
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

		if(TMath::Abs(center-a1)>1 && TMath::Abs(center-a2)>1 && TMath::Abs(center-b1)>1 && TMath::Abs(center-b2)>1)
			contour->SetBinContent(xbin, ybin, (a1+a2+b1+b2)/4.);
	}

	for(int xbin=1; xbin<=nbinstheta; xbin++) for(int ybin=1; ybin<=nbinsmass; ybin++)
	{
		double val = contour->GetBinContent(xbin, ybin);
		contour->SetBinContent(xbin, ybin, val-min+0.005);
	}

	cout << "min  " << min << endl;
	cout << "tan^2 theta min:  " << contour->GetXaxis()->GetBinCenter(minx) << endl;
	cout << "ms min:  " << contour->GetYaxis()->GetBinCenter(miny) << endl;
	contour->GetZaxis()->SetRangeUser(0, 100);
	contour->GetXaxis()->SetTitle("tan^{2}#theta_{12}");
	contour->GetYaxis()->SetTitle("#Delta m^{2}_{12}");
	contour->SetTitle("JP contour in full region");
	contour->SetStats(kFALSE);
	//contour->SetContour(4);
	//contour->SetContourLevel(1, min+2.29582);
	//contour->SetContourLevel(2, min+5.99146);
	//contour->SetContourLevel(3, min+11.8290);
	contour->SaveAs("JPfull.root");
	return 0;
}
