/*
	 Draw sterile.

	 Weiran, Apr. 6, 2018.

*/

#include <iostream>
#include <fstream>
#include "TH2D.h"
#include <string>
#include "TMath.h"

const int nbinstheta = 300;
const int nbinsmass = 300;

using namespace std;
int Draw()
{
	double xx[nbinstheta+1];
	double yy[nbinstheta+1];
	for(int i=0; i<nbinstheta+1; i++) {
		double log_tan_2_theta = 7./nbinstheta * i - 6.;
		double tan_2_theta = pow(10, log_tan_2_theta);
		xx[i] = tan_2_theta;
		double log_ms = -12. + 9. /nbinsmass * i;
		double ms = pow(10, log_ms);
		yy[i] = ms;
	}

	TH2D* contour = new TH2D("RADsterile", "", nbinstheta, xx, nbinsmass, yy);
	double min = 30000;
	for(int xbin = 1; xbin <= nbinstheta; xbin++)
	{
		//if(xbin==275) continue;
		string path = "./tmp/";
		char corenum[5];
		sprintf(corenum, "%d", xbin-1);
		string core = corenum;
		ifstream fin((path+"core"+core+".dat").c_str());
		for(int ybin=1; ybin<=nbinsmass; ybin++)
		{
			int XBin, YBin;
			double chi2 = 30000;
			int binnum = contour->GetBin(xbin, ybin);
			//int binnum = contour->GetBin(XBin+1, YBin+1);
			contour->SetBinContent(binnum, chi2);
			fin >> XBin >> YBin >> chi2;
			if(chi2<=min) min = chi2;
			contour->SetBinContent(binnum, chi2);
		}
		fin.close();
		cout << "Finish " << xbin << "  thetas." << endl;

	}

	for(int xbin=2; xbin<nbinstheta; xbin++) for(int ybin = 2; ybin<nbinsmass; ybin++)
	{
		double center = contour->GetBinContent(xbin, ybin);
		double a1 = contour->GetBinContent(xbin-1, ybin);
		double a2 = contour->GetBinContent(xbin+1, ybin);
		double b1 = contour->GetBinContent(xbin, ybin-1);
		double b2 = contour->GetBinContent(xbin, ybin+1);

		//if(TMath::Abs(center-a1)>0.3 && TMath::Abs(center-a2)>0.3 && TMath::Abs(center-b1)>0.3 && TMath::Abs(center-b2)>0.3)
		//	contour->SetBinContent(xbin, ybin, (a1+a2+b1+b2)/4.);
	}

	cout << "min  " << min << endl;
	contour->GetZaxis()->SetRangeUser(min, min+15);
	contour->GetXaxis()->SetTitle("sin^{2}#theta_{12}");
	contour->GetYaxis()->SetTitle("#Delta m^{2}_{12}");
	contour->SetTitle("All solar combined contour in LMA");
	contour->SetStats(kFALSE);
	//contour->SetContour(4);
	//contour->SetContourLevel(1, min+2.29582);
	//contour->SetContourLevel(2, min+5.99146);
	//contour->SetContourLevel(3, min+11.8290);
	contour->SaveAs("test.root");
	return 0;
}
