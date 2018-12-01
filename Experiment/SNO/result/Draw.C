/*
	 Draw out SK's contour.

	 Weiran, Mar. 27, 2018.

*/

#include <iostream>
#include <fstream>
#include "TH2D.h"
#include <string>

int nbinstheta = 300;
int nbinsmass = 300;

using namespace std;
int main()
{
	TH2D* contour = new TH2D("SNOcontour", "", nbinstheta, -4, 1, nbinsmass, -12, -3);
	for(int xbin = 1; xbin <= nbinstheta; xbin++)
	{
		string path = "./tmp/";
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
			contour->SetBinContent(binnum, chi2);
		}
		fin.close();
		cout << "Finish " << xbin << "  thetas." << endl;
	}
	double min = contour->GetMinimum();
	contour->GetZaxis()->SetRangeUser(min, min+15);
	contour->SetStats(kFALSE);
	contour->GetXaxis()->SetRangeUser(-2,0);
	contour->GetYaxis()->SetRangeUser(-8,-3.5);
	contour->GetXaxis()->SetTitle("tan^{2}#theta_{12}");
	contour->GetYaxis()->SetTitle("#Delta m^{2}_{12}");
	contour->SetTitle("SNO-only 2#mu analysis");
//	contour->SetContour(4);
//	contour->SetContourLevel(1, min+2.29582);
//	contour->SetContourLevel(2, min+5.99146);
//	contour->SetContourLevel(3, min+11.8290);
	contour->SaveAs("SNO2nu.root");
	return 0;
}
