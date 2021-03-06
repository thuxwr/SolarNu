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
int main()
{
	TH2D* contour = new TH2D("SNOcontour", "", nbinstheta, 0.1, 0.9, nbinsmass, 1e-5, 15e-5);
	double min = 100;
	int minx, miny;
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
	cout << "min  " << min << endl;
	cout << "tan^2 theta min:  " << contour->GetXaxis()->GetBinCenter(minx) << endl;
	cout << "ms min:  " << contour->GetYaxis()->GetBinCenter(miny) << endl;
	contour->GetZaxis()->SetRangeUser(min, min+15);
	contour->GetXaxis()->SetTitle("tan^{2}#theta_{12}");
	contour->GetYaxis()->SetTitle("#Delta m^{2}_{12}");
	contour->SetTitle("SNO contour in LMA region");
	contour->SetStats(kFALSE);
	contour->SetContour(4);
	contour->SetContourLevel(1, min+2.29582);
	contour->SetContourLevel(2, min+5.99146);
	contour->SetContourLevel(3, min+11.8290);
	contour->SaveAs("SNOLMA2nu.root");
	return 0;
}
