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
	double xx[nbinstheta+1];
	for(int i=0; i<nbinstheta+1; i++) xx[i] = (0.1+0.8*i/nbinstheta)/(1.1+0.8*i/nbinstheta);

	TH2D* contour = new TH2D("SKcontour", "", nbinstheta, xx, nbinsmass, 1e-5, 15e-5);
	double min = 200;
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
			if(chi2<=min) min = chi2;
			int binnum = contour->GetBin(XBin+1, YBin+1);
			contour->SetBinContent(binnum, chi2);
		}
		fin.close();
		//cout << "Finish " << xbin << "  thetas." << endl;

	}
	contour->SetStats(kFALSE);
	cout << "min  " << min << endl;
	contour->GetZaxis()->SetRangeUser(min, min+15);
	contour->GetXaxis()->SetTitle("sin^{2}#theta_{12}");
	contour->GetYaxis()->SetTitle("#Delta m^{2}_{12}");
	contour->SetTitle("SK I~IV combined contour in LMA");
	//contour->SetContour(4);
	//contour->SetContourLevel(1, min+2.29582);
	//contour->SetContourLevel(2, min+5.99146);
	//contour->SetContourLevel(3, min+11.8290);
	contour->SaveAs("SKLMAnn.root");
	return 0;
}
