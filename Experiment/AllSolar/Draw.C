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
	TH2D* contour = new TH2D("AllSolar", "", nbinstheta, -4, 1, nbinsmass, -12, -3);
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
	contour->SaveAs("contour.root");
	return 0;
}
