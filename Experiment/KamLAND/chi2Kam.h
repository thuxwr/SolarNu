/*
	 Get KamLAND chi2 for given mass and theta.

	 Weiran, Mar. 28, 2018.

*/
#include "TH2D.h"
#include "TFile.h"
#include <string>
#include <iostream>

using namespace std;

class chi2Kam
{
	public:
		chi2Kam()
		{
			string nu = getenv("neutrino");
			if(nu=="")
			{
				cout << "Environment variable 'neutrino' undefined." << endl;
				exit(0);
			}

			string filepath = nu + "/Experiment/KamLAND/KamLAND.root";
			file = new TFile(filepath.c_str(),"READ");
			KamLAND = (TH2D*)file->Get("KamContour");
		}
		~chi2Kam(){};

		double chi2(double tan_2_theta, double ms)
		{
			return KamLAND->Interpolate(tan_2_theta, ms);
		}

	private:
		TFile* file;
		TH2D* KamLAND;
};
