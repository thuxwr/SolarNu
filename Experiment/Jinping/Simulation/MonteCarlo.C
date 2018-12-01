/*
	 Do Monte-Carlo simulation for 2kt fiducial mass and 1500 days livetime.
	 Background from reactors are not included since rate is too low.

	 Input theory from GS98 and AGS09 respectively.
	 One should notice that Boron 8 neutrino flux is from SNO NC measurement, regardless of model.

	 Weiran, Mar. 16, 2018.

*/

#include <string>
#include <sstream>
#include "TH1D.h"
#include "TFile.h"
#include "TRandom3.h"
#include "../../../Background/Beta/BetaSpec.h"
#include "../../../Background/ExtGamma/ExtGammaSpec.h"
#include "../../../Solar/SolarNu.h"
#include "../../../Target/SolarElec.h"
#include "../../../Detector/FastSim.h"

using namespace std;

struct Expt
{
	char Name[20];
	double Mass; //Target mass after all cuts in 100ton.
	double RunTime; //day

	double C14R; //(day*100ton)
	double Kr85R;
	double Bi210R;
	double C11R;

	double C10R;
	double Tl208R;
	double Be11R;
	double ExtTl208R; // day*100ton*total surface area
};

double day100ton = 24*60*60*3.307e31;
string Daynight[2] = {"day","night"};

void AddHist( TH1D* cur, TH1D* add);

int main(int argc, char** argv)
{
	/* Input oscillation pars. */
	double sin_2_theta, ms;
	double runtime = 1500;
	string outpath;
	if(argc==1) //No input pars, use solar best-fit or global best-fit.
	{
		sin_2_theta = 0.327; ms = 4.8e-5;
		//sin_2_theta = 0.334; ms = 4.8e-5; //SK's best.
		// sin_2_theta = 0.307; ms = 7.5e-5;
	}
	else if(argc==3) //Input mixing angle and mass splitting.
	{
		stringstream ss1;
		ss1 << argv[1];
		ss1 >> sin_2_theta;
		stringstream ss2;
		ss2 << argv[2];
		ss2 >> ms;
	}
	else if(argc==2) //Input a path for output.
	{
		sin_2_theta = 0.327; ms = 4.8e-5;
		stringstream ss;
		ss << argv[1];
		ss >> outpath;
	}
	else if(argc==4) //Input mixing angle, mass splitting and runtime
	{
		stringstream ss1;
		ss1 << argv[1];
		ss1 >> sin_2_theta;
		stringstream ss2; 
		ss2 << argv[2];
		ss2 >> ms;
		stringstream ss3;
		ss3 << argv[3];
		ss3 >> runtime;
	}
		
	else if(argc==5) //Input path, mixing angle, mass splitting and runtime.
	{
		stringstream ss1;
		ss1 << argv[1];
		ss1 >> sin_2_theta;
		stringstream ss2; 
		ss2 << argv[2];
		ss2 >> ms;
		stringstream ss3;
		ss3 << argv[3];
		ss3 >> outpath;
		stringstream ss4;
		ss4 << argv[4];
		ss4 >> runtime;
	}

	else { 
		cout << "Number of input parameters invalid!" << endl;
		exit(0);
	}

	cout << "Generate MC sample with fiducial mass 2kt, runtime " << runtime << " years." << endl;

	int DetectModel = 2; //500PE. 
	Expt Jinping = {"Jinping2k", 20, runtime * 365, 3.46e6, 1, 25.0, 0.15, 0.0031, 0.084, 0.00016, 1.26};
	SolarElec elec(3); //3 generations
	SolarNu solar;

	/* Background. */
	BetaSpec* BetaBkg[7];
	BetaBkg[0] = new BetaSpec("Kr85", 36, 85, -1);
	BetaBkg[1] = new BetaSpec("Bi210",83, 210,-1);
	BetaBkg[2] = new BetaSpec("C11",  6,  11,  1);
	BetaBkg[3] = new BetaSpec("C14",  6,  14, -1);
	BetaBkg[4] = new BetaSpec("C10",  6,  10,  1);
	BetaBkg[5] = new BetaSpec("Tl208",81, 208,-1);
	BetaBkg[6] = new BetaSpec("Be11", 4,  11, -1);

	ExtGammaSpec* ExtGamBkg;
	ExtGamBkg = new ExtGammaSpec("ExtTl208", 0.4);

	FastSim sim;

	TH1D* BackgT[8];
	for(int bkg=0; bkg<7; bkg++) BackgT[bkg] = new TH1D(*(BetaBkg[bkg]->hVisTotal));
	BackgT[7] = new TH1D(*(ExtGamBkg->hVisTotal));

	TH1D* BackgD[8];
	for(int bkg=0; bkg<8; bkg++) 
	{
		sim.Convert(BackgD[bkg], BackgT[bkg], DetectModel);
		BackgD[bkg]->Scale(Jinping.Mass * Jinping.RunTime * 0.5); //Day and night livetime = RunTime * 0.5
		string name = Jinping.Name;
		if(bkg<7) name = name + BetaBkg[bkg]->isotope;
		else name = name + ExtGamBkg->isotope;
		BackgD[bkg]->SetName(name.c_str());
		BackgD[bkg]->SetTitle(name.c_str());
		BackgD[bkg]->GetYaxis()->SetTitle("Events/5keV");
	}

	BackgD[0]->Scale(Jinping.Kr85R);
	BackgD[1]->Scale(Jinping.Bi210R);
	BackgD[2]->Scale(Jinping.C11R);
	BackgD[3]->Scale(Jinping.C14R);
	BackgD[4]->Scale(Jinping.C10R);
	BackgD[5]->Scale(Jinping.Tl208R);
	BackgD[6]->Scale(Jinping.Be11R);
	BackgD[7]->Scale(Jinping.ExtTl208R);

	string name = Jinping.Name;
	double binwidth = 5e-3;
	double emin = 0, emax = 20;
	int nbin = int(floor((emax-emin)/binwidth))+1;
	TH1D* BackgDetSum = new TH1D((name+"BackgDetSum").c_str(),"",nbin, emin, nbin*binwidth + emin);

	for(int bkg=0; bkg<8; bkg++) AddHist(BackgDetSum, BackgD[bkg]);


	/* Signal. For GS and AGS separatedly. */
	for(int model=1; model<=2; model++)
	{
		TH1D* SolarD[2];
		TH1D* PredDet[2];
		TH1D* Sim[2];
		for(int daynight=0; daynight<2; daynight++)
		{
			/* B8 flux for Monte Carlo is from SNO NC measurement, regardless of model. */
			double flux[9];
			for(int comp=1; comp<=9; comp++) flux[comp-1] = solar.GetFlux(model, comp);
			flux[5] = 5.25e-4;
			TH1D* SolarT = elec.GetElecTotal(sin_2_theta, ms, daynight, flux);
			//TH1D* SolarT = elec.GetElecTotal(daynight, model, sin_2_theta, ms);

			sim.Convert(SolarD[daynight], SolarT, DetectModel);

			SolarD[daynight]->SetName((name+"SolarDet"+Daynight[daynight]).c_str());
			SolarD[daynight]->Scale(5e-3);
			SolarD[daynight]->Scale(day100ton * Jinping.Mass * Jinping.RunTime * 0.5); //0.5 * RunTime = day livetime

			PredDet[daynight] = new TH1D((name+"PredDet"+Daynight[daynight]).c_str(),name.c_str(),nbin, emin, nbin*binwidth + emin);
			AddHist(PredDet[daynight], BackgDetSum);
			AddHist(PredDet[daynight], SolarD[daynight]);

			Sim[daynight] = new TH1D((name+"Sim"+Daynight[daynight]).c_str(),name.c_str(),nbin, emin, nbin*binwidth + emin);

			TRandom3 rnd(0); //The seed is computed using UUID.
			for(int bin=1; bin<PredDet[daynight]->GetNbinsX(); bin++)
			{
				double content = PredDet[daynight]->GetBinContent(bin);
				Sim[daynight]->SetBinContent(bin, content);
				//if(content>10) Sim[daynight]->SetBinContent(bin, content);
				//else Sim[daynight]->SetBinContent(bin, rnd.PoissonD(content));
				//if(content<1e9) Sim[daynight]->SetBinContent(bin, rnd.Poisson(content));
				//else Sim[daynight]->SetBinContent(bin, rnd.PoissonD(content));
			}

			delete SolarT;
		}

		string nu = getenv("neutrino");
		//string filename = nu + "/Experiment/Jinping/Simulation/result/" + solar.GetModelName(model) +"-"+ sim.GetModeName(DetectModel) + ".root";
		string filename = nu + "/Experiment/Jinping/Simulation/result/solar/" + solar.GetModelName(model) +"-"+ sim.GetModeName(DetectModel) + ".root";
		if(argc==2 || argc==5) filename = outpath + "/" + solar.GetModelName(model) + "-" + sim.GetModeName(DetectModel) + ".root";
		TFile* file = new TFile(filename.c_str(),"RECREATE");

		for(int bkg=0; bkg<8; bkg++) BackgD[bkg]->Write();
		BackgDetSum->Write();

		for(int daynight=0; daynight<2; daynight++)
		{
			SolarD[daynight]->Write();
			delete SolarD[daynight];
			PredDet[daynight]->Write();
			delete PredDet[daynight];
			Sim[daynight]->Write();
			delete Sim[daynight];
		}


			
	}
	return 0;
}


void AddHist(TH1D* cur, TH1D* add)
{
	for(int bin=1; bin<=add->GetNbinsX(); bin++)
	{
		double addi = add->GetBinContent(bin);
		double curr = cur->GetBinContent(bin);
		cur->SetBinContent(bin, curr+addi);
	}
}


