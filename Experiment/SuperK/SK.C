#include <iostream>
#include "SK.h"
#include <string>
#include <fstream>
#include <sstream>
#include "TMath.h"

using namespace std;

SK::SK()
{
	string nu = getenv("neutrino");
	if(nu=="")
	{
		cout << "Environment variable 'neutrino' undefined." << endl;
		exit(0);
	}

	SKpath = nu + "/Experiment/SuperK";

	PhaseName[0] = "SK1";
	PhaseName[1] = "SK2";
	PhaseName[2] = "SK3";
	PhaseName[3] = "SK4";

	for(int phase=1; phase<=4; phase++)
	{
		string line;

		/* Get correlated systematic error. */
		ifstream CorreErrorFile((SKpath+"/data/Systematic/"+PhaseName[phase-1]+"-SysErr-Corre.dat").c_str());
		nbinsCorre[phase-1] = 0;
		while(getline(CorreErrorFile, line))
		{
			if(line[0]=='#') continue;
			istringstream iss(line);
			iss >> CorreEnergyLowEdge[phase-1][nbinsCorre[phase-1]] 
					>> CorreEnergyHighEdge[phase-1][nbinsCorre[phase-1]]
				  >> SpecErrorHigh[phase-1][nbinsCorre[phase-1]]      
					>> SpecErrorLow[phase-1][nbinsCorre[phase-1]]
					>> ScaleErrorHigh[phase-1][nbinsCorre[phase-1]]     
					>> ScaleErrorLow[phase-1][nbinsCorre[phase-1]]
					>> ResoErrorHigh[phase-1][nbinsCorre[phase-1]]      
					>> ResoErrorLow[phase-1][nbinsCorre[phase-1]];
			nbinsCorre[phase-1]++;
		}

		/* Get uncorrelated systematic error. */
		ifstream UncorreErrorFile((SKpath+"/data/Systematic/"+PhaseName[phase-1]+"-SysErr-Uncorre.dat").c_str());
		nbinsUncorre[phase-1] = 0;
		while(getline(UncorreErrorFile, line))
		{
			if(line[0]=='#') continue;
			istringstream iss(line);
			iss >> UncorreEnergyLowEdge[phase-1][nbinsUncorre[phase-1]]
				  >> UncorreEnergyHighEdge[phase-1][nbinsUncorre[phase-1]]
					>> UncorreErrorHigh[phase-1][nbinsUncorre[phase-1]]
					>> UncorreErrorLow[phase-1][nbinsUncorre[phase-1]];
			nbinsUncorre[phase-1]++;
		}

		/* Get detected rate. */
		ifstream RateFile((SKpath+"/data/Rate/"+PhaseName[phase-1]+"-rate.dat").c_str());
		nbinsRate[phase-1] = 0;
		while(getline(RateFile, line))
		{
			if(line[0]=='#') continue;
			istringstream iss(line);
			iss >> RateEnergyLowEdge[phase-1][nbinsRate[phase-1]]
					>> RateEnergyHighEdge[phase-1][nbinsRate[phase-1]]
					>> ObservedRateAll[phase-1][nbinsRate[phase-1]]
					>> RateAllErrorHigh[phase-1][nbinsRate[phase-1]]
					>> RateAllErrorLow[phase-1][nbinsRate[phase-1]]
					>> ObservedRateDay[phase-1][nbinsRate[phase-1]]
					>> RateDayErrorHigh[phase-1][nbinsRate[phase-1]]
					>> RateDayErrorLow[phase-1][nbinsRate[phase-1]]
					>> ObservedRateNight[phase-1][nbinsRate[phase-1]]
					>> RateNightErrorHigh[phase-1][nbinsRate[phase-1]]
					>> RateNightErrorLow[phase-1][nbinsRate[phase-1]]
					>> ExpectedB8Rate[phase-1][nbinsRate[phase-1]]
					>> ExpectedhepRate[phase-1][nbinsRate[phase-1]];
			nbinsRate[phase-1]++;
		}

		/* Get zenith spectra. */


	}
}

double SK::GetCorreSysError(double energy, int source, int phase, int sign)
{
	double error = 0;
	double Energy; //Convert to total energy for SK1~SK3.
	if(phase==4) Energy = energy + 0.010;
	else Energy = energy + 0.510;
	for(int bin=0; bin<nbinsCorre[phase-1]; bin++)
	{
		if(CorreEnergyHighEdge[phase-1][bin]<=Energy) continue;
		switch(source)
		{
			case 1:
				if(sign==0) error = SpecErrorHigh[phase-1][bin];
				if(sign==1) error = SpecErrorLow[phase-1][bin];
				break;
			case 2:
				if(sign==0) error = ScaleErrorHigh[phase-1][bin];
				if(sign==1) error = ScaleErrorLow[phase-1][bin];
				break;
			case 3:
				if(sign==0) error = ResoErrorHigh[phase-1][bin];
				if(sign==1) error = ResoErrorLow[phase-1][bin];
				break;
			default:
				cout << "Incorrect source for correlated systematic error." << endl;
				break;
		}
		break;
	}
	return error * 0.01;
}

double SK::GetUncorreSysError(double energy, int phase, int sign)
{
	double error = 0;
	double Energy; //Convert to total energy for SK1~SK3.
	if(phase==4) Energy = energy + 0.010;
	else Energy = energy + 0.510;
	for(int bin=0; bin<nbinsUncorre[phase-1]; bin++)
	{
		if(UncorreEnergyHighEdge[phase-1][bin]<=Energy) continue;
		if(sign==0) error = UncorreErrorHigh[phase-1][bin];
		if(sign==1) error = UncorreErrorLow[phase-1][bin];
		break;
	}
	return error;
}

double SK::GetRate(int bin, int daynight, int phase)
{
	double rate = 0;
	switch(daynight)
	{
		case 0: rate = ObservedRateDay[phase-1][bin-1]; break;
		case 1: rate = ObservedRateNight[phase-1][bin-1]; break;
		case 2: rate = ObservedRateAll[phase-1][bin-1]; break;
		default: cout << "Incorrect daynight option." << endl; break;
	}
	return rate;
}

double SK::GetStatError(int bin, int daynight, int phase, int sign)
{
	double error = 0;
	switch(daynight)
	{
		case 0:
			if(sign==0) error = RateDayErrorHigh[phase-1][bin-1];
			if(sign==1) error = RateDayErrorLow[phase-1][bin-1];
			break;
		case 1:
			if(sign==0) error = RateNightErrorHigh[phase-1][bin-1];
			if(sign==1) error = RateNightErrorLow[phase-1][bin-1];
			break;
		case 2:
			if(sign==0) error = RateAllErrorHigh[phase-1][bin-1];
			if(sign==1) error = RateAllErrorLow[phase-1][bin-1];
			break;
		default:
			cout << "Incorrect daynight option." << endl;
			break;
	}
	return error;
}

int SK::GetNbinsX(int phase)
{
	return nbinsRate[phase-1];
}

double* SK::GetXbins(int phase)
{
	double* xbin = new double[nbinsRate[phase-1]+1];
	for(int bin=0; bin<nbinsRate[phase-1]; bin++)
	{
		if(phase==4)
			xbin[bin] = RateEnergyLowEdge[phase-1][bin] - 0.010;
		else
			xbin[bin] = RateEnergyLowEdge[phase-1][bin] - 0.510; //Convert total energy to kinetic energy.
	}
	xbin[nbinsRate[phase-1]] = RateEnergyHighEdge[phase-1][nbinsRate[phase-1]-1];
	return xbin;
}

double SK::GetFluxSysError(int phase) //From SK-IV, Table V, subtotal.
{
	double error = 0;
	switch(phase)
	{
		case 1: error = 2.8; break;
		case 2: error = 4.8; break;
		case 3: error = 1.6; break;
		case 4: error = 1.2; break;
		default: break;
	}
	error *= 0.01;
	return error;
}

double SK::GetUncorreSysError(double energy, int phase)
{
	double up = GetUncorreSysError(energy, phase, 0);
	double low = GetUncorreSysError(energy, phase, 1);
	double avg = (up-low)/2.;
	return avg;
}

double SK::GetStatError(int bin, int daynight, int phase)
{
	double up = GetStatError(bin, daynight, phase, 0);
	double low = GetStatError(bin, daynight, phase, 1);
	double avg = (up-low)/2.;
	return avg;
}

double SK::GetRateMCB8(int bin, int phase)
{
	double rate = ExpectedB8Rate[phase-1][bin-1];
	return rate;
}

double SK::GetRateMChep(int bin, int phase)
{
	double rate = ExpectedhepRate[phase-1][bin-1];
	return rate;
}

double SK::GetRateMC(int bin, int phase)
{
	double rateB8 = GetRateMCB8(bin, phase);
	double ratehep = GetRateMChep(bin, phase);
	double rate = rateB8 + ratehep;
	return rate;
}

TH1D* SK::GetSpec(int phase, int daynight)
{
	int nbins = GetNbinsX(phase);
	TH1D* spec = new TH1D("", "", nbins, GetXbins(phase));
	for(int bin=1; bin<=nbins; bin++)
	{
		double kinenergy = spec->GetBinCenter(bin);
		double rate = GetRate(bin, daynight, phase);
		double mcrate;
		if(phase==4) mcrate = GetRateMC(bin, phase);
		else mcrate = GetRateMCB8(bin, phase)*5.25/5.79  + GetRateMChep(bin,phase);
		double ratio = rate/mcrate;
		/* Error here for night rate=0. */
		if(rate==0)
		{
			spec->SetBinContent(bin, 0);
			spec->SetBinError(bin, GetStatError(bin, daynight, phase)/mcrate);
			continue;
		}
		double error1 = GetStatError(bin, daynight, phase)/rate;
		double error2 = GetUncorreSysError(kinenergy, phase) * 0.01;
		double error = sqrt(error1*error1 + error2*error2);
		spec->SetBinContent(bin, ratio);
		spec->SetBinError(bin, error*ratio);
	}
	return spec;
}

