/*
	 Test survival probability in all regions including LMA, SMA, LOW, Quasi-VAC and VAC.
	 Energy in [0,20]MeV is uniformly divided to 2000 points.

	 Input parameters: argv[1] = core number, argv[2] = energy.
	 Output file: energy \t survprob \n.
	 Weiran, Mar. 2, 2018
*/

#include "../LMA/SurvProb.h"
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

using namespace std;

int main(int argc, char** argv)
{
	double sin_2_t12 = 0.327;
	double ms12 = 4.8e-5;
	double sin_2_t13 = 0.0219;
	double ms13 = 2.5e-3;
	int daynight = 1;

	stringstream ss1;
	ss1 << argv[1];
	stringstream ss2;
	ss2 << argv[2];

	int core;
	ss1 >> core;
	double energy;
	ss2 >> energy;

	char corenum[5]; //At most 99999 cores for one time. I think this is enough.
	sprintf(corenum, "%d", core);
	string Core = corenum;
	string filepath = "./tmp/core" + Core + ".dat"; 

	ofstream fout(filepath.c_str(), ios::app); //Write results to this file.

	fout << left << setw(5) << energy << "\t";

	SurvProb surv(3);
	surv.SetupProb(energy, sin_2_t12, sin_2_t13, ms12, ms13);
	double prob = surv.GetProbFromCalculation(daynight,1,6,0);
	fout << prob << endl;
	fout.close();

	return 0;
}



	
