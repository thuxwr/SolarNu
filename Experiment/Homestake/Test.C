#include "ClConstrain.h"
#include "../../Solar/SolarNu.h"
#include <iostream>

using namespace std;

int main()
{
	ClConstrain cl;
	SolarNu solar;
	cl.SetupParameter(0.307, 7.5e-5);
	double flux[9];
	for(int comp=1; comp<=9; comp++)
	{
		flux[comp-1] = solar.GetFlux(1,comp);
	}
	flux[5] *= 2;
	cout << "B8 flux:  " << flux[5] << endl;
	cout << "chi2:   " << cl.chi2(flux) << endl;
	return 0;
}

