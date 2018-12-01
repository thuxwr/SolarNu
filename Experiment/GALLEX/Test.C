#include "GaConstrain.h"
#include "../../Solar/SolarNu.h"
#include <iostream>

using namespace std;

int main()
{
	GaConstrain ga;
	SolarNu solar;
	ga.SetupParameter(0.307, 7.5e-5);
	double flux[9];
	for(int comp=1; comp<=9; comp++)
	{
		flux[comp-1] = solar.GetFlux(1,comp);
	}
	flux[4] *= 2;
	cout << "Be7 flux:  " << flux[4] << endl;
	cout << "chi2:   " << ga.chi2(flux) << endl;
	return 0;
}

