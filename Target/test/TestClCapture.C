#include "../ClCapture.h"
#include "../../Solar/SolarNu.h"

using namespace std;

int main()
{
	ClCapture capt;
	SolarNu solar;
	double flux[9];
	for(int comp=1; comp<=9; comp++)
	{
		flux[comp-1] = solar.GetFlux(1, comp);
	}
	capt.SetCapture(0.307, 7.5e-5, flux);
	cout << "Cl rate:   " << capt.GetFlux() << endl;
	cout << "Cl systematic error:    " << capt.GetError() << endl;

	return 0;
}
