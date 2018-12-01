#include "../SNO.h"
#include <iostream>

using namespace std;

int main()
{
	SNO sno(4);
	sno.SetupParameter(1.08e-6, 1.11e-12);
	//for(double flux=3e6; flux<10e6; flux+=1e5) cout << "flux   " << flux << "   chi2  "  << sno.chi2(flux) << endl;
//	cout << sno.chi2(5.3e6) << endl;
//	cout << sno.chi2(5.4e6) << endl;
//	cout << sno.chi2(5.5e6) << endl;
	return 0;
}
