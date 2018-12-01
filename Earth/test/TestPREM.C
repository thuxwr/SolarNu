#include "../PREM.h"
#include "TGraph.h"
#include <iostream>

using namespace std;

int main()
{
	PREM prem;
	double cos_zenith[7] = {1,0,-0.2,-0.4,-0.6,-0.8,-1};
	int nseg;
	double L[3000], edens[3000], ndens[3000];

	for(int i=0; i<7; i++)
	{
		cout << "cos zenith:  " << cos_zenith[i] << endl;
		prem.Intersect(cos_zenith[i], nseg, L, edens, ndens);
		cout << "nseg:   " << nseg << endl;
		for(int j=0; j<nseg; j++)
		{
			cout << "L:" << L[j] << "  edens:" << edens[j] << "  ndens:" << ndens[j] << endl;
		}
		cout << endl;
		cout << endl;
		cout << endl;
	}

	return 0;
}


