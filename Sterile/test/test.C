#include "../Sterile.h"
#include <iostream>

using namespace std;

int main()
{
	Sterile sterile;
	double ratio = sterile.ncratio(1e-3, 1e-5);
	cout << ratio << endl;
	return 0;
}

