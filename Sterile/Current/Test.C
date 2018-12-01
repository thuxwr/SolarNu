#include "../../Experiment/SNO/SNO.h"
#include "../Sterile.h"
#include <iostream>

using namespace std;

Sterile sterile;
SNO sno(4);

int main()
{
	sno.SetupParameter(1e-5, 1e-5);
	sterile.ncratio(1e-5, 1e-5);
	return 0;
}
