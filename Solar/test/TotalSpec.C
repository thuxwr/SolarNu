#include "../SolarNu.h"
#include "TH1D.h"

using namespace std;

int main()
{
	SolarNu solar;
	TH1D* TotalSpec = solar.GetTotalSpec(1);
//	TotalSpec->SaveAs("test2.root");
	return 0;
}
