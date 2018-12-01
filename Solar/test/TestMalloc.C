#include "../SolarNu.h"
#include "TH1D.h"

int main()
{
	SolarNu solar;
//	for(int i=0; i<10000000; i++)
	{
		TH1D* th = solar.GetSpec(1,1,1);
		delete th;
	}
	return 0;
}

