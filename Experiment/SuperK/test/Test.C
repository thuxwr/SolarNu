#include "../SK.h"
#include "TH1D.h"

using namespace std;

int main()
{
	SK superK;
	for(int phase=1; phase<=4; phase++)
	{
		double s = 0;
		TH1D* th = superK.GetSpec(phase);
		for(int bin=1; bin<=superK.GetNbinsX(phase); bin++)
		{


