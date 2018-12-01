#include "../SurvProb.h"
#include "TGraph.h"

int main()
{
	SurvProb su(3);
	for(int i=0; i<=100000000; i++)
	{
		TGraph* gr = su.GetProb(0.307, 7e-5, 1, "e");
		delete gr;
	}
	return 0;
}
