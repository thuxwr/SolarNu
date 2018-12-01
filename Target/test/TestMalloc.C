#include "../GaCapture.h"

int main()
{
	GaCapture ga;
	double flux[9]={1,1,1,1,1,1,1,1,1};
	for(int i=0; i<100000000; i++)
	{
		ga.SetCapture(0.307, 7e-5,flux);
	}
	return 0;
}
