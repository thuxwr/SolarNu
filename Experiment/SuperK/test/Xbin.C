#include "../SK.h"
#include <iostream>

using namespace std;

int main()
{
	SK SuperK;
	int phase=1;
	int size = SuperK.GetNbinsX(phase);
	for(int i=0; i<size; i++) 
		cout << SuperK.GetStatError(i+1,2,1,0) << endl;
	return 0;
}
