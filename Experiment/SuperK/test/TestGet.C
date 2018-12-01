#include "TH1D.h"
#include "../SK.h"

int main()
{
	SK SuperK;
	TH1D* th = SuperK.GetSpec(1);
	th->SaveAs("th.root");
	return 0;
}
