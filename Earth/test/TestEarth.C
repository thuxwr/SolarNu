#include "../Earth.h"
#include <iostream>

using namespace std;

int main()
{
	Earth e;
	double L[12],Ne[12],Nn[12];
	int n;
	double pi = 3.14159265359;
	//e.Intersect(pi,n,L,Ne,Nn);
	//cout << "Zenith angle:  " << pi << endl;
	//cout << n << " segments." << endl;
	//for(int i=0; i<n; i++)
	//{
	//	cout << L[i] << "  " << Ne[i] << "  " << Nn[i] << endl;
	//}
	e.Intersect(pi*0.7,n,L,Ne,Nn);
	cout << "Zenith angle:  " << pi*0.55 << endl;
	cout << n << " segments." << endl;
	for(int i=0; i<n; i++)
	{
		cout << L[i] << "  " << Ne[i] << "  " << Nn[i] << endl;
	}
	return 0;
}
