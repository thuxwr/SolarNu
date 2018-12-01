/*
	 Calculate livetime for CJPL, assuming data-taking time is an integer in unit year. True livetime distribution should be re-calculated after all data being taken.
	 The latitude of Jinping is 28.2.

	 Weiran, Feb 25, 2018

	 Change x axis to cos(theta_z) in order to make JP consistent with livetime in SNO.
	 Apr. 8, 2018.

*/

#include <iostream>

using namespace std;

void livetime()
{	
	double zen = 0; //Solar zenith angle
//	double lat = 28.2; //local latitude of CJPL
	double lat = 36.42; //latitude of SK
//	double lat = 46.5; //latitude of SNO
	lat = lat * 3.14159265359 / 180.; //change Unit from deg to rad.
	double dec = 0; //declination of the Sun
	double h = 0; //hour angle
	int N = 0; //number of days, calculated from Jan 1.

	TH1D* livetime = new TH1D("livetime","livetime",360,-1,1);

	/* Average */
	for(N=0;N<365;N++)
	{
		/*For approximation, here I calculate the average of the whole year. 
			It should be changed due to the true experiment days.*/
		//dec = -asin(0.39779*cos(0.0172028*(N+10)+0.0334056*sin(0.0172028*(N-2))));
		dec = -23.44 * 3.14159265359 / 180. * cos(360./365. * 3.14159265359/180. * (N+10));
		
		for(int sec=0; sec<24*60*60; sec++) //One day
		{
			h = (sec/3600. -12)*15*3.14159265359 / 180.;
			double cos_zen = sin(lat)*sin(dec)+cos(lat)*cos(dec)*cos(h);
			zen = acos(cos_zen);
			livetime->Fill(cos_zen);
		}
	}

	double integ = livetime->Integral("width");
	livetime->Scale(1/integ);
	livetime->GetXaxis()->SetTitle("cos(#theta_{z})");
	livetime->GetYaxis()->SetTitle("livetime");
	livetime->SetTitle("SuperK livetime distribution");
	livetime->SaveAs("../../SuperK/data/Zenith/livetime.root");
}

		
